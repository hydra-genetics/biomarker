import csv
import gzip
import os
import shutil
import tempfile
import unittest


class TestUnitUtils(unittest.TestCase):
    """Hand-computed checks of prs_daniel.py on the fixtures in .tests/unit.

    Every marker in the score file exercises one feature of the script:
        1:100        heterozygous SNP                       1 x  0.5
        1:200        phased homozygous SNP                  2 x -0.2
        1:250        inside a hom-ref gVCF block            0 x  0.7 (called, not missing)
        1:300        absent from the gVCF                   NA
        1:400        no-call '.|1' inside a hom-ref block   NA (2 x 0.9 with --r-compat)
        2:500        indel, CAA/CAAA in the gVCF vs C/CA    1 x  1.0
        3:100+3:200  multi-allelic pair, both 1/1           capped at 2 x 1.0, one row kept
        4:100+4:200  proxy missing, replaced site 1/1       2 x  2.0 (proxy weight), one row kept
    """

    def setUp(self):
        self.score_file = ".tests/unit/bcprs_score_file.csv"
        self.vcf = ".tests/unit/bcprs_HD832.HES45_T.g.vcf"
        self.tempdir = tempfile.mkdtemp()
        self.detailed = os.path.join(self.tempdir, "HD832.HES45_T.detailed_scores.csv")
        self.total = os.path.join(self.tempdir, "HD832.HES45_T.total_scores.csv")

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def _run(self, vcf=None, extra=()):
        from prs_daniel import main

        rc = main(["-s", self.score_file, "-v", vcf or self.vcf, "-d", self.detailed, "-t", self.total] + list(extra))
        self.assertEqual(rc, 0)
        with open(self.total) as handle:
            total = next(csv.DictReader(handle))
        with open(self.detailed) as handle:
            detailed = {row["RSID"]: row for row in csv.DictReader(handle)}
        return {key: None if value == "NA" else float(value) for key, value in total.items()}, detailed

    def _assert_total(self, total, expected):
        for key, value in expected.items():
            self.assertAlmostEqual(total[key], value, places=12, msg=key)

    def test_total_scores(self):
        total, _ = self._run()
        self._assert_total(
            total,
            {
                "SumPRS": 7.1,
                "CountNA": 2,
                "MaxNeg": 7.1,
                "MaxPos": 9.5,
                "MaxRange": 2.4,
                "ExpectedScore": 7.4,
                "DiffSumVsExpected": -0.3,
            },
        )

    def test_detailed_scores(self):
        _, detailed = self._run()
        effect_alleles = {rsid: row["Number_of_effect_alleles"] for rsid, row in detailed.items()}
        self.assertEqual(
            effect_alleles,
            {
                "rs1": "1",
                "rs2": "2",
                "rs3": "0",
                "rs4": "NA",
                "rs5": "NA",
                "rs6": "1",
                "rs7": "2",
                "rs9": "2",
            },
        )
        # The block-derived reference call is reported as such
        self.assertEqual(detailed["rs3"]["VCF_Ref"], "R")
        # The multi-allelic pair and the proxy pair are each counted once
        self.assertNotIn("rs8", detailed)
        self.assertNotIn("rs10", detailed)

    def test_gzipped_vcf(self):
        gzipped = os.path.join(self.tempdir, "HD832.HES45_T.g.vcf.gz")
        with open(self.vcf, "rb") as plain, gzip.open(gzipped, "wb") as packed:
            packed.write(plain.read())
        total, _ = self._run(vcf=gzipped)
        self._assert_total(total, {"SumPRS": 7.1, "CountNA": 2})

    def test_r_compat(self):
        # prs.Daniel.R misreads '.|1' as a homozygote and lets the block override the no-call
        total, detailed = self._run(extra=["--r-compat"])
        self._assert_total(total, {"SumPRS": 8.9, "CountNA": 1})
        self.assertEqual(detailed["rs5"]["Number_of_effect_alleles"], "2")

    def test_no_proxy_dedup(self):
        # Both proxy rows are kept: the missing proxy is NA and the replaced site scores 2 x 5.0
        total, detailed = self._run(extra=["--no-proxy-dedup"])
        self._assert_total(total, {"SumPRS": 13.1, "CountNA": 3})
        self.assertEqual(detailed["rs9"]["Number_of_effect_alleles"], "NA")
        self.assertEqual(detailed["rs10"]["Number_of_effect_alleles"], "2")

    def test_incomplete_genotypes(self):
        from prs_daniel import is_incomplete_genotype

        for genotype in [None, "", ".", "./.", ".|.", "./1", ".|1", "1/.", "1|.", "0/.", "./2", "2|.", "./././."]:
            self.assertTrue(is_incomplete_genotype(genotype), genotype)
        for genotype in ["0/0", "0/1", "1/1", "0|1", "1|0", "1|1", "1/2", "2/2"]:
            self.assertFalse(is_incomplete_genotype(genotype), genotype)
        # The R list of missing genotypes only catches unphased spellings
        self.assertTrue(is_incomplete_genotype("./1", r_compat=True))
        self.assertFalse(is_incomplete_genotype(".|1", r_compat=True))

    def test_outputs_written_to_given_paths(self):
        self.detailed = os.path.join(self.tempdir, "a", "b", "detailed.csv")
        self.total = os.path.join(self.tempdir, "c", "total.csv")
        self._run()
        self.assertTrue(os.path.isfile(self.detailed))
        self.assertTrue(os.path.isfile(self.total))
        self.assertEqual(sorted(os.listdir(self.tempdir)), ["a", "c"])

    def test_bad_output_paths_refused(self):
        from prs_daniel import main

        with self.assertRaises(SystemExit):
            main(["-s", self.score_file, "-v", self.vcf, "-d", self.detailed, "-t", self.detailed])
        with self.assertRaises(SystemExit):
            main(["-s", self.score_file, "-v", self.vcf, "-d", self.tempdir, "-t", self.total])
