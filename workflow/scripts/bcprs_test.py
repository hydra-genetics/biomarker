import csv
import gzip
import os
import runpy
import shutil
import sys
import tempfile
import unittest
from types import SimpleNamespace


class TestUnitUtils(unittest.TestCase):
    """Hand-computed checks of bcprs.py on the fixtures in .tests/unit.

    Every marker in the score file exercises one feature of the script:
        1:100        heterozygous SNP                       1 x  0.5
        1:200        phased homozygous SNP                  2 x -0.2
        1:250        inside a hom-ref gVCF block            0 x  0.7 (called, not missing)
        1:300        absent from the gVCF                   NA
        1:400        no-call '.|1' inside a hom-ref block   NA (2 x 0.9 with r_compat)
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

    def _run(self, vcf=None, **settings):
        from bcprs import bcprs

        bcprs(self.score_file, vcf or self.vcf, self.detailed, self.total, **settings)
        return self._read_outputs()

    def _read_outputs(self):
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
        total, detailed = self._run(r_compat=True)
        self._assert_total(total, {"SumPRS": 8.9, "CountNA": 1})
        self.assertEqual(detailed["rs5"]["Number_of_effect_alleles"], "2")

    def test_no_proxy_dedup(self):
        # Both proxy rows are kept: the missing proxy is NA and the replaced site scores 2 x 5.0
        total, detailed = self._run(no_proxy_dedup=True)
        self._assert_total(total, {"SumPRS": 13.1, "CountNA": 3})
        self.assertEqual(detailed["rs9"]["Number_of_effect_alleles"], "NA")
        self.assertEqual(detailed["rs10"]["Number_of_effect_alleles"], "2")

    def test_incomplete_genotypes(self):
        from bcprs import is_incomplete_genotype

        for genotype in [None, "", ".", "./.", ".|.", "./1", ".|1", "1/.", "1|.", "0/.", "./2", "2|.", "./././."]:
            self.assertTrue(is_incomplete_genotype(genotype), genotype)
        for genotype in ["0/0", "0/1", "1/1", "0|1", "1|0", "1|1", "1/2", "2/2"]:
            self.assertFalse(is_incomplete_genotype(genotype), genotype)
        # The R list of missing genotypes only catches unphased spellings
        self.assertTrue(is_incomplete_genotype("./1", r_compat=True))
        self.assertFalse(is_incomplete_genotype(".|1", r_compat=True))

    def test_snakemake_entry_point(self):
        # Run the script the way the rule's script: directive does, with an injected snakemake object
        log = os.path.join(self.tempdir, "HD832.HES45_T.log")
        snakemake = SimpleNamespace(
            input=SimpleNamespace(score_file=self.score_file, vcf=self.vcf),
            output=SimpleNamespace(detailed=self.detailed, total=self.total),
            params=SimpleNamespace(no_proxy_dedup=False, r_compat=True),
            log=[log],
        )
        stdout, stderr = sys.stdout, sys.stderr
        try:
            runpy.run_path(
                os.path.join(os.path.dirname(os.path.abspath(__file__)), "bcprs.py"),
                init_globals={"snakemake": snakemake},
                run_name="__main__",
            )
        finally:
            if sys.stdout is not stdout:
                sys.stdout.close()
            sys.stdout, sys.stderr = stdout, stderr
        total, _ = self._read_outputs()
        self._assert_total(total, {"SumPRS": 8.9, "CountNA": 1})
        with open(log) as handle:
            self.assertIn("SumPRS = 8.9 (CountNA = 1)", handle.read())
