import tempfile
import os
import unittest
import gzip


class TestUnitUtils(unittest.TestCase):
    def setUp(self):
        self.vcf = ".tests/unit/HD832.HES45_T.background_annotation.vcf.gz"
        self.artifacts = open(".tests/unit/artifact_panel_chr1.tsv")
        self.background_panel = ".tests/unit/background_panel_HES45.tsv"
        self.filter_nr_observations = 1
        self.dp_limit = 200
        self.vd_limit = 10
        self.af_lower_limit = 0.05
        self.af_upper_limit = 0.45
        self.gnomad_limit = 0.0001
        self.db1000g_limit = 0.0001
        self.background_sd_limit = 5
        self.nssnv_tmb_correction = 0.78
        self.nssnv_ssnv_tmb_correction = 0.57

        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        pass

    def _test_tmb(self, test_table, variants):
        for variant in variants:
            columns = variant.strip().split("\t")
            try:
                self.assertEqual(test_table[columns[0]], columns[1])
            except AssertionError as e:
                print("Failed TMB calculation of: " + str(variant))
                raise e

    def test_tmb(self):
        from tmb import tmb

        out_tmb = open(os.path.join(self.tempdir, "HD832.HES45_T.TMB.txt"), "w")

        # Run scarHDR
        tmb(
            self.vcf, self.artifacts, self.background_panel, out_tmb,
            self.filter_nr_observations, self.dp_limit, self.vd_limit, self.af_lower_limit,
            self.af_upper_limit, self.gnomad_limit, self.db1000g_limit, self.background_sd_limit,
            self.nssnv_tmb_correction, self.nssnv_ssnv_tmb_correction,
        )
        out_tmb.close()

        result_file = open(os.path.join(self.tempdir, "HD832.HES45_T.TMB.txt"))

        header = True
        i = 0
        result = []
        for line in result_file:
            result.append(line)
            i += 1
            if i == 4:
                break

        test_table = {
            "nsSNV TMB:": "0.78",
            "nsSNV variants:": "1",
            "TMB:": "0.57",
            "SNV in coding regions:": "1",
        }

        self._test_tmb(test_table, result)
