import tempfile
import os
import unittest
import gzip


class TestUnitUtils(unittest.TestCase):
    def setUp(self):
        self.vcf = ".tests/unit/HD832.HES45_T.background_annotation.vcf.gz"
        self.dp_limit = 200
        self.vd_limit = 10
        self.af_lower_limit = 0.05
        self.af_upper_limit = 0.95
        self.af_germline_lower_limit = 0.47
        self.af_germline_upper_limit = 0.53
        self.gnomad_limit = 0.0001
        self.db1000g_limit = 0.0001
        self.nr_avg_germline_snvs = 2.0
        self.nssnv_tmb_correction = 0.84

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

        tmb(
            self.vcf, out_tmb, self.dp_limit, self.vd_limit,
            self.af_lower_limit, self.af_upper_limit, self.af_germline_lower_limit, self.af_germline_upper_limit,
            self.gnomad_limit, self.db1000g_limit, self.nr_avg_germline_snvs, self.nssnv_tmb_correction,
        )
        out_tmb.close()

        result_file = open(os.path.join(self.tempdir, "HD832.HES45_T.TMB.txt"))

        header = True
        i = 0
        result = []
        for line in result_file:
            result.append(line)
            print(line)
            i += 1
            if i == 2:
                break

        test_table = {
            "TMB:": "0.84",
            "Number of variants:": "3",
        }

        self._test_tmb(test_table, result)
