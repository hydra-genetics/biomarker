import tempfile
import os
import unittest
import gzip


class TestUnitUtils(unittest.TestCase):
    def setUp(self):
        self.in_hrd = open(".tests/unit/HD832.HES45_T_HRDresults.txt")

        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        pass

    def _test_fix_scarhrd(self, test_data, variants):
        for variant in variants:
            columns = variant.strip().split("\t")
            try:
                self.assertEqual(variant, test_data)
            except AssertionError as e:
                print("Failed scarHRD of: " + str(variant))
                raise e

    def test_fix_scarhrd(self):
        from fix_scarhrd_output import fix_scarhrd

        out_hrd = open(os.path.join(self.tempdir, "HD832.HES45_T.scarhrd_cnvkit_score.txt"), "w")

        # Run scarHDR
        fix_scarhrd(self.in_hrd, out_hrd)

        result_file = open(os.path.join(self.tempdir, "HD832.HES45_T.scarhrd_cnvkit_score.txt"))

        header = True
        result = []
        for line in result_file:
            if header:
                header = False
                continue
            result.append(line)

        test_data = "6\t1\t2\t3\n"  # HRD score

        self._test_fix_scarhrd(test_data, result)
