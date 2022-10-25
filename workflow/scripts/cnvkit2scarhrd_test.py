import tempfile
import os
import unittest
import gzip


class TestUnitUtils(unittest.TestCase):
    def setUp(self):
        self.in_seg = ".tests/unit/HD832.HES45_T.loh.cns"

        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        pass

    def _test_cnvkit2scarhrd(self, test_table, variants):
        for variant in variants:
            columns = variant.strip().split("\t")
            try:
                self.assertEqual(
                    test_table["{}\t{}\t{}".format(columns[0], columns[1], columns[2])],
                    "{}\t{}\t{}\t{}".format(columns[4], columns[5], columns[6], columns[7]),
                )
            except AssertionError as e:
                print("Failed cnvkit 2 scarhrd conversion of: " + str(variant))
                raise e

    def test_cnvkit2scarhrd(self):
        from cnvkit2scarhrd import cnvkit_2_scarhrd

        out_seg = open(os.path.join(self.tempdir, "HD832.HES45_T.scarhrd.cns"), "w")

        # Convert cnvkit segmentation to be compatible with scarHRD
        cnvkit_2_scarhrd(self.in_seg, out_seg)
        out_seg.close()

        result_file = open(os.path.join(self.tempdir, "HD832.HES45_T.scarhrd.cns"))

        header = True
        result = []
        for line in result_file:
            print(line)
            if header:
                header = False
                continue
            result.append(line)

        test_table = {
            "HD832.HES45_T\tchr1\t150500": '0\t0\t0\tNA',  # deletion
            "HD832.HES45_T\tchr1\t935853": '3\t2\t1\tNA',  # duplication without cn1
        }

        self._test_cnvkit2scarhrd(test_table, result)
