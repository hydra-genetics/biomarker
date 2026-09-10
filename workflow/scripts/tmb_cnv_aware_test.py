import tempfile
import os
import unittest


class TestUnitUtils(unittest.TestCase):
    def setUp(self):
        self.vcf = ".tests/unit/HD832.HES45_T.background_annotation.vcf.gz"
        self.cns = ".tests/unit/HD832.HES45_T.loh.cns"
        self.artifacts = ""
        self.background_panel = ""
        self.filter_genes = ".tests/unit/tmb_filter_genes.txt"
        self.filter_nr_observations = 1
        self.filter_regions = []
        self.dp_limit = 200
        self.vd_limit = 20
        self.af_lower_limit = 0.05
        self.af_upper_limit = 0.95
        self.af_germline_lower_limit = 0.47
        self.af_germline_upper_limit = 0.53
        self.gnomad_limit = 0.0001
        self.db1000g_limit = 0.0001
        self.background_sd_limit = 5
        self.nr_avg_germline_snvs = 2.0
        self.nssnv_tmb_correction = 0.84
        self.variant_type_list = ["missense_variant", "stop_gained", "stop_lost"]
        self.cn_baf_tolerance = 0.05
        self.pmean_limit = 15
        self.nm_limit = 4

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
        """
        Same fixtures as tmb_test.py::test_tmb - none of the three counted
        variants (all on chr1, positions 27088708/27088710/27088711) are
        covered by a segment in the real .loh.cns fixture (both its chr1
        segments end at or before 2460407), and both segments have an empty
        baf anyway, so the CN-aware germline check must fall back to the
        same fixed diploid window as the original tmb() - and none of the
        three variants have low PMEAN (31.1) or high NM (1.5), so the new
        QC gates don't remove them either. Result must match tmb_test.py
        exactly: same inputs, same effective filtering, same output.
        """
        from tmb_cnv_aware import tmb

        out_tmb = open(os.path.join(self.tempdir, "HD832.HES45_T.TMB.txt"), "w")

        tmb(
            self.vcf, self.artifacts, self.background_panel, out_tmb, self.filter_genes, self.filter_nr_observations,
            self.filter_regions, self.dp_limit, self.vd_limit, self.af_lower_limit, self.af_upper_limit,
            self.af_germline_lower_limit, self.af_germline_upper_limit, self.gnomad_limit, self.db1000g_limit,
            self.background_sd_limit, self.nr_avg_germline_snvs, self.nssnv_tmb_correction, self.variant_type_list,
            self.cns, self.cn_baf_tolerance, self.pmean_limit, self.nm_limit,
        )
        out_tmb.close()

        result_file = open(os.path.join(self.tempdir, "HD832.HES45_T.TMB.txt"))

        header = True
        i = 0
        result = []
        for line in result_file:
            result.append(line)
            i += 1
            if i == 2:
                break

        test_table = {
            "TMB:": "0.84",
            "Number of variants:": "3",
        }

        self._test_tmb(test_table, result)

    def test_tmb_pmean_nm_gate(self):
        """
        Sanity-check that the new PMEAN/NM gates actually do something: with
        a pmean_limit above the real fixture's PMEAN (31.1) all three
        previously-counted variants must now be excluded.
        """
        from tmb_cnv_aware import tmb

        out_tmb = open(os.path.join(self.tempdir, "HD832.HES45_T.pmean_gate.TMB.txt"), "w")

        tmb(
            self.vcf, self.artifacts, self.background_panel, out_tmb, self.filter_genes, self.filter_nr_observations,
            self.filter_regions, self.dp_limit, self.vd_limit, self.af_lower_limit, self.af_upper_limit,
            self.af_germline_lower_limit, self.af_germline_upper_limit, self.gnomad_limit, self.db1000g_limit,
            self.background_sd_limit, self.nr_avg_germline_snvs, self.nssnv_tmb_correction, self.variant_type_list,
            self.cns, self.cn_baf_tolerance, 40, self.nm_limit,
        )
        out_tmb.close()

        result_file = open(os.path.join(self.tempdir, "HD832.HES45_T.pmean_gate.TMB.txt"))
        i = 0
        result = []
        for line in result_file:
            result.append(line)
            i += 1
            if i == 2:
                break

        test_table = {
            "TMB:": "0",
            "Number of variants:": "0",
        }
        self._test_tmb(test_table, result)

    def test_read_cns_segments(self):
        from tmb_cnv_aware import read_cns_segments

        cns_dict = read_cns_segments(self.cns)

        try:
            self.assertEqual([150500, 934993, None], cns_dict["chr1"][0])
            self.assertEqual([935853, 2460407, None], cns_dict["chr1"][1])
        except AssertionError as e:
            print(f"Failed reading cns segments. {cns_dict}")
            raise e

    def test_lookup_local_baf(self):
        from tmb_cnv_aware import read_cns_segments, lookup_local_baf

        cns_dict = read_cns_segments(self.cns)

        # Covered by a segment, but that segment's baf is empty
        baf = lookup_local_baf(cns_dict, "chr1", 200000)
        try:
            self.assertIsNone(baf)
        except AssertionError as e:
            print(f"Failed looking up empty-baf segment. {baf}")
            raise e

        # Not covered by any segment
        baf = lookup_local_baf(cns_dict, "chr1", 27088708)
        try:
            self.assertIsNone(baf)
        except AssertionError as e:
            print(f"Failed looking up uncovered position. {baf}")
            raise e

        # A synthetic segment with a real, fitted baf
        synthetic_cns = os.path.join(self.tempdir, "synthetic.loh.cns")
        with open(synthetic_cns, "w") as f:
            f.write("chromosome\tstart\tend\tgene\tlog2\tbaf\tci_hi\tci_lo\tcn\tcn1\tcn2\tdepth\tprobes\tweight\n")
            f.write("chr2\t1000\t2000\t-\t0.58\t0.75\t0.6\t0.5\t3\t2\t1\t500\t10\t8\n")
        cns_dict2 = read_cns_segments(synthetic_cns)
        baf = lookup_local_baf(cns_dict2, "chr2", 1500)
        try:
            self.assertEqual(0.75, baf)
        except AssertionError as e:
            print(f"Failed looking up populated-baf segment. {baf}")
            raise e

    def test_is_likely_germline(self):
        from tmb_cnv_aware import is_likely_germline

        # No local baf - fall back to the fixed diploid window
        try:
            self.assertTrue(is_likely_germline(0.50, None, 0.47, 0.53, 0.05))
            self.assertFalse(is_likely_germline(0.75, None, 0.47, 0.53, 0.05))
        except AssertionError as e:
            print("Failed fixed-window fallback case")
            raise e

        # Local baf=0.75 (a gained/LOH segment) - germline SNPs cluster near
        # 0.75 or its mirror 0.25, not near 0.5
        try:
            self.assertTrue(is_likely_germline(0.75, 0.75, 0.47, 0.53, 0.05))
            self.assertTrue(is_likely_germline(0.24, 0.75, 0.47, 0.53, 0.05))
            self.assertFalse(is_likely_germline(0.50, 0.75, 0.47, 0.53, 0.05))
        except AssertionError as e:
            print("Failed baf/1-baf clustering case")
            raise e
