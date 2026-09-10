
import gzip


def read_bedfile(filter_regions_dict, in_bed):
    for line in in_bed:
        columns = line.strip().split("\t")
        chrom = columns[0]
        if chrom not in filter_regions_dict:
            filter_regions_dict[chrom] = []
        filter_regions_dict[chrom].append([int(columns[1]), int(columns[2])])
    return filter_regions_dict


def read_cns_segments(cns_path):
    """
    Read a CNVkit purity/BAF-aware `.loh.cns` file (produced by `cnvkit call
    --purity <tc> --vcf <germline.vcf>`) into a per-chromosome sorted list of
    [start, end, baf]. `baf` is CNVkit's per-segment fitted B-allele-frequency
    (folded to [0.5, 1.0], on the same raw AF scale as any individual
    variant's AF - not further purity/ploidy-rescaled the way cn/cn1/cn2 are),
    or None when the segment had too few germline SNPs for CNVkit to fit one
    (routinely empty in a targeted panel, not a rare edge case).

    param cns_path: path to a CNVkit .loh.cns file
    return: dict {chrom: [[start, end, baf_or_None], ...]}, sorted by start
    """
    cns_dict = {}
    with open(cns_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        chrom_i = header.index("chromosome")
        start_i = header.index("start")
        end_i = header.index("end")
        baf_i = header.index("baf")
        for line in f:
            if not line.strip():
                continue
            columns = line.rstrip("\n").split("\t")
            chrom = columns[chrom_i]
            start = int(columns[start_i])
            end = int(columns[end_i])
            baf = float(columns[baf_i]) if columns[baf_i] != "" else None
            cns_dict.setdefault(chrom, []).append([start, end, baf])
    for chrom in cns_dict:
        cns_dict[chrom].sort(key=lambda seg: seg[0])
    return cns_dict


def lookup_local_baf(cns_dict, chrom, pos):
    """
    Find the fitted baf of the CNVkit segment covering a genomic position.

    param cns_dict: dict created by read_cns_segments
    param chrom: chromosome name (must match the .loh.cns file's naming)
    param pos: 1-based genomic position
    return: baf (float), or None if uncovered or the segment's baf was empty
    """
    for start, end, baf in cns_dict.get(chrom, []):
        if start <= pos <= end:
            return baf
    return None


def is_likely_germline(af, baf, af_germline_lower_limit, af_germline_upper_limit, baf_tolerance):
    """
    Decide whether an observed AF looks like a germline call. Without local
    CN/BAF information this is the original fixed-window heuristic (assumes a
    plain diploid heterozygous SNP at ~50% AF). When CNVkit fitted a reliable
    baf for the covering segment, a germline SNP there clusters near baf or
    its mirror (1 - baf) instead of near 0.5, so both are checked.

    param af: observed variant allele fraction
    param baf: local segment's fitted CNVkit baf, or None if unavailable
    param af_germline_lower_limit: lower bound of the diploid germline window
    param af_germline_upper_limit: upper bound of the diploid germline window
    param baf_tolerance: +/- tolerance around baf / (1 - baf)
    return: True if af looks like a germline call
    """
    if baf is None:
        return af_germline_lower_limit <= af <= af_germline_upper_limit
    return abs(af - baf) <= baf_tolerance or abs(af - (1 - baf)) <= baf_tolerance


def tmb(
    vcf, artifacts_filename, background_panel_filename, output_tmb, filter_genes_filename, filter_nr_observations,
    filter_regions, dp_limit, ad_limit, af_lower_limit, af_upper_limit, af_germline_lower_limit,
    af_germline_upper_limit, gnomad_limit, db1000g_limit, background_sd_limit, nr_avg_germline_snvs,
    nssnv_tmb_correction, variant_type_list, cns_filename, cn_baf_tolerance, pmean_limit, nm_limit
):

    '''Problematic regions'''
    filter_regions_dict = {}
    for filter_bed_file in filter_regions:
        filter_regions_dict = read_bedfile(filter_regions_dict, open(filter_bed_file))

    FFPE_SNV_artifacts = {}
    if artifacts_filename == "":
        filter_nr_observations = 1
    else:
        artifacts = open(artifacts_filename)
        header = True
        for line in artifacts:
            columns = line.strip().split("\t")
            if header:
                header = False
                continue
            chrom = columns[0]
            pos = columns[1]
            key = chrom + "_" + pos
            type = columns[2]
            if type != "SNV":
                continue
            max_observations = 0
            i = 0
            for observation in columns[3:]:
                if i % 3 == 0:
                    if int(observation) > max_observations:
                        max_observations = int(observation)
                i += 1
            FFPE_SNV_artifacts[key] = max_observations

    '''Background'''
    gvcf_panel_dict = {}
    if background_panel_filename != "":
        background_panel = open(background_panel_filename)
        next(background_panel)
        for line in background_panel:
            columns = line.strip().split()
            chrom = columns[0]
            pos = columns[1]
            key = chrom + "_" + pos
            median = float(columns[2])
            sd = float(columns[3])
            gvcf_panel_dict[key] = [median, sd]

    '''Filter genes'''
    filter_gene_dict = {}
    if filter_genes_filename != "":
        filter_genes = open(filter_genes_filename)
        for line in filter_genes:
            filter_gene_dict[line.strip()] = ""

    '''Local copy-number/BAF, for CN-aware germline detection'''
    cns_dict = read_cns_segments(cns_filename)

    nr_SNV_TMB = 0
    header = True
    prev_pos = ""
    prev_chrom = ""
    TMB_SNV = []
    vep_dict = {}
    with gzip.open(vcf, 'rt') as vcf_infile:
        file_content = vcf_infile.read().split("\n")
        for line in file_content:
            if line == "":
                continue
            if header:
                if line[:6] == "#CHROM":
                    header = False
                if line[:14] == "##INFO=<ID=CSQ":
                    vep_string = line.split("Format: ")[1]
                    vep_list = vep_string.split("|")
                    i = 0
                    for v in vep_list:
                        vep_dict[v] = i
                        i += 1
                continue
            lline = line.strip().split("\t")
            chrom = lline[0]
            pos = lline[1]
            if chrom == prev_chrom and pos == prev_pos:
                continue
            prev_pos = pos
            prev_chrom = chrom
            key = chrom + "_" + pos
            ref = lline[3]
            alt = lline[4]
            filter = lline[6]
            INFO = lline[7]
            INFO_list = INFO.split(";")
            AF_index = -1
            Caller_index = 0
            PMEAN_index = -1
            NM_index = -1
            i = 0
            for info in INFO_list:
                if info[:3] == "AF=":
                    AF_index = i
                if info[:8] == "CALLERS=":
                    Caller_index = i
                if info[:6] == "PMEAN=":
                    PMEAN_index = i
                if info[:3] == "NM=":
                    NM_index = i
                i += 1
            if AF_index == -1:
                continue
            AF = float(INFO_list[AF_index][3:])
            PMEAN = float(INFO_list[PMEAN_index][6:]) if PMEAN_index != -1 else None
            NM = float(INFO_list[NM_index][3:]) if NM_index != -1 else None
            Callers = INFO_list[Caller_index]
            nr_Callers = len(Callers.split(","))
            if nr_Callers < 2:
                continue
            VEP_INFO = INFO.split("CSQ=")[1].split("|")
            Variant_type = VEP_INFO[vep_dict["Consequence"]].split("&")
            db1000G = VEP_INFO[vep_dict["AF"]]
            if db1000G == "":
                db1000G = 0
            else:
                db1000G = float(db1000G)
            GnomAD = VEP_INFO[vep_dict["MAX_AF"]]
            if GnomAD == "":
                GnomAD = 0
            else:
                GnomAD = float(GnomAD)
            FORMAT = lline[8].split(":")
            AD_index = 0
            DP_index = 0
            VD_index = 0
            i = 0
            for f in FORMAT:
                if f == "AD":
                    AD_index = i
                elif f == "DP":
                    DP_index = i
                elif f == "VD":
                    VD_index = i
                i += 1
            DATA = lline[9].split(":")
            AD = DATA[AD_index].split(",")
            if len(AD) == 2:
                VD = int(AD[1])
            else:
                VD = int(DATA[VD_index])
            DP = int(DATA[DP_index])

            # Only SNVs
            Observations = 0
            if not(len(ref) == 1 and len(alt) == 1):
                continue

            # Artifact observations
            if key in FFPE_SNV_artifacts:
                Observations = FFPE_SNV_artifacts[key]

            # Gene name
            gene_name = VEP_INFO[vep_dict["SYMBOL"]]

            # Local copy-number/BAF-aware germline check
            local_baf = lookup_local_baf(cns_dict, chrom, int(pos))
            looks_germline = is_likely_germline(
                AF, local_baf, af_germline_lower_limit, af_germline_upper_limit, cn_baf_tolerance
            )

            # Read-quality signature: mean read position and mean mismatches per read
            ok_pmean = PMEAN is None or PMEAN >= pmean_limit
            ok_nm = NM is None or NM <= nm_limit

            # TMB
            if (DP > dp_limit and VD > ad_limit and AF >= af_lower_limit and AF <= af_upper_limit and
                    not looks_germline and ok_pmean and ok_nm and
                    GnomAD <= gnomad_limit and db1000G <= db1000g_limit and
                    Observations < filter_nr_observations and INFO.find("Complex") == -1) and gene_name not in filter_gene_dict:
                panel_median = 1000
                panel_sd = 1000
                pos_sd = 1000
                key2 = key[3:]
                if key2 in gvcf_panel_dict:
                    panel_median = gvcf_panel_dict[key2][0]
                    panel_sd = gvcf_panel_dict[key2][1]
                if panel_sd > 0.0:
                    pos_sd = (AF - panel_median) / panel_sd
                if background_panel_filename == "" or pos_sd > background_sd_limit:
                    ok_variant_type = False
                    for vt in Variant_type:
                        if vt in variant_type_list:
                            ok_variant_type = True
                    print(Variant_type, ok_variant_type, variant_type_list)
                    if ok_variant_type:
                        filter_region = False
                        if chrom in filter_regions_dict:
                            for region in filter_regions_dict[chrom]:
                                if int(pos) >= region[0] and int(pos) <= region[1]:
                                    filter_region = True
                                    break
                        if not filter_region:
                            nr_SNV_TMB += 1
                            TMB_SNV.append([line, panel_median, panel_sd, AF, pos_sd])

    TMB = (nr_SNV_TMB - nr_avg_germline_snvs) * nssnv_tmb_correction
    if TMB < 0:
        TMB = 0
    output_tmb.write("TMB:\t" + str(TMB) + "\n")
    output_tmb.write("Number of variants:\t" + str(nr_SNV_TMB) + "\n")
    for variant in TMB_SNV:
        if background_panel_filename == "":
            output_tmb.write(variant[0] + "\n")
        else:
            output_tmb.write(
                variant[0] + "\t" + "{:.4f}".format(variant[1]) + "\t" + "{:.4f}".format(variant[2]) + "\t" +
                "\t" + "{:.4f}".format(variant[3]) + "\t" + "{:.2f}".format(variant[4]) + "\n"
            )


if __name__ == "__main__":
    log = snakemake.log_fmt_shell(stdout=False, stderr=True)

    tmb(
        snakemake.input.vcf,
        snakemake.params.artifacts,
        snakemake.params.background_panel,
        open(snakemake.output.tmb, "w"),
        snakemake.params.filter_genes,
        snakemake.params.filter_nr_observations,
        snakemake.params.filter_regions,
        snakemake.params.dp_limit,
        snakemake.params.vd_limit,
        snakemake.params.af_lower_limit,
        snakemake.params.af_upper_limit,
        snakemake.params.af_germline_lower_limit,
        snakemake.params.af_germline_upper_limit,
        snakemake.params.gnomad_limit,
        snakemake.params.db1000g_limit,
        snakemake.params.background_sd_limit,
        snakemake.params.nr_avg_germline_snvs,
        snakemake.params.nssnv_tmb_correction,
        snakemake.params.variant_type_list,
        snakemake.input.cns,
        snakemake.params.cn_baf_tolerance,
        snakemake.params.pmean_limit,
        snakemake.params.nm_limit,
    )
