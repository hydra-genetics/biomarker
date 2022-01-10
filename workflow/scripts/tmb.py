
import gzip

vcf = snakemake.input.vcf
artifacts = open(snakemake.input.artifacts)
background_panel_filename = snakemake.input.background_panel
background_run = open(snakemake.input.background_run)
output_tmb = open(snakemake.output.tmb, "w")
filter_nr_observations = snakemake.params.filter_nr_observations
DP_limit = snakemake.params.DP_limit
VD_limit = snakemake.params.VD_limit
AF_lower_limit = snakemake.params.AF_lower_limit
AF_upper_limit = snakemake.params.AF_upper_limit
GnomAD_limit = snakemake.params.GnomAD_limit
db1000G_limit = snakemake.params.db1000G_limit
Background_sd_limit = snakemake.params.Background_sd_limit
nsSNV_TMB_correction = snakemake.params.nsSNV_TMB_correction
nsSNV_sSNV_TMB_correction = snakemake.params.nsSNV_sSNV_TMB_correction


FFPE_SNV_artifacts = {}
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
    for nr_variant in columns[3:]:
        if int(nr_variant) > max_observations:
            max_observations = int(nr_variant)
    FFPE_SNV_artifacts[key] = max_observations


'''Background'''
gvcf_panel_dict = {}
gvcf_run_dict = {}
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
next(background_run)
for line in background_run:
    columns = line.strip().split()
    chrom = columns[0]
    pos = columns[1]
    key = chrom + "_" + pos
    median = float(columns[2])
    gvcf_run_dict[key] = median


nr_nsSNV_TMB = 0
nr_sSNV_TMB = 0
header = True
prev_pos = ""
prev_chrom = ""
TMB_nsSNV = []
TMB_sSNV = []
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
        if INFO[:3] == "AA=":
            continue
        INFO_list = INFO.split(";")
        AF_index = 0
        Caller_index = 0
        i = 0
        for info in INFO_list:
            if info[:3] == "AF=":
                AF_index = i
            if info[:8] == "CALLERS=":
                Caller_index = i
            i += 1
        AF = float(INFO_list[AF_index][3:])
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

        # Artifact observations
        Observations = 0
        if len(ref) == 1 and len(alt) == 1:
            if key in FFPE_SNV_artifacts:
                Observations = FFPE_SNV_artifacts[key]
        else:
            continue

        # TMB
        if (DP > DP_limit and VD > VD_limit and AF >= AF_lower_limit and AF <= AF_upper_limit and
                GnomAD <= GnomAD_limit and db1000G <= db1000G_limit and
                Observations < filter_nr_observations and INFO.find("Complex") == -1):
            if len(ref) == 1 and len(alt) == 1:
                panel_median = 1000
                panel_sd = 1000
                run_median = 1000
                pos_sd = 1000
                key2 = key[3:]
                if key2 in gvcf_panel_dict:
                    panel_median = gvcf_panel_dict[key2][0]
                    panel_sd = gvcf_panel_dict[key2][1]
                if key2 in gvcf_run_dict:
                    run_median = gvcf_run_dict[key2]
                if panel_sd > 0.0:
                    pos_sd = (AF - panel_median) / panel_sd
                if pos_sd > Background_sd_limit:
                    if ("missense_variant" in Variant_type or
                            "stop_gained" in Variant_type or
                            "stop_lost" in Variant_type):
                        nr_nsSNV_TMB += 1
                        TMB_nsSNV.append([line, panel_median, panel_sd, run_median, AF, pos_sd])
                    elif "synonymous_variant" in Variant_type:
                        nr_sSNV_TMB += 1
                        TMB_sSNV.append([line, panel_median, panel_sd, run_median, AF, pos_sd])

nsTMB = nr_nsSNV_TMB * nsSNV_TMB_correction
total_TMB = (nr_sSNV_TMB + nr_nsSNV_TMB) * nsSNV_sSNV_TMB_correction
output_tmb.write("nsSNV TMB:\t" + str(nsTMB) + "\n")
output_tmb.write("nsSNV variants:\t" + str(nr_nsSNV_TMB) + "\n")
output_tmb.write("TMB:\t" + str(total_TMB) + "\n")
output_tmb.write("SNV in coding regions:\t" + str(nr_sSNV_TMB + nr_nsSNV_TMB) + "\nList of variants:\n")
for TMB in TMB_nsSNV:
    output_tmb.write(
        TMB[0].strip() + "\t" + "{:.4f}".format(TMB[1]) + "\t" + "{:.4f}".format(TMB[2]) + "\t" + "{:.4f}".format(TMB[3]) +
        "\t" + "{:.4f}".format(TMB[4]) + "\t" + "{:.2f}".format(TMB[5]) + "\n"
    )
for TMB in TMB_sSNV:
    output_tmb.write(
        TMB[0].strip() + "\t" + "{:.4f}".format(TMB[1]) + "\t" + "{:.4f}".format(TMB[2]) + "\t" + "{:.4f}".format(TMB[3]) +
        "\t" + "{:.4f}".format(TMB[4]) + "\t" + "{:.2f}".format(TMB[5]) + "\n"
    )
