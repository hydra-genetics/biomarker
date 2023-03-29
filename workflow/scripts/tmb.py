
import gzip


def tmb(
    vcf, output_tmb, dp_limit, ad_limit, af_lower_limit, af_upper_limit,
    af_germline_lower_limit, af_germline_upper_limit, gnomad_limit, db1000g_limit,
    nr_avg_germline_snvs, nssnv_tmb_correction
):

    nr_nsSNV_TMB = 0
    nr_sSNV_TMB = 0
    header = True
    prev_pos = ""
    prev_chrom = ""
    TMB_nsSNV = []
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
            i = 0
            for info in INFO_list:
                if info[:3] == "AF=":
                    AF_index = i
                if info[:8] == "CALLERS=":
                    Caller_index = i
                i += 1
            if AF_index == -1:
                continue
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

            # Only SNVs
            if not (len(ref) == 1 and len(alt) == 1):
                continue

            # TMB
            if (DP > dp_limit and VD > ad_limit and AF >= af_lower_limit and AF <= af_upper_limit and
                    (AF <= af_germline_lower_limit or AF >= af_germline_upper_limit) and
                    GnomAD <= gnomad_limit and db1000G <= db1000g_limit and
                    INFO.find("Complex") == -1):
                if ("missense_variant" in Variant_type or
                        "stop_gained" in Variant_type or
                        "stop_lost" in Variant_type):
                    nr_nsSNV_TMB += 1
                    TMB_nsSNV.append(line)

    nsTMB = (nr_nsSNV_TMB - nr_avg_germline_snvs) * nssnv_tmb_correction
    if nsTMB < 0:
        nsTMB = 0
    output_tmb.write("TMB:\t" + str(nsTMB) + "\n")
    output_tmb.write("Number of variants:\t" + str(nr_nsSNV_TMB) + "\n")
    output_tmb.write("List of variants:\n")
    for TMB in TMB_nsSNV:
        output_tmb.write(TMB)


if __name__ == "__main__":
    log = snakemake.log_fmt_shell(stdout=False, stderr=True)

    tmb(
        snakemake.input.vcf,
        open(snakemake.output.tmb, "w"),
        snakemake.params.dp_limit,
        snakemake.params.vd_limit,
        snakemake.params.af_lower_limit,
        snakemake.params.af_upper_limit,
        snakemake.params.af_germline_lower_limit,
        snakemake.params.af_germline_upper_limit,
        snakemake.params.gnomad_limit,
        snakemake.params.db1000g_limit,
        snakemake.params.nr_avg_germline_snvs,
        snakemake.params.nssnv_tmb_correction,
    )
