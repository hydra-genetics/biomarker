
gatk_cnv_seg_filename = snakemake.input.seg
gatk_cnv_seg = open(gatk_cnv_seg_filename)
scrahrdseg = open(snakemake.output.seg, "w")


scrahrdseg.write("SampleID\tChromosome\tStart_position\tEnd_position\ttotal_cn\tA_cn\tB_cn\tploidy\n")

sample_name = gatk_cnv_seg_filename.split("/")[-1].split(".loh")[0]

header = True
for line in gatk_cnv_seg:
    columns = line.strip().split("\t")
    if header:
        if columns[0] != "#CHROM":
            continue
        header = False
        continue
    chrom = columns[0]
    start_pos = columns[1]
    end_pos = columns[7].split(";END=")[1].split(";")[0]
    total_cn_exact = float(columns[7].split(";CORR_CN=")[1].split(";")[0])
    total_cn = 2
    if total_cn_exact > 0.3:
        total_cn = int(round(total_cn_exact, 0))
    elif total_cn_exact < -0.4:
        total_cn = int(round(total_cn_exact, 0))
    A_cn = total_cn - 1
    B_cn = 1
    ploidy = "NA"
    # ignore Y chrom with 0 copies, probably female ##
    if chrom == "chrY" and total_cn == 0:
        continue
    if A_cn == "" or B_cn == "":
        A_cn = total_cn - 1
        B_cn = 1
    scrahrdseg.write(f"{sample_name}\t{chrom}\t{start_pos}\t{end_pos}\t{total_cn}\t{A_cn}\t{B_cn}\t{ploidy}\n")
