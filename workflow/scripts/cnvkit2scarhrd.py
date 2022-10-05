
cnvkitseg_filename = snakemake.input.seg
cnvkitseg = open(cnvkitseg_filename)
scrahrdseg = open(snakemake.output.seg, "w")


scrahrdseg.write("SampleID\tChromosome\tStart_position\tEnd_position\ttotal_cn\tA_cn\tB_cn\tploidy\n")

sample_name = cnvkitseg_filename.split("/")[-1].split(".")[0]

header = True
header_dict = {}
for line in cnvkitseg:
    columns = line.strip().split("\t")
    if header:
        i = 0
        for c in columns:
            header_dict[c] = i
            i += 1
        header = False
        continue
    chrom = columns[header_dict["chromosome"]]
    start_pos = columns[header_dict["start"]]
    ens_pos = columns[header_dict["end"]]
    total_cn = columns[header_dict["cn"]]
    A_cn = columns[header_dict["cn1"]]
    B_cn = columns[header_dict["cn2"]]
    ploidy = "NA"
    # ignore Y chrom with 0 copies, probably female ##
    if chrom == "chrY" and total_cn == 0:
        continue
    if A_cn == "" or B_cn == "":
        A_cn = total_cn - 1
        B_cn = 1
    scrahrdseg.write(f"{sample}\t{chrom}\t{start_pos}\t{end_pos}\t{total_cn}\t{A_cn}\t{B_cn}\t{ploidy}\n")
