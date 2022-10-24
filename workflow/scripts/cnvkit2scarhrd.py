

def cnvkit_2_scarhrd(cnvkitseg_filename, scrahrdseg):

    cnvkitseg = open(cnvkitseg_filename)
    scrahrdseg.write("SampleID\tChromosome\tStart_position\tEnd_position\ttotal_cn\tA_cn\tB_cn\tploidy\n")
    sample_name = cnvkitseg_filename.split("/")[-1].split(".loh")[0]

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
        end_pos = columns[header_dict["end"]]
        total_cn = int(columns[header_dict["cn"]])
        A_cn = columns[header_dict["cn1"]]
        B_cn = columns[header_dict["cn2"]]
        ploidy = "NA"
        # ignore Y chrom with 0 copies, probably female ##
        if chrom == "chrY" and total_cn == 0:
            continue
        if A_cn == "" or B_cn == "":
            A_cn = total_cn - 1
            B_cn = 1
        scrahrdseg.write(f"{sample_name}\t{chrom}\t{start_pos}\t{end_pos}\t{total_cn}\t{A_cn}\t{B_cn}\t{ploidy}\n")


if __name__ == "__main__":
    log = snakemake.log_fmt_shell(stdout=False, stderr=True)

    cnvkit_2_scarhrd(
        snakemake.input.seg,
        open(snakemake.output.seg, "w")
    )
