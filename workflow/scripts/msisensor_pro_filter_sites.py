def read_msi_sites(msi_sites_filename):
    msi_sites = open(msi_sites_filename)
    msi_sites_dict = {}
    for line in msi_sites:
        columns = line.strip().split("\t")
        chrom = columns[0]
        start_pos = int(columns[1])
        end_pos = int(columns[2])
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
        pos = start_pos
        while pos <= end_pos:
            chrom_pos = f"{chrom}_{pos}"
            msi_sites_dict[chrom_pos] = ""
            pos += 1
    return msi_sites_dict


def write_msi_sites(in_PoN, out_PoN, msi_sites_dict):
    header = True
    for line in in_PoN:
        if header:
            out_PoN.write(line)
            header = False
            continue
        columns = line.strip().split("\t")
        chrom = columns[0]
        pos = columns[1]
        chrom_pos = f"{chrom}_{pos}"
        if chrom_pos in msi_sites_dict or msi_sites_dict == {}:
            out_PoN.write(line)


def filter_sites(in_PoN, out_PoN, msi_sites_filename):
    msi_sites_dict = {}
    if msi_sites_filename != "":
        msi_sites_dict = read_msi_sites(msi_sites_filename)
    write_msi_sites(in_PoN, out_PoN, msi_sites_dict)


if __name__ == "__main__":
    filter_sites(
        open(snakemake.input.PoN),
        open(snakemake.output.PoN, "w"),
        snakemake.params.msi_sites_bed,
    )
