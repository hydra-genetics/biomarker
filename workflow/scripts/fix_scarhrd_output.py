
import subprocess
import os


def fix_scarhrd(in_hrd, out_hrd):

    out_hrd.write("HRD-score\tHRD\tTelomeric_AI\tLST\n")
    header = True
    for line in in_hrd:
        if header:
            header = False
            continue
        columns = line.strip().split("\t")
        out_hrd.write(f"{columns[4]}\t{columns[1]}\t{columns[2]}\t{columns[3]}\n")
    out_hrd.close()
    in_hrd.close()


if __name__ == "__main__":
    log = snakemake.log_fmt_shell(stdout=False, stderr=True)

    fix_scarhrd(
        open(snakemake.input.hrd),
        open(snakemake.output.hrd, "w"),
    )
