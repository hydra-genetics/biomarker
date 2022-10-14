
import subprocess
import os

seg_cnvkit = snakemake.input.seg_cnvkit
reference_name = snakemake.params.reference_name
seqz = snakemake.params.seqz
output_name_fixed_cnvkit = snakemake.output.hrd_cnvkit
outputdir_cnvkit = os.path.dirname(output_name_fixed_cnvkit)
output_cnvkit = open(output_name_fixed_cnvkit, "w")

# scarHRD on cnvkit
cmd = f"Rscript -e 'scarHRD::scar_score(\"{seg_cnvkit}\", reference=\"{reference_name}\", seqz=\"{seqz}\", \
    outputdir=\"{outputdir_cnvkit}\")'"
print(cmd)
subprocess.run(cmd, shell=True)


sample_name = seg_cnvkit.split("/")[-1].split(".scarhrd")[0]
R_output_name = outputdir_cnvkit + "/" + sample_name + "_HRDresults.txt"
R_input = open(R_output_name)

output_cnvkit.write("HRD-score\tHRD\tTelomeric_AI\tLST\n")
header = True
for line in R_input:
    if header:
        header = False
        continue
    columns = line.strip().split("\t")
    output_cnvkit.write(f"{columns[4]}\t{columns[1]}\t{columns[2]}\t{columns[3]}\n")
output_cnvkit.close()
R_input.close()
