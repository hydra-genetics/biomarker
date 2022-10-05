
import subprocess
import os

seg = snakemake.input.seg
reference_name = snakemake.params.reference_name
seqz = snakemake.params.seqz
output_name_fixed = snakemake.output.hrd
outputdir = os.path.dirname(output_name_fixed)
output = open(output_name_fixed, "w")


cmd = f"Rscript -e 'scarHRD::scar_score(\"{seg}\", reference=\"{reference_name}\", seqz=FALSE, outputdir=\"{outputdir}\")'"
print(cmd)
subprocess.run(cmd, shell=True)


sample_name = seg.split("/")[-1].split(".scarhrd")[0]
R_output_name = outputdir + "/" + sample_name + "_HRDresults.txt"
R_input = open(R_output_name)

output.write("HRD-score\tHRD\tTelomeric_AI\tLST\n")
header = True
for line in R_input:
    if header:
        header = False
        continue
    columns = line.strip().split("\t")
    output.write(f"{columns[4]}\t{columns[1]}\t{columns[2]}\t{columns[3]}\n")
output.close()
R_input.close()
