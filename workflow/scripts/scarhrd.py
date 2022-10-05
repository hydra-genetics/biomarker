
import subprocess

seg = snakemake.input.seg
reference_name = snakemake.params.reference_name
seqz = snakemake.params.seqz


cmd = f"Rscript -e 'scarHRD::scar_score(\"{seg}\", reference=\"{reference_name}\", seqz=FALSE)'"
print(cmd)
subprocess.run(cmd, shell=True)
