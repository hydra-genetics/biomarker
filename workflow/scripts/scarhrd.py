
import subprocess

seq = snakemake.input.seg
reference_name = snakemake.params.reference
seqz = snakemake.params.seqz



cmd = f"Rscript -e 'scarHRD::scar_score(\"{seg}\", reference=\"{reference_name}\", seqz={seqzIn})'"
print(cmd)
subprocess.run(cmd, shell=True)
