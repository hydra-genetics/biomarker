__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2026, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule fragle:
    input:
        bam=lambda wildcards: get_input_aligned_bam(wildcards, config)[0],
        bai=lambda wildcards: get_input_aligned_bam(wildcards, config)[1],
    output:
        csv=temp("biomarker/fragle/{sample}_{type}/Fragle.csv"),
        data=temp("biomarker/fragle/{sample}_{type}/data.pkl"),
        ht=temp("biomarker/fragle/{sample}_{type}/HT.csv"),
        lt=temp("biomarker/fragle/{sample}_{type}/LT.csv"),
        off_target_bam=temp("biomarker/fragle/{sample}_{type}/off_target_bams/{sample}_{type}.bam"),
        off_target_bai=temp("biomarker/fragle/{sample}_{type}/off_target_bams/{sample}_{type}.bam.bai"),
    params:
        outdir=lambda wildcards, output: os.path.dirname(os.path.dirname(output.csv)),
        design_bed=config.get("fragle", {}).get("design_bed", ""),
        genome_build=config.get("fragle", {}).get("genome_build", "hg19"),
        mode=config.get("fragle", {}).get("mode", "T"),
        threads_calc=lambda wildcards, threads: max(1, threads // 2),
    log:
        "biomarker/fragle/{sample}_{type}.log",
    benchmark:
        repeat(
            "biomarker/fragle/{sample}_{type}.benchmark.tsv",
            config.get("fragle", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fragle", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fragle", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fragle", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("fragle", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("fragle", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("fragle", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fragle", {}).get("container", config["default_container"])
    message:
        "{rule}: Run FRAGLE ctDNA fraction estimation on {input.bam}"
    shell:
        "(python /opt/Fragle_app/main.py "
        "--input {input.bam} "
        "--output {params.outdir} "
        "--mode {params.mode} "
        "--genome_build {params.genome_build} "
        "--target_bed {params.design_bed} "
        "--cpu {params.threads_calc} "
        "--threads {params.threads_calc}) &> {log}"
