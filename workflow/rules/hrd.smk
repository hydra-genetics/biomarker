__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule hrd:
    input:
        segment="cnv_sv/cnvkit_call/{sample}_{type}.loh.cns",
    output:
        hrd=temp("biomarker/hrd/{sample}_{type}.hrd_score.txt"),
    log:
        "biomarker/hrd/{sample}_{type}.hrd_score.txt.log",
    benchmark:
        repeat("biomarker/hrd/{sample}_{type}.hrd_score.txt.benchmark.tsv", config.get("hrd", {}).get("benchmark_repeats", 1))
    threads: config.get("hrd", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("hrd", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("hrd", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("hrd", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("hrd", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("hrd", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("hrd", {}).get("container", config["default_container"])
    conda:
        "../envs/hrd.yaml"
    message:
        "{rule}: calculate hrd in hrd/{rule}/{wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/hrd.py"
