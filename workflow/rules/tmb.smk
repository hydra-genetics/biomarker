__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule tmb:
    input:
        vcf="annotation/background_annotation/{sample}_{type}.background_annotation.vcf.gz",
    output:
        tmb=temp("biomarker/tmb/{sample}_{type}.TMB.txt"),
    params:
        dp_limit=config.get("tmb", {}).get("dp_limit", 200),
        vd_limit=config.get("tmb", {}).get("vd_limit", 10),
        af_lower_limit=config.get("tmb", {}).get("af_lower_limit", 0.05),
        af_upper_limit=config.get("tmb", {}).get("af_upper_limit", 0.95),
        af_germline_lower_limit=config.get("tmb", {}).get("af_germline_lower_limit", 0.47),
        af_germline_upper_limit=config.get("tmb", {}).get("af_germline_upper_limit", 0.53),
        gnomad_limit=config.get("tmb", {}).get("gnomad_limit", 0.0001),
        db1000g_limit=config.get("tmb", {}).get("db1000g_limit", 0.0001),
        nr_avg_germline_snvs=config.get("tmb", {}).get("nr_avg_germline_snvs", 2.0),
        nssnv_tmb_correction=config.get("tmb", {}).get("nssnv_tmb_correction", 0.84),
    log:
        "biomarker/tmb/{sample}_{type}.TMB.txt.log",
    benchmark:
        repeat("biomarker/tmb/{sample}_{type}.TMB.txt.benchmark.tsv", config.get("tmb", {}).get("benchmark_repeats", 1))
    threads: config.get("tmb", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("tmb", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("tmb", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("tmb", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("tmb", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("tmb", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("tmb", {}).get("container", config["default_container"])
    conda:
        "../envs/tmb.yaml"
    message:
        "{rule}: calculate TMB in tmb/{rule}/{wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/tmb.py"
