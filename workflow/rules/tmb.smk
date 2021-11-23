# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule tmb:
    input:
        vcf="snv_indels/ensemble_vcf/{sample}_{type}.ensembled.vep_annotated.vcf",
        artifacts=config["reference"]["artifacts"],
        background_panel=config["reference"]["background"],
        background_run=lambda wildcards: "annotation/calculate_seqrun_background/%s_seqrun_background.tsv"
        % get_run(units, wildcards),
        gvcf="snv_indels/mutect2_gvcf/{sample}_{type}.merged.gvcf.gz",
    output:
        tmb=temp("biomarker/tmb/{sample}_{type}.TMB.txt"),
    log:
        "biomarker/tmb/{sample}_{type}.log",
    benchmark:
        repeat("biomarker/tmb/{sample}_{type}.benchmark.tsv", config.get("tmb", {}).get("benchmark_repeats", 1))
    threads: config.get("tmb", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("tmb", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("tmb", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("tmb", {}).get("container", config["default_container"])
    conda:
        "../envs/tmb.yaml"
    message:
        "{rule}: Calculate TMB in tmb/{rule}/{wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/tmb.py"
