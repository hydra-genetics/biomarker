
# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule optitype:
    input:
        fastq1="prealignment/merged/{sample}_{type}_fastq1.fastq.gz",
        fastq2="prealignment/merged/{sample}_{type}_fastq2.fastq.gz",
    output:
        coverage_plot=temp("biomarker/optitype/{sample}_{type}/{sample}_{type}_coverage_plot.pdf"),
        hla_type=temp("biomarker/optitype/{sample}_{type}/{sample}_{type}_result.tsv"),
        out_dir=temp(directory("biomarker/optitype/{sample}_{type}/")),
    params:
        extra=config.get("optitype", {}).get("extra", ""),
        out_prefix="{sample}_{type}",
        sample_type=config.get("optitype", {}).get("sample_type", "-d"),
        enumeration=config.get("optitype", {}).get("enumeration", "4"),
    log:
        "biomarker/optitype/{sample}_{type}/{sample}_{type}_result.tsv.log",
    benchmark:
        repeat(
            "biomarker/optitype/{sample}_{type}/{sample}_{type}_result.tsv.benchmark.tsv",
            config.get("optitype", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("optitype", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("optitype", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("optitype", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("optitype", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("optitype", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("optitype", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("optitype", {}).get("container", config["default_container"])
    conda:
        "../envs/optitype.yaml"
    message:
        "{rule}: Determine HLA-type in biomarker/{rule}/{wildcards.sample}_{wildcards.type}/{wildcards.sample}_{wildcards.type}"
    shell:
        "(python /usr/local/bin/OptiType/OptiTypePipeline.py "
        "-i {input.fastq1} {intput.fastq2} "
        "{params.sample_type} "
        "--enumerate {params.enumerate} "
        "-p {params.out_prefix} "
        "-o {params.out_dir}) &> {log}"
