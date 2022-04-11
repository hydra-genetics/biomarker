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
        coverage_plot=temp("biomarker/optitype/{sample}_{type}/{sample}_{type}_hla_type_coverage_plot.pdf"),
        hla_type=temp("biomarker/optitype/{sample}_{type}/{sample}_{type}_hla_type_result.tsv"),
        out_dir=temp(directory("biomarker/optitype/{sample}_{type}/")),
    params:
        extra=config.get("optitype", {}).get("extra", ""),
        out_prefix="{sample}_{type}_hla_type",
        sample_type=config.get("optitype", {}).get("sample_type", "-d"),
        enumeration=config.get("optitype", {}).get("enumeration", "4"),
    log:
        "biomarker/optitype/{sample}_{type}/{sample}_{type}_result.tsv.log",
    benchmark:
        repeat(
            "biomarker/optitype/{sample}_{type}/{sample}_{type}_result.tsv.benchmark.tsv",
            config.get("optitype", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("optitype", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("optitype", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("optitype", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("optitype", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("optitype", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("optitype", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("optitype", {}).get("container", config["default_container"])
    conda:
        "../envs/optitype.yaml"
    message:
        "{rule}: Determine HLA-type in {output.hla_type}"
    shell:
        "(OptiTypePipeline.py "
        "-i {input.fastq1} {input.fastq2} "
        "{params.sample_type} "
        "--enumerate {params.enumeration} "
        "-p {params.out_prefix} "
        "-o {output.out_dir}) &> {log}"
