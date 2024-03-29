__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule msisensor_pro:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        PoN=config.get("msisensor_pro", {}).get("PoN", ""),
    output:
        msi_score=temp("biomarker/msisensor_pro/{sample}_{type}"),
        msi_all=temp("biomarker/msisensor_pro/{sample}_{type}_all"),
        msi_dis=temp("biomarker/msisensor_pro/{sample}_{type}_dis"),
        msi_unstable=temp("biomarker/msisensor_pro/{sample}_{type}_unstable"),
    params:
        extra=config.get("msisensor_pro", {}).get("extra", "-c 50 -b 2"),
        out_prefix="biomarker/msisensor_pro/{sample}_{type}",
    log:
        "biomarker/msisensor_pro/{sample}_{type}.log",
    benchmark:
        repeat(
            "biomarker/msisensor_pro/{sample}_{type}.benchmark.tsv", config.get("msisensor_pro", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("msisensor_pro", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("msisensor_pro", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("msisensor_pro", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("msisensor_pro", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("msisensor_pro", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("msisensor_pro", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("msisensor_pro", {}).get("container", config["default_container"])
    message:
        "{rule}: calculate MSI-status in msisensor_pro/{rule}/{wildcards.sample}_{wildcards.type}"
    shell:
        "(msisensor-pro pro {params.extra} -d {input.PoN} -t {input.bam} -o {params.out_prefix}) &> {log}"
