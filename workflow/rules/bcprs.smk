__author__ = "Per Lundmark"
__copyright__ = "Copyright 2026, Per Lundmark"
__email__ = "per.lundmark@scilifelab.uu.se"
__license__ = "GPL-3"


rule bcprs:
    input:
        vcf="snv_indels/deepvariant/{sample}_{type}.g.fix_af.bcprs.vcf.gz",
        score_file=config.get("bcprs", {}).get("score_file", ""),
    output:
        detailed=temp("biomarker/bcprs/{sample}_{type}.detailed_scores.csv"),
        total=temp("biomarker/bcprs/{sample}_{type}.total_scores.csv"),
    params:
        no_proxy_dedup=config.get("bcprs", {}).get("no_proxy_dedup", False),
        r_compat=config.get("bcprs", {}).get("r_compat", False),
    log:
        "biomarker/bcprs/{sample}_{type}.log",
    benchmark:
        repeat("biomarker/bcprs/{sample}_{type}.benchmark.tsv", config.get("bcprs", {}).get("benchmark_repeats", 1))
    threads: config.get("bcprs", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bcprs", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcprs", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcprs", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcprs", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bcprs", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bcprs", {}).get("container", config["default_container"])
    message:
        "{rule}: calculate polygenic risk score for {input.vcf} using {input.score_file}"
    script:
        "../scripts/bcprs.py"
