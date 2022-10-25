__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2022, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnvkit2scarhrd:
    input:
        seg="cnv_sv/cnvkit_call/{sample}_{type}.loh.cns",
    output:
        seg=temp("biomarker/cnvkit2scarhrd/{sample}_{type}.scarhrd.cns"),
    log:
        "biomarker/cnvkit2scarhrd/{sample}_{type}.scarhrd.cns.log",
    benchmark:
        repeat(
            "biomarker/cnvkit2scarhrd/{sample}_{type}.scarhrd.cns.benchmark.tsv",
            config.get("cnvkit2scarhrd", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit2scarhrd", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnvkit2scarhrd", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit2scarhrd", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit2scarhrd", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnvkit2scarhrd", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit2scarhrd", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnvkit2scarhrd", {}).get("container", config["default_container"])
    conda:
        "../envs/scarhrd.yaml"
    message:
        "{rule}: convert cnvkit segmentation files to scarhrd input: {output.seg}"
    script:
        "../scripts/cnvkit2scarhrd.py"


rule scarhrd:
    input:
        seg_cnvkit="biomarker/cnvkit2scarhrd/{sample}_{type}.scarhrd.cns",
    output:
        hrd="biomarker/scarhrd/{sample}_{type}/{sample}_{type}_HRDresults.txt",
        hrd_dir=directory("biomarker/scarhrd/{sample}_{type}/"),
    params:
        reference_name=config.get("scarhrd", {}).get("reference_name", "grch37"),
        seqz=config.get("scarhrd", {}).get("seqz", False),
    log:
        "biomarker/scarhrd/{sample}_{type}_HRDresults.txt.log",
    benchmark:
        repeat(
            "biomarker/scarhrd/{sample}_{type}_HRDresults.txt.benchmark.tsv",
            config.get("scarhrd", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("scarhrd", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("scarhrd", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("scarhrd", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("scarhrd", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("scarhrd", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("scarhrd", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("scarhrd", {}).get("container", config["default_container"])
    conda:
        "../envs/scarhrd.yaml"
    message:
        "{rule}: calculate hrd on {input.seg_cnvkit}"
    shell:
        "(Rscript -e 'scarHRD::scar_score( "
        '"{input.seg_cnvkit}", '
        'reference="{params.reference_name}", '
        'seqz="{params.seqz}", '
        'outputdir="{output.hrd_dir}")\') &> {log}'


rule fix_scarhrd_output:
    input:
        hrd="biomarker/scarhrd/{sample}_{type}/{sample}_{type}_HRDresults.txt",
    output:
        hrd="biomarker/scarhrd/{sample}_{type}.scarhrd_cnvkit_score.txt",
    log:
        "biomarker/scarhrd/{sample}_{type}.scarhrd_score.txt.log",
    benchmark:
        repeat(
            "biomarker/scarhrd/{sample}_{type}.scarhrd_score.txt.benchmark.tsv",
            config.get("fix_scarhrd_output", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fix_scarhrd_output", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fix_scarhrd_output", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fix_scarhrd_output", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("fix_scarhrd_output", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("fix_scarhrd_output", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("fix_scarhrd_output", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fix_scarhrd_output", {}).get("container", config["default_container"])
    conda:
        "../envs/scarhrd.yaml"
    message:
        "{rule}: fix scarhrd output into {output.hrd}"
    script:
        "../scripts/fix_scarhrd_output.py"
