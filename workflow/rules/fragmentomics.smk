__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2025, Jonas Almlöf"
__email__ = "jonas.scilifelab.uu.se"
__license__ = "GPL-3"


rule fragmentomics_fragment_length_patient_score:
    input:
        bam="alignment/bwa_mem_realign_consensus_reads/{sample}_{type}.umi.bam",
        bai="alignment/bwa_mem_realign_consensus_reads/{sample}_{type}.umi.bam.bai",
        reference_set=config.get("fragmentomics_fragment_length_patient_score", {}).get("reference_set", ""),
    output:
        patient_score=temp("biomarker/fragmentomics_fragment_length_patient_score/{sample}_{type}.fragment_length_patient_score.txt"),
    log:
        "biomarker/fragmentomics_fragment_length_patient_score/{sample}_{type}.fragment_length_patient_score.txt.log",
    benchmark:
        repeat(
            "biomarker/fragmentomics_fragment_length_patient_score/{sample}_{type}.fragment_length_patient_score.txt.benchmark.tsv",
            config.get("fragmentomics_fragment_length_patient_score", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fragmentomics_fragment_length_patient_score", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fragmentomics_fragment_length_patient_score", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fragmentomics_fragment_length_patient_score", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("fragmentomics_fragment_length_patient_score", {}).get(
            "partition", config["default_resources"]["partition"]
        ),
        threads=config.get("fragmentomics_fragment_length_patient_score", {}).get(
            "threads", config["default_resources"]["threads"]
        ),
        time=config.get("fragmentomics_fragment_length_patient_score", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fragmentomics_fragment_length_patient_score", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates a fragment length score and store in {output.patient_score}"
    script:
        "../scripts/patient_score.R"