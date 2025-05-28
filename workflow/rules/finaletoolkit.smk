__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2025, Jonas Almlöf"
__email__ = "jonas.almlöf@scilifelab.uu.se"
__license__ = "GPL-3"


rule finaletoolkit_end_motifs:
    input:
        bam="alignment/bwa_mem_realign_consensus_reads/{sample}_{type}.umi.bam",
        bai="alignment/bwa_mem_realign_consensus_reads/{sample}_{type}.umi.bam.bai",
        reference=config.get("reference", {}).get("fasta_2bit", ""),
    output:
        end_motifs="biomarker/finaletoolkit_end_motifs/{sample}_{type}.end-motifs.tsv",
    params:
        extra=config.get("finaletoolkit_end_motifs", {}).get("extra", ""),
    log:
        "biomarker/finaletoolkit_end_motifs/{sample}_{type}.end-motifs.tsv.log",
    benchmark:
        repeat(
            "biomarker/finaletoolkit_end_motifs/{sample}_{type}.end-motifs.tsv.benchmark.tsv",
            config.get("finaletoolkit_end_motifs", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("finaletoolkit_end_motifs", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("finaletoolkit_end_motifs", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("finaletoolkit_end_motifs", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("finaletoolkit_end_motifs", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("finaletoolkit_end_motifs", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("finaletoolkit_end_motifs", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("finaletoolkit_end_motifs", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculate end motif frequencies in plasma samples into {output.end_motifs}"
    shell:
        "finaletoolkit end-motifs "
        "{input.bam} "
        "{input.reference} "
        "-o {output.end_motifs} "
        "{params.extra} >& {log}"
