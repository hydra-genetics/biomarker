__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2025, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule finaletoolkit_end_motifs:
    input:
        bam="alignment/bwa_mem_realign_consensus_reads/{sample}_{type}.umi.bam",
        bai="alignment/bwa_mem_realign_consensus_reads/{sample}_{type}.umi.bam.bai",
        reference=config.get("reference", {}).get("fasta_2bit", ""),
    output:
        end_motifs=temp("biomarker/finaletoolkit_end_motifs/{sample}_{type}.end-motifs.tsv"),
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
        "{params.extra} "
        "-o {output.end_motifs} "
        "{params.extra} >& {log}"


rule finaletoolkit_interval_end_motifs:
    input:
        bam="alignment/bwa_mem_realign_consensus_reads/{sample}_{type}.umi.bam",
        bai="alignment/bwa_mem_realign_consensus_reads/{sample}_{type}.umi.bam.bai",
        design_bed=config.get("reference", {}).get("design_bed", ""),
        reference=config.get("reference", {}).get("fasta_2bit", ""),
    output:
        interval_end_motifs=temp("biomarker/finaletoolkit_interval_end_motifs/{sample}_{type}.interval-end-motifs.tsv"),
    params:
        extra=config.get("finaletoolkit_interval_end_motifs", {}).get("extra", ""),
    log:
        "biomarker/finaletoolkit_interval_end_motifs/{sample}_{type}.interval-end-motifs.tsv.log",
    benchmark:
        repeat(
            "biomarker/finaletoolkit_interval_end_motifs/{sample}_{type}.interval-end-motifs.tsv.benchmark.tsv",
            config.get("finaletoolkit_interval_end_motifs", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("finaletoolkit_interval_end_motifs", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("finaletoolkit_interval_end_motifs", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("finaletoolkit_interval_end_motifs", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("finaletoolkit_interval_end_motifs", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("finaletoolkit_interval_end_motifs", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("finaletoolkit_interval_end_motifs", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("finaletoolkit_interval_end_motifs", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculate end motif frequencies for each region in bed file in plasma samples into {output.interval_end_motifs}"
    shell:
        "finaletoolkit interval-end-motifs "
        "{input.bam} "
        "{input.reference} "
        "{input.design_bed} "
        "{params.extra} "
        "-o {output.interval_end_motifs} > {log}"


rule finaletoolkit_mds:
    input:
        end_motifs="biomarker/finaletoolkit_end_motifs/{sample}_{type}.end-motifs.tsv",
    output:
        mds=temp("biomarker/finaletoolkit_mds/{sample}_{type}.mds.txt"),
    log:
        "biomarker/finaletoolkit_mds/{sample}_{type}.mds.txt.log",
    benchmark:
        repeat(
            "biomarker/finaletoolkit_mds/{sample}_{type}.mds.txt.benchmark.tsv",
            config.get("finaletoolkit_mds", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("finaletoolkit_mds", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("finaletoolkit_mds", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("finaletoolkit_mds", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("finaletoolkit_mds", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("finaletoolkit_mds", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("finaletoolkit_mds", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("finaletoolkit_mds", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculate mds based on end motif frequencies in plasma samples into {output.mds}"
    shell:
        "finaletoolkit mds "
        "{input.end_motifs} "
        "> {output.mds} 2> {log}"


rule finaletoolkit_interval_mds:
    input:
        interval_end_motifs="biomarker/finaletoolkit_interval_end_motifs/{sample}_{type}.interval-end-motifs.tsv",
    output:
        interval_mds=temp("biomarker/finaletoolkit_interval_mds/{sample}_{type}.interval-mds.txt"),
    log:
        "biomarker/finaletoolkit_interval_mds/{sample}_{type}.interval-mds.txt.log",
    benchmark:
        repeat(
            "biomarker/finaletoolkit_interval_mds/{sample}_{type}.interval-mds.txt.benchmark.tsv",
            config.get("finaletoolkit_interval_mds", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("finaletoolkit_interval_mds", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("finaletoolkit_interval_mds", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("finaletoolkit_interval_mds", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("finaletoolkit_interval_mds", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("finaletoolkit_interval_mds", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("finaletoolkit_interval_mds", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("finaletoolkit_interval_mds", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculate mds based on interval end motif frequencies in plasma samples into {output.interval_mds}"
    shell:
        "finaletoolkit interval-mds "
        "{input.interval_end_motifs} "
        "{output.interval_mds} > {log}"


rule finaletoolkit_frag_length_bins:
    input:
        bam="alignment/bwa_mem_realign_consensus_reads/{sample}_{type}.umi.bam",
        bai="alignment/bwa_mem_realign_consensus_reads/{sample}_{type}.umi.bam.bai",
    output:
        frag_len_bins=temp("biomarker/finaletoolkit_frag_length_bins/{sample}_{type}.frag-length-bins.tsv"),
    params:
        extra=config.get("finaletoolkit_frag_length_bins", {}).get("extra", ""),
    log:
        "biomarker/finaletoolkit_frag_length_bins/{sample}_{type}.frag-length-bins.tsv.log",
    benchmark:
        repeat(
            "biomarker/finaletoolkit_frag_length_bins/{sample}_{type}.frag-length-bins.tsv.benchmark.tsv",
            config.get("finaletoolkit_frag_length_bins", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("finaletoolkit_frag_length_bins", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("finaletoolkit_frag_length_bins", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("finaletoolkit_frag_length_bins", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("finaletoolkit_frag_length_bins", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("finaletoolkit_frag_length_bins", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("finaletoolkit_frag_length_bins", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("finaletoolkit_frag_length_bins", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculate histogram with fragment lengths into {output.frag_len_bins}"
    shell:
        "finaletoolkit frag-length-bins "
        "-o {output.frag_len_bins} "
        "{params.extra} "
        "{input.bam} >& {log}"
