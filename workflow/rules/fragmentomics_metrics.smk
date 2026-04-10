# Fragmentomics Metrics rules for Hydra Genetics Biomarker module


rule fragmentomics_metrics_get_bed_from_bam:
    input:
        bam="alignment/bwa_mem_realign_consensus_reads/{sample}_{type}.umi.bam",
        bai="alignment/bwa_mem_realign_consensus_reads/{sample}_{type}.umi.bam.bai",
    output:
        bed_gz=temp("biomarker/fragmentomics_metrics_get_bed_from_bam/{sample}_{type}.bed.gz"),
    params:
        canonical_cds_bed=config.get("fragmentomics_metrics", {}).get(
            "canonical_cds_bed", "resources/UCSC_hg19_canonical_cds.bed"
        ),
        min_mapq=config.get("fragmentomics_metrics", {}).get("min_mapq", 30),
        sort_mem=lambda wildcards, resources: f"{int(resources.mem_mb * 0.8)}M" if getattr(resources, "mem_mb", None) else "8G",
    log:
        "biomarker/fragmentomics_metrics_get_bed_from_bam/{sample}_{type}.log",
    benchmark:
        repeat(
            "biomarker/fragmentomics_metrics_get_bed_from_bam/{sample}_{type}.benchmark.tsv",
            config.get("fragmentomics_metrics_get_bed_from_bam", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fragmentomics_metrics_get_bed_from_bam", {}).get("threads", 6)
    resources:
        mem_mb=config.get("fragmentomics_metrics_get_bed_from_bam", {}).get("mem_mb", 30000),
        mem_per_cpu=config.get("fragmentomics_metrics_get_bed_from_bam", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("fragmentomics_metrics_get_bed_from_bam", {}).get(
            "partition", config["default_resources"]["partition"]
        ),
        threads=config.get("fragmentomics_metrics_get_bed_from_bam", {}).get("threads", 6),
        time=config.get("fragmentomics_metrics_get_bed_from_bam", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fragmentomics_metrics", {}).get("container", config["default_container"])
    message:
        "{rule}: Extracts BED from BAM for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        set -o pipefail
        (samtools view -@ {threads} -b -q {params.min_mapq} {input.bam} \
        | bedtools bamtobed -i stdin \
        | cut -f 1-6 \
        | bedtools intersect -wa -wb -a stdin -b {params.canonical_cds_bed} \
        | awk '{{ if (NF < 14) {{ print "ERROR: Expected at least 14 fields from intersect, got " NF > "/dev/stderr"; exit 1; }} print }}' \
        | cut -f 1,2,3,10,11,12,13,14 \
        | sort --parallel={threads} -S {params.sort_mem} -k1,1 -k2,2n \
        | gzip > {output}) > {log} 2>&1
        """


rule fragmentomics_metrics_filter_comm_panel_genes:
    input:
        bed_gz="biomarker/fragmentomics_metrics_get_bed_from_bam/{sample}_{type}.bed.gz",
    output:
        bed=temp("biomarker/fragmentomics_metrics_filter_comm_panel_genes/{sample}_{type}.filtered.bed"),
    log:
        "biomarker/fragmentomics_metrics_filter_comm_panel_genes/{sample}_{type}.log",
    benchmark:
        repeat(
            "biomarker/fragmentomics_metrics_filter_comm_panel_genes/{sample}_{type}.benchmark.tsv",
            config.get("fragmentomics_metrics_filter_comm_panel_genes", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fragmentomics_metrics_filter_comm_panel_genes", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fragmentomics_metrics_filter_comm_panel_genes", {}).get(
            "mem_mb", config["default_resources"]["mem_mb"]
        ),
        mem_per_cpu=config.get("fragmentomics_metrics_filter_comm_panel_genes", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("fragmentomics_metrics_filter_comm_panel_genes", {}).get(
            "partition", config["default_resources"]["partition"]
        ),
        threads=config.get("fragmentomics_metrics_filter_comm_panel_genes", {}).get(
            "threads", config["default_resources"]["threads"]
        ),
        time=config.get("fragmentomics_metrics_filter_comm_panel_genes", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fragmentomics_metrics", {}).get("container", config["default_container"])
    message:
        "{rule}: Filters commercial panel genes for {wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/fragmentomics_metrics/filter_comm_panels.R"


rule fragmentomics_metrics_gzip_comm_panel_genes:
    input:
        bed="biomarker/fragmentomics_metrics_filter_comm_panel_genes/{sample}_{type}.filtered.bed",
    output:
        bed_gz_filtered=temp("biomarker/fragmentomics_metrics_filter_comm_panel_genes/{sample}_{type}.filtered.bed.gz"),
    log:
        "biomarker/fragmentomics_metrics_gzip_comm_panel_genes/{sample}_{type}.log",
    benchmark:
        repeat(
            "biomarker/fragmentomics_metrics_gzip_comm_panel_genes/{sample}_{type}.benchmark.tsv",
            config.get("fragmentomics_metrics_gzip_comm_panel_genes", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fragmentomics_metrics_gzip_comm_panel_genes", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fragmentomics_metrics_gzip_comm_panel_genes", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fragmentomics_metrics_gzip_comm_panel_genes", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("fragmentomics_metrics_gzip_comm_panel_genes", {}).get(
            "partition", config["default_resources"]["partition"]
        ),
        threads=config.get("fragmentomics_metrics_gzip_comm_panel_genes", {}).get(
            "threads", config["default_resources"]["threads"]
        ),
        time=config.get("fragmentomics_metrics_gzip_comm_panel_genes", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fragmentomics_metrics", {}).get("container", config["default_container"])
    message:
        "{rule}: Gzips filtered commercial panel genes for {wildcards.sample}_{wildcards.type}"
    shell:
        "gzip -c {input.bed} > {output.bed_gz_filtered} 2> {log}"


rule fragmentomics_metrics_get_SE_fragstats:
    input:
        bed_gz_filtered="biomarker/fragmentomics_metrics_filter_comm_panel_genes/{sample}_{type}.filtered.bed.gz",
    output:
        se_fragstats=temp("biomarker/fragmentomics_metrics_get_SE_fragstats/{sample}_{type}.txt.gz"),
    log:
        "biomarker/fragmentomics_metrics_get_SE_fragstats/{sample}_{type}.log",
    benchmark:
        repeat(
            "biomarker/fragmentomics_metrics_get_SE_fragstats/{sample}_{type}.benchmark.tsv",
            config.get("fragmentomics_metrics_get_SE_fragstats", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fragmentomics_metrics_get_SE_fragstats", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fragmentomics_metrics_get_SE_fragstats", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fragmentomics_metrics_get_SE_fragstats", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("fragmentomics_metrics_get_SE_fragstats", {}).get(
            "partition", config["default_resources"]["partition"]
        ),
        threads=config.get("fragmentomics_metrics_get_SE_fragstats", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("fragmentomics_metrics_get_SE_fragstats", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fragmentomics_metrics", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates SE fragstats for {wildcards.sample}_{wildcards.type}"
    priority: 20
    shell:
        """
        (zcat {input.bed_gz_filtered} \
        | awk -F"\t" 'BEGIN{{OFS="\t"}}{{print $$6, $$7, ($$3-$$2)}}' \
        | sort -k1,1 -k2,2n -k3,3n \
        | uniq -c \
        | awk 'BEGIN{{OFS="\t"}}{{print $$2, $$3, $$4, $$1}}' \
        | gzip > {output}) > {log} 2>&1
        """


rule fragmentomics_metrics_calculate_SE:
    input:
        se_fragstats="biomarker/fragmentomics_metrics_get_SE_fragstats/{sample}_{type}.txt.gz",
    output:
        se_metrics=temp("biomarker/fragmentomics_metrics_calculate_SE/{sample}_{type}.SE.tsv"),
    log:
        "biomarker/fragmentomics_metrics_calculate_SE/{sample}_{type}.log",
    benchmark:
        repeat(
            "biomarker/fragmentomics_metrics_calculate_SE/{sample}_{type}.benchmark.tsv",
            config.get("fragmentomics_metrics_calculate_SE", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fragmentomics_metrics_calculate_SE", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fragmentomics_metrics_calculate_SE", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fragmentomics_metrics_calculate_SE", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("fragmentomics_metrics_calculate_SE", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("fragmentomics_metrics_calculate_SE", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("fragmentomics_metrics_calculate_SE", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fragmentomics_metrics", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates SE metrics for {wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/fragmentomics_metrics/calculate_SE.R"


rule fragmentomics_metrics_get_depth_fragstats:
    input:
        bed_gz_filtered="biomarker/fragmentomics_metrics_filter_comm_panel_genes/{sample}_{type}.filtered.bed.gz",
    output:
        depth_fragstats=temp("biomarker/fragmentomics_metrics_get_depth_fragstats/{sample}_{type}.txt"),
    log:
        "biomarker/fragmentomics_metrics_get_depth_fragstats/{sample}_{type}.log",
    benchmark:
        repeat(
            "biomarker/fragmentomics_metrics_get_depth_fragstats/{sample}_{type}.benchmark.tsv",
            config.get("fragmentomics_metrics_get_depth_fragstats", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fragmentomics_metrics_get_depth_fragstats", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fragmentomics_metrics_get_depth_fragstats", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fragmentomics_metrics_get_depth_fragstats", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("fragmentomics_metrics_get_depth_fragstats", {}).get(
            "partition", config["default_resources"]["partition"]
        ),
        threads=config.get("fragmentomics_metrics_get_depth_fragstats", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("fragmentomics_metrics_get_depth_fragstats", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fragmentomics_metrics", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates depth fragstats for {wildcards.sample}_{wildcards.type}"
    priority: 20
    shell:
        """
        (zcat {input.bed_gz_filtered} \
        | cut -f 6,7 \
        | sort -k1,1 -k2,2n \
        | uniq -c \
        | awk -F"\t" 'BEGIN{{OFS="\t"}}{{print $$2, $$3, $$1 }}' > {output.depth_fragstats}) > {log} 2>&1
        """


rule fragmentomics_metrics_calculate_normalized_depth:
    input:
        depth_fragstats="biomarker/fragmentomics_metrics_get_depth_fragstats/{sample}_{type}.txt",
    output:
        depth=temp("biomarker/fragmentomics_metrics_calculate_normalized_depth/{sample}_{type}.depth.tsv"),
    params:
        exon_sizes=config.get("fragmentomics_metrics", {}).get("exon_sizes", "resources/UCSC_hg38_exon_sizes.tsv"),
    log:
        "biomarker/fragmentomics_metrics_calculate_normalized_depth/{sample}_{type}.log",
    benchmark:
        repeat(
            "biomarker/fragmentomics_metrics_calculate_normalized_depth/{sample}_{type}.benchmark.tsv",
            config.get("fragmentomics_metrics_calculate_normalized_depth", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fragmentomics_metrics_calculate_normalized_depth", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fragmentomics_metrics_calculate_normalized_depth", {}).get(
            "mem_mb", config["default_resources"]["mem_mb"]
        ),
        mem_per_cpu=config.get("fragmentomics_metrics_calculate_normalized_depth", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("fragmentomics_metrics_calculate_normalized_depth", {}).get(
            "partition", config["default_resources"]["partition"]
        ),
        threads=config.get("fragmentomics_metrics_calculate_normalized_depth", {}).get(
            "threads", config["default_resources"]["threads"]
        ),
        time=config.get("fragmentomics_metrics_calculate_normalized_depth", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fragmentomics_metrics", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates normalized depth for {wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/fragmentomics_metrics/calculate_depth.R"


rule fragmentomics_metrics_calculate_frag_bins:
    input:
        se_fragstats="biomarker/fragmentomics_metrics_get_SE_fragstats/{sample}_{type}.txt.gz",
    output:
        fragbins=temp("biomarker/fragmentomics_metrics_calculate_frag_bins/{sample}_{type}.fragbins.tsv"),
    log:
        "biomarker/fragmentomics_metrics_calculate_frag_bins/{sample}_{type}.log",
    benchmark:
        repeat(
            "biomarker/fragmentomics_metrics_calculate_frag_bins/{sample}_{type}.benchmark.tsv",
            config.get("fragmentomics_metrics_calculate_frag_bins", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fragmentomics_metrics_calculate_frag_bins", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fragmentomics_metrics_calculate_frag_bins", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fragmentomics_metrics_calculate_frag_bins", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("fragmentomics_metrics_calculate_frag_bins", {}).get(
            "partition", config["default_resources"]["partition"]
        ),
        threads=config.get("fragmentomics_metrics_calculate_frag_bins", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("fragmentomics_metrics_calculate_frag_bins", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fragmentomics_metrics", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates fragment bins for {wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/fragmentomics_metrics/calculate_frag_bins.R"


rule fragmentomics_metrics_calculate_small_frags:
    input:
        se_fragstats="biomarker/fragmentomics_metrics_get_SE_fragstats/{sample}_{type}.txt.gz",
    output:
        smallfrag=temp("biomarker/fragmentomics_metrics_calculate_small_frags/{sample}_{type}.smallfrag.tsv"),
    log:
        "biomarker/fragmentomics_metrics_calculate_small_frags/{sample}_{type}.log",
    benchmark:
        repeat(
            "biomarker/fragmentomics_metrics_calculate_small_frags/{sample}_{type}.benchmark.tsv",
            config.get("fragmentomics_metrics_calculate_small_frags", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fragmentomics_metrics_calculate_small_frags", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fragmentomics_metrics_calculate_small_frags", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fragmentomics_metrics_calculate_small_frags", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("fragmentomics_metrics_calculate_small_frags", {}).get(
            "partition", config["default_resources"]["partition"]
        ),
        threads=config.get("fragmentomics_metrics_calculate_small_frags", {}).get(
            "threads", config["default_resources"]["threads"]
        ),
        time=config.get("fragmentomics_metrics_calculate_small_frags", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fragmentomics_metrics", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates small fragment metrics for {wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/fragmentomics_metrics/calculate_small_frags.R"


rule fragmentomics_metrics_calculate_full_gene_depth:
    input:
        depth_fragstats="biomarker/fragmentomics_metrics_get_depth_fragstats/{sample}_{type}.txt",
    output:
        fullgenedepth=temp("biomarker/fragmentomics_metrics_calculate_full_gene_depth/{sample}_{type}.fullgenedepth.tsv"),
    params:
        exon_sizes=config.get("fragmentomics_metrics", {}).get("exon_sizes", "resources/UCSC_hg38_exon_sizes.tsv"),
    log:
        "biomarker/fragmentomics_metrics_calculate_full_gene_depth/{sample}_{type}.log",
    benchmark:
        repeat(
            "biomarker/fragmentomics_metrics_calculate_full_gene_depth/{sample}_{type}.benchmark.tsv",
            config.get("fragmentomics_metrics_calculate_full_gene_depth", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fragmentomics_metrics_calculate_full_gene_depth", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fragmentomics_metrics_calculate_full_gene_depth", {}).get(
            "mem_mb", config["default_resources"]["mem_mb"]
        ),
        mem_per_cpu=config.get("fragmentomics_metrics_calculate_full_gene_depth", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("fragmentomics_metrics_calculate_full_gene_depth", {}).get(
            "partition", config["default_resources"]["partition"]
        ),
        threads=config.get("fragmentomics_metrics_calculate_full_gene_depth", {}).get(
            "threads", config["default_resources"]["threads"]
        ),
        time=config.get("fragmentomics_metrics_calculate_full_gene_depth", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fragmentomics_metrics", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates full gene depth for {wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/fragmentomics_metrics/calculate_full_gene_depth.R"


rule fragmentomics_metrics_get_left_4mer:
    input:
        sample_reads="biomarker/fragmentomics_metrics_filter_comm_panel_genes/{sample}_{type}.filtered.bed.gz",
    output:
        left_4mers=temp("biomarker/fragmentomics_metrics_get_left_4mer/{sample}_{type}_left_4mers.txt.gz"),
    params:
        reference=config.get("reference", {}).get("fasta", ""),
    log:
        "biomarker/fragmentomics_metrics_get_left_4mer/{sample}_{type}.log",
    benchmark:
        repeat(
            "biomarker/fragmentomics_metrics_get_left_4mer/{sample}_{type}.benchmark.tsv",
            config.get("fragmentomics_metrics_get_left_4mer", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fragmentomics_metrics_get_left_4mer", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fragmentomics_metrics_get_left_4mer", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fragmentomics_metrics_get_left_4mer", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("fragmentomics_metrics_get_left_4mer", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("fragmentomics_metrics_get_left_4mer", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("fragmentomics_metrics_get_left_4mer", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fragmentomics_metrics", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates left 4-mer motifs for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        (zcat {input.sample_reads} \
        | awk -F "\t" 'BEGIN{{OFS=FS}}; {{print $$1, $$2, $$2+4, $$4, $$5, $$6, $$7, $$8}}' \
        | bedtools getfasta -fi {params.reference} -bed stdin -bedOut \
        | awk -F "\t" 'BEGIN{{OFS=FS}}; {{print $$6, $$7, toupper($$9)}}' \
        | sort -k1,1 -k2,2n -k3,3 \
        | uniq -c \
        | gzip > {output}) > {log} 2>&1
        """


rule fragmentomics_metrics_get_right_4mer:
    input:
        sample_reads="biomarker/fragmentomics_metrics_filter_comm_panel_genes/{sample}_{type}.filtered.bed.gz",
    output:
        right_4mers=temp("biomarker/fragmentomics_metrics_get_right_4mer/{sample}_{type}_right_4mers.txt.gz"),
    params:
        reference=config.get("reference", {}).get("fasta", ""),
    log:
        "biomarker/fragmentomics_metrics_get_right_4mer/{sample}_{type}.log",
    benchmark:
        repeat(
            "biomarker/fragmentomics_metrics_get_right_4mer/{sample}_{type}.benchmark.tsv",
            config.get("fragmentomics_metrics_get_right_4mer", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fragmentomics_metrics_get_right_4mer", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fragmentomics_metrics_get_right_4mer", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fragmentomics_metrics_get_right_4mer", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("fragmentomics_metrics_get_right_4mer", {}).get(
            "partition", config["default_resources"]["partition"]
        ),
        threads=config.get("fragmentomics_metrics_get_right_4mer", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("fragmentomics_metrics_get_right_4mer", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fragmentomics_metrics", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates right 4-mer motifs for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        (zcat {input.sample_reads} \
        | awk -F "\t" 'BEGIN{{OFS=FS}}; {{print $$1, $$3-4, $$3, $$4, $$5, "-", $$6, $$7, $$8}}' \
        | bedtools getfasta -fi {params.reference} -bed stdin -bedOut -s \
        | awk -F "\t" 'BEGIN{{OFS=FS}}; {{print $$7, $$8, toupper($$10)}}' \
        | sort -k1,1 -k2,2n -k3,3 \
        | uniq -c \
        | gzip > {output}) > {log} 2>&1
        """


rule fragmentomics_metrics_calculate_MDS:
    input:
        left="biomarker/fragmentomics_metrics_get_left_4mer/{sample}_{type}_left_4mers.txt.gz",
        right="biomarker/fragmentomics_metrics_get_right_4mer/{sample}_{type}_right_4mers.txt.gz",
    output:
        mds=temp("biomarker/fragmentomics_metrics_calculate_MDS/{sample}_{type}_mds.txt"),
    log:
        "biomarker/fragmentomics_metrics_calculate_MDS/{sample}_{type}.log",
    benchmark:
        repeat(
            "biomarker/fragmentomics_metrics_calculate_MDS/{sample}_{type}.benchmark.tsv",
            config.get("fragmentomics_metrics_calculate_MDS", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fragmentomics_metrics_calculate_MDS", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fragmentomics_metrics_calculate_MDS", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fragmentomics_metrics_calculate_MDS", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("fragmentomics_metrics_calculate_MDS", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("fragmentomics_metrics_calculate_MDS", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("fragmentomics_metrics_calculate_MDS", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fragmentomics_metrics", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates MDS entropy for {wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/fragmentomics_metrics/calculate_mds.R"


rule fragmentomics_metrics_overlap_TFBS:
    input:
        sample_reads="biomarker/fragmentomics_metrics_filter_comm_panel_genes/{sample}_{type}.filtered.bed.gz",
    output:
        tfbs_frag_count=temp("biomarker/fragmentomics_metrics_overlap_TFBS/{sample}_{type}_TFBS_frag_count.txt.gz"),
    params:
        tfbs=config.get("fragmentomics_metrics", {}).get("tfbs_midpoints", ""),
    log:
        "biomarker/fragmentomics_metrics_overlap_TFBS/{sample}_{type}.log",
    benchmark:
        repeat(
            "biomarker/fragmentomics_metrics_overlap_TFBS/{sample}_{type}.benchmark.tsv",
            config.get("fragmentomics_metrics_overlap_TFBS", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fragmentomics_metrics_overlap_TFBS", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fragmentomics_metrics_overlap_TFBS", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fragmentomics_metrics_overlap_TFBS", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("fragmentomics_metrics_overlap_TFBS", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("fragmentomics_metrics_overlap_TFBS", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("fragmentomics_metrics_overlap_TFBS", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fragmentomics_metrics", {}).get("container", config["default_container"])
    message:
        "{rule}: Overlaps TFBS midpoints for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        (zcat {input.sample_reads} \
        | bedtools intersect -a stdin -b {params.tfbs} -wa -wb \
        | awk -F"\t" 'BEGIN {{OFS="\t"}} {{print $$3-$$2, $$12}}' \
        | sort -k2,2 -k1,1n \
        | uniq -c \
        | awk -F" " 'BEGIN{{OFS="\t"}}{{print $$2, $$3, $$1}}' \
        | gzip  > {output}) > {log} 2>&1
        """


rule fragmentomics_metrics_calculate_TFBS_entropy:
    input:
        tfbs_frag_count="biomarker/fragmentomics_metrics_overlap_TFBS/{sample}_{type}_TFBS_frag_count.txt.gz",
    output:
        tfbs_entropy=temp("biomarker/fragmentomics_metrics_calculate_TFBS_entropy/{sample}_{type}_TFBS_entropy.txt"),
    log:
        "biomarker/fragmentomics_metrics_calculate_TFBS_entropy/{sample}_{type}.log",
    benchmark:
        repeat(
            "biomarker/fragmentomics_metrics_calculate_TFBS_entropy/{sample}_{type}.benchmark.tsv",
            config.get("fragmentomics_metrics_calculate_TFBS_entropy", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fragmentomics_metrics_calculate_TFBS_entropy", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fragmentomics_metrics_calculate_TFBS_entropy", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fragmentomics_metrics_calculate_TFBS_entropy", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("fragmentomics_metrics_calculate_TFBS_entropy", {}).get(
            "partition", config["default_resources"]["partition"]
        ),
        threads=config.get("fragmentomics_metrics_calculate_TFBS_entropy", {}).get(
            "threads", config["default_resources"]["threads"]
        ),
        time=config.get("fragmentomics_metrics_calculate_TFBS_entropy", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fragmentomics_metrics", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates TFBS entropy for {wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/fragmentomics_metrics/calculate_TFBS_entropy.R"


rule fragmentomics_metrics_overlap_ATAC:
    input:
        sample_reads="biomarker/fragmentomics_metrics_filter_comm_panel_genes/{sample}_{type}.filtered.bed.gz",
    output:
        atac_frag_count=temp("biomarker/fragmentomics_metrics_overlap_ATAC/{sample}_{type}_ATAC_frag_count.txt.gz"),
    params:
        atac=config.get("fragmentomics_metrics", {}).get("atac_peaks", ""),
    log:
        "biomarker/fragmentomics_metrics_overlap_ATAC/{sample}_{type}.log",
    benchmark:
        repeat(
            "biomarker/fragmentomics_metrics_overlap_ATAC/{sample}_{type}.benchmark.tsv",
            config.get("fragmentomics_metrics_overlap_ATAC", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fragmentomics_metrics_overlap_ATAC", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fragmentomics_metrics_overlap_ATAC", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fragmentomics_metrics_overlap_ATAC", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("fragmentomics_metrics_overlap_ATAC", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("fragmentomics_metrics_overlap_ATAC", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("fragmentomics_metrics_overlap_ATAC", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fragmentomics_metrics", {}).get("container", config["default_container"])
    message:
        "{rule}: Overlaps ATAC peaks for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        (zcat {input.sample_reads} \
        | bedtools intersect -a stdin -b {params.atac} -wa -wb \
        | awk -F"\t" 'BEGIN{{OFS="\t"}}; {{split($$12, arr, "_"); print arr[1], $$3-$$2}}' \
        | sort -k1,1 -k2,2n \
        | uniq -c \
        | awk -F" " 'BEGIN{{OFS="\t"}}{{print $$2, $$3, $$1}}' \
        | gzip > {output}) > {log} 2>&1
        """


rule fragmentomics_metrics_calculate_ATAC_entropy:
    input:
        atac_frag_count="biomarker/fragmentomics_metrics_overlap_ATAC/{sample}_{type}_ATAC_frag_count.txt.gz",
    output:
        atac_entropy=temp("biomarker/fragmentomics_metrics_calculate_ATAC_entropy/{sample}_{type}_ATAC_entropy.txt"),
    log:
        "biomarker/fragmentomics_metrics_calculate_ATAC_entropy/{sample}_{type}.log",
    benchmark:
        repeat(
            "biomarker/fragmentomics_metrics_calculate_ATAC_entropy/{sample}_{type}.benchmark.tsv",
            config.get("fragmentomics_metrics_calculate_ATAC_entropy", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fragmentomics_metrics_calculate_ATAC_entropy", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fragmentomics_metrics_calculate_ATAC_entropy", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fragmentomics_metrics_calculate_ATAC_entropy", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("fragmentomics_metrics_calculate_ATAC_entropy", {}).get(
            "partition", config["default_resources"]["partition"]
        ),
        threads=config.get("fragmentomics_metrics_calculate_ATAC_entropy", {}).get(
            "threads", config["default_resources"]["threads"]
        ),
        time=config.get("fragmentomics_metrics_calculate_ATAC_entropy", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fragmentomics_metrics", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates ATAC entropy for {wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/fragmentomics_metrics/calculate_ATAC_entropy.R"


rule fragmentomics_metrics_build_feature_tables:
    input:
        se=[
            f"biomarker/fragmentomics_metrics_calculate_SE/{sample}_{type}.SE.tsv"
            for sample in get_samples(samples)
            for type in get_unit_types(units, sample)
        ],
        depth=[
            f"biomarker/fragmentomics_metrics_calculate_normalized_depth/{sample}_{type}.depth.tsv"
            for sample in get_samples(samples)
            for type in get_unit_types(units, sample)
        ],
        frag_bins=[
            f"biomarker/fragmentomics_metrics_calculate_frag_bins/{sample}_{type}.fragbins.tsv"
            for sample in get_samples(samples)
            for type in get_unit_types(units, sample)
        ],
        small_frags=[
            f"biomarker/fragmentomics_metrics_calculate_small_frags/{sample}_{type}.smallfrag.tsv"
            for sample in get_samples(samples)
            for type in get_unit_types(units, sample)
        ],
        full_gene_depth=[
            f"biomarker/fragmentomics_metrics_calculate_full_gene_depth/{sample}_{type}.fullgenedepth.tsv"
            for sample in get_samples(samples)
            for type in get_unit_types(units, sample)
        ],
        mds=[
            f"biomarker/fragmentomics_metrics_calculate_MDS/{sample}_{type}_mds.txt"
            for sample in get_samples(samples)
            for type in get_unit_types(units, sample)
        ],
        TFBS_entropy=[
            f"biomarker/fragmentomics_metrics_calculate_TFBS_entropy/{sample}_{type}_TFBS_entropy.txt"
            for sample in get_samples(samples)
            for type in get_unit_types(units, sample)
        ],
        ATAC_entropy=[
            f"biomarker/fragmentomics_metrics_calculate_ATAC_entropy/{sample}_{type}_ATAC_entropy.txt"
            for sample in get_samples(samples)
            for type in get_unit_types(units, sample)
        ],
    output:
        se=temp("biomarker/fragmentomics_metrics_build_feature_tables/se.rds"),
        depth=temp("biomarker/fragmentomics_metrics_build_feature_tables/depth.rds"),
        frag_bins=temp("biomarker/fragmentomics_metrics_build_feature_tables/frag_bins.rds"),
        small_frags=temp("biomarker/fragmentomics_metrics_build_feature_tables/small_frags.rds"),
        full_gene_depth=temp("biomarker/fragmentomics_metrics_build_feature_tables/full_gene_depth.rds"),
        mds=temp("biomarker/fragmentomics_metrics_build_feature_tables/mds.rds"),
        TFBS_entropy=temp("biomarker/fragmentomics_metrics_build_feature_tables/TFBS_entropy.rds"),
        ATAC_entropy=temp("biomarker/fragmentomics_metrics_build_feature_tables/ATAC_entropy.rds"),
    log:
        "biomarker/fragmentomics_metrics_build_feature_tables/build_feature_tables.log",
    benchmark:
        repeat(
            "biomarker/fragmentomics_metrics_build_feature_tables/build_feature_tables.benchmark.tsv",
            config.get("fragmentomics_metrics_build_feature_tables", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fragmentomics_metrics_build_feature_tables", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fragmentomics_metrics_build_feature_tables", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fragmentomics_metrics_build_feature_tables", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("fragmentomics_metrics_build_feature_tables", {}).get(
            "partition", config["default_resources"]["partition"]
        ),
        threads=config.get("fragmentomics_metrics_build_feature_tables", {}).get(
            "threads", config["default_resources"]["threads"]
        ),
        time=config.get("fragmentomics_metrics_build_feature_tables", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fragmentomics_metrics", {}).get("container", config["default_container"])
    message:
        "{rule}: Builds feature tables for all samples"
    script:
        "../scripts/fragmentomics_metrics/build_feature_tables.R"


rule fragmentomics_metrics_extract_first_exon:
    input:
        se="biomarker/fragmentomics_metrics_build_feature_tables/se.rds",
        depth="biomarker/fragmentomics_metrics_build_feature_tables/depth.rds",
        frag_bins="biomarker/fragmentomics_metrics_build_feature_tables/frag_bins.rds",
        small_frags="biomarker/fragmentomics_metrics_build_feature_tables/small_frags.rds",
        full_gene_depth="biomarker/fragmentomics_metrics_build_feature_tables/full_gene_depth.rds",
        mds="biomarker/fragmentomics_metrics_build_feature_tables/mds.rds",
        TFBS_entropy="biomarker/fragmentomics_metrics_build_feature_tables/TFBS_entropy.rds",
        ATAC_entropy="biomarker/fragmentomics_metrics_build_feature_tables/ATAC_entropy.rds",
    output:
        se_E1="biomarker/fragmentomics_metrics_extract_first_exon/se_E1.rds",
        depth_E1="biomarker/fragmentomics_metrics_extract_first_exon/depth_E1.rds",
        small_frags_E1="biomarker/fragmentomics_metrics_extract_first_exon/small_frags_E1.rds",
        mds_E1="biomarker/fragmentomics_metrics_extract_first_exon/mds_E1.rds",
    params:
        strand_mapping=config.get("fragmentomics_metrics", {}).get("exon_strand_mapping", ""),
    log:
        "biomarker/fragmentomics_metrics_extract_first_exon/extract_first_exon.log",
    benchmark:
        repeat(
            "biomarker/fragmentomics_metrics_extract_first_exon/extract_first_exon.benchmark.tsv",
            config.get("fragmentomics_metrics_extract_first_exon", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fragmentomics_metrics_extract_first_exon", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fragmentomics_metrics_extract_first_exon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fragmentomics_metrics_extract_first_exon", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("fragmentomics_metrics_extract_first_exon", {}).get(
            "partition", config["default_resources"]["partition"]
        ),
        threads=config.get("fragmentomics_metrics_extract_first_exon", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("fragmentomics_metrics_extract_first_exon", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fragmentomics_metrics", {}).get("container", config["default_container"])
    message:
        "{rule}: Extracts first exons for feature tables"
    script:
        "../scripts/fragmentomics_metrics/filter_first_exons.R"


rule fragmentomics_metrics_build_combined_ft:
    input:
        se="biomarker/fragmentomics_metrics_build_feature_tables/se.rds",
        depth="biomarker/fragmentomics_metrics_build_feature_tables/depth.rds",
        frag_bins="biomarker/fragmentomics_metrics_build_feature_tables/frag_bins.rds",
        small_frags="biomarker/fragmentomics_metrics_build_feature_tables/small_frags.rds",
        full_gene_depth="biomarker/fragmentomics_metrics_build_feature_tables/full_gene_depth.rds",
        mds="biomarker/fragmentomics_metrics_build_feature_tables/mds.rds",
        TFBS_entropy="biomarker/fragmentomics_metrics_build_feature_tables/TFBS_entropy.rds",
        ATAC_entropy="biomarker/fragmentomics_metrics_build_feature_tables/ATAC_entropy.rds",
    output:
        combined_rds="biomarker/fragmentomics_metrics_build_combined_ft/all_combined.rds",
    log:
        "biomarker/fragmentomics_metrics_build_combined_ft/build_combined_ft.log",
    benchmark:
        repeat(
            "biomarker/fragmentomics_metrics_build_combined_ft/build_combined_ft.benchmark.tsv",
            config.get("fragmentomics_metrics_build_combined_ft", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fragmentomics_metrics_build_combined_ft", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fragmentomics_metrics_build_combined_ft", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fragmentomics_metrics_build_combined_ft", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("fragmentomics_metrics_build_combined_ft", {}).get(
            "partition", config["default_resources"]["partition"]
        ),
        threads=config.get("fragmentomics_metrics_build_combined_ft", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("fragmentomics_metrics_build_combined_ft", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fragmentomics_metrics", {}).get("container", config["default_container"])
    message:
        "{rule}: Builds combined feature table for all samples"
    script:
        "../scripts/fragmentomics_metrics/generate_combined_ft.R"
