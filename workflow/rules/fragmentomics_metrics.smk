# Fragmentomics Metrics rules for Hydra Genetics Biomarker module


rule fragmentomics_metrics_get_bed_from_bam:
    input:
        bam="alignment/bwa_mem_realign_consensus_reads/{sample}_{type}.umi.bam",
        bai="alignment/bwa_mem_realign_consensus_reads/{sample}_{type}.umi.bam.bai",
    output:
        "biomarker/fragmentomics_metrics/data/{sample}_{type}.bed.gz",
    params:
        canonical_cds_bed=config.get("fragmentomics_metrics", {}).get(
            "canonical_cds_bed", "resources/UCSC_hg19_canonical_cds.bed"
        ),
        min_mapq=config.get("fragmentomics_metrics", {}).get("min_mapq", 30),
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
        config.get("fragmentomics_metrics_get_bed_from_bam", {}).get("container", config["default_container"])
    message:
        "{rule}: Extracts BED from BAM for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        (samtools view -@ {threads} -b -q {params.min_mapq} {input.bam} \
        | bedtools bamtobed -i stdin \
        | bedtools intersect -wa -wb -a stdin -b {params.canonical_cds_bed} \
        | awk 'BEGIN{{OFS=\"\\t\"}}{{print $$1,$$2,$$3,$$4,$$5,$$12,$$13,$$6}}' \
        | sort --parallel={threads} -S 30G -k1,1 -k2,2n \
        | gzip > {output}) > {log} 2>&1
        """


rule fragmentomics_metrics_filter_comm_panel_genes:
    input:
        "biomarker/fragmentomics_metrics/data/{sample}_{type}.bed.gz",
    output:
        temp("biomarker/fragmentomics_metrics/data/{sample}_{type}.bed"),
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
        config.get("fragmentomics_metrics_filter_comm_panel_genes", {}).get("container", config["default_container"])
    message:
        "{rule}: Filters commercial panel genes for {wildcards.sample}_{wildcards.type}"
    shell:
        "Rscript workflow/scripts/fragmentomics_metrics/filter_comm_panels.R {input} > {log} 2>&1"


rule fragmentomics_metrics_gzip_comm_panel_genes:
    input:
        "biomarker/fragmentomics_metrics/data/{sample}_{type}.bed",
    output:
        "biomarker/fragmentomics_metrics/data/{sample}_{type}.bed.gz_filtered",
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
        config.get("fragmentomics_metrics_gzip_comm_panel_genes", {}).get("container", config["default_container"])
    message:
        "{rule}: Gzips filtered commercial panel genes for {wildcards.sample}_{wildcards.type}"
    shell:
        "gzip -c {input} > {output} 2> {log}"


rule fragmentomics_metrics_get_SE_fragstats:
    input:
        "biomarker/fragmentomics_metrics/data/{sample}_{type}.bed.gz_filtered",
    output:
        "biomarker/fragmentomics_metrics/fragstats/SE_files/{sample}_{type}.txt.gz",
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
        config.get("fragmentomics_metrics_get_SE_fragstats", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates SE fragstats for {wildcards.sample}_{wildcards.type}"
    priority: 20
    shell:
        """
        (zcat {input} \
        | awk '{{print $6, \"\\t\", $7, \"\\t\", ($3-$2)}}' \
        | sort -k1,1 -k2,2n -k3,3n \
        | uniq -c \
        | awk '{{print $2, \"\\t\", $3, \"\\t\", $4, \"\\t\", $1}}' \
        | gzip > {output}) > {log} 2>&1
        """


rule fragmentomics_metrics_calculate_SE:
    input:
        "biomarker/fragmentomics_metrics/fragstats/SE_files/{sample}_{type}.txt.gz",
    output:
        "biomarker/fragmentomics_metrics/metrics/se/{sample}_{type}.SE.tsv",
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
        config.get("fragmentomics_metrics_calculate_SE", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates SE metrics for {wildcards.sample}_{wildcards.type}"
    shell:
        "Rscript workflow/scripts/fragmentomics_metrics/calculate_SE.R {input} {output} > {log} 2>&1"


rule fragmentomics_metrics_get_depth_fragstats:
    input:
        "biomarker/fragmentomics_metrics/data/{sample}_{type}.bed.gz_filtered",
    output:
        "biomarker/fragmentomics_metrics/fragstats/depth_files/{sample}_{type}.txt",
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
        config.get("fragmentomics_metrics_get_depth_fragstats", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates depth fragstats for {wildcards.sample}_{wildcards.type}"
    priority: 20
    shell:
        """
        (zcat {input} \
        | cut -f 6,7 \
        | sort -k1,1 -k2,2n \
        | uniq -c \
        | awk '{{print $2, \"\\t\", $3, \"\\t\", $1 }}' > {output}) > {log} 2>&1
        """


rule fragmentomics_metrics_calculate_normalized_depth:
    input:
        "biomarker/fragmentomics_metrics/fragstats/depth_files/{sample}_{type}.txt",
    output:
        "biomarker/fragmentomics_metrics/metrics/depth/{sample}_{type}.depth.tsv",
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
        config.get("fragmentomics_metrics_calculate_normalized_depth", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates normalized depth for {wildcards.sample}_{wildcards.type}"
    shell:
        "Rscript workflow/scripts/fragmentomics_metrics/calculate_depth.R {params.exon_sizes} {input} {output} > {log} 2>&1"


rule fragmentomics_metrics_calculate_frag_bins:
    input:
        "biomarker/fragmentomics_metrics/fragstats/SE_files/{sample}_{type}.txt.gz",
    output:
        "biomarker/fragmentomics_metrics/metrics/frag_bins/{sample}_{type}.fragbins.tsv",
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
        config.get("fragmentomics_metrics_calculate_frag_bins", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates fragment bins for {wildcards.sample}_{wildcards.type}"
    shell:
        "Rscript workflow/scripts/fragmentomics_metrics/calculate_frag_bins.R {input} {output} > {log} 2>&1"


rule fragmentomics_metrics_calculate_small_frags:
    input:
        "biomarker/fragmentomics_metrics/fragstats/SE_files/{sample}_{type}.txt.gz",
    output:
        "biomarker/fragmentomics_metrics/metrics/small_frags/{sample}_{type}.smallfrag.tsv",
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
        config.get("fragmentomics_metrics_calculate_small_frags", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates small fragment metrics for {wildcards.sample}_{wildcards.type}"
    shell:
        "Rscript workflow/scripts/fragmentomics_metrics/calculate_small_frags.R {input} {output} > {log} 2>&1"


rule fragmentomics_metrics_calculate_full_gene_depth:
    input:
        "biomarker/fragmentomics_metrics/fragstats/depth_files/{sample}_{type}.txt",
    output:
        "biomarker/fragmentomics_metrics/metrics/full_gene_depth/{sample}_{type}.fullgenedepth.tsv",
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
        config.get("fragmentomics_metrics_calculate_full_gene_depth", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates full gene depth for {wildcards.sample}_{wildcards.type}"
    shell:
        "Rscript workflow/scripts/fragmentomics_metrics/calculate_full_gene_depth.R {params.exon_sizes} {input} {output} > {log} 2>&1"


rule fragmentomics_metrics_get_left_4mer:
    input:
        sample_reads="biomarker/fragmentomics_metrics/data/{sample}_{type}.bed.gz_filtered",
    output:
        temp("biomarker/fragmentomics_metrics/metrics/mds/{sample}_{type}_left_4mers.txt.gz"),
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
        config.get("fragmentomics_metrics_get_left_4mer", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates left 4-mer motifs for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        (zcat {input.sample_reads} \
        | awk -F \"\\t\" '{{OFS=FS}}; {{print $1, $2, $2+4, $4, $5, $6, $7, $8}}' \
        | bedtools getfasta -fi {params.reference} -bed stdin -bedOut \
        | awk -F \"\\t\" '{{OFS=FS}}; {{print $6, $7, toupper($9)}}' \
        | sort -k1,1 -k2,2n -k3,3 \
        | uniq -c \
        | gzip > {output}) > {log} 2>&1
        """


rule fragmentomics_metrics_get_right_4mer:
    input:
        sample_reads="biomarker/fragmentomics_metrics/data/{sample}_{type}.bed.gz_filtered",
    output:
        temp("biomarker/fragmentomics_metrics/metrics/mds/{sample}_{type}_right_4mers.txt.gz"),
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
        config.get("fragmentomics_metrics_get_right_4mer", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates right 4-mer motifs for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        (zcat {input.sample_reads} \
        | awk -F \"\\t\" '{{OFS=FS}}; {{print $1, $3-4, $3, $4, $5, \"-\", $6, $7, $8}}' \
        | bedtools getfasta -fi {params.reference} -bed stdin -bedOut -s \
        | awk -F \"\\t\" '{{OFS=FS}}; {{print $7, $8, toupper($10)}}' \
        | sort -k1,1 -k2,2n -k3,3 \
        | uniq -c \
        | gzip > {output}) > {log} 2>&1
        """


rule fragmentomics_metrics_calculate_MDS:
    input:
        left="biomarker/fragmentomics_metrics/metrics/mds/{sample}_{type}_left_4mers.txt.gz",
        right="biomarker/fragmentomics_metrics/metrics/mds/{sample}_{type}_right_4mers.txt.gz",
    output:
        "biomarker/fragmentomics_metrics/metrics/mds/{sample}_{type}_mds.txt",
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
        config.get("fragmentomics_metrics_calculate_MDS", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates MDS entropy for {wildcards.sample}_{wildcards.type}"
    shell:
        "Rscript workflow/scripts/fragmentomics_metrics/calculate_mds.R {input.left} {input.right} {output} > {log} 2>&1"


rule fragmentomics_metrics_overlap_TFBS:
    input:
        sample_reads="biomarker/fragmentomics_metrics/data/{sample}_{type}.bed.gz_filtered",
    output:
        temp("biomarker/fragmentomics_metrics/metrics/TFBS_entropy/{sample}_{type}_TFBS_frag_count.txt.gz"),
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
        config.get("fragmentomics_metrics_overlap_TFBS", {}).get("container", config["default_container"])
    message:
        "{rule}: Overlaps TFBS midpoints for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        (zcat {input.sample_reads} \
        | bedtools intersect -a stdin -b {params.tfbs} -wa -wb \
        | awk -F\"\\t\" 'BEGIN {{OFS=FS}} {{print $3-$2, $12}}' \
        | sort -k2,2 -k1,1n \
        | uniq -c \
        | gzip  > {output}) > {log} 2>&1
        """


rule fragmentomics_metrics_calculate_TFBS_entropy:
    input:
        "biomarker/fragmentomics_metrics/metrics/TFBS_entropy/{sample}_{type}_TFBS_frag_count.txt.gz",
    output:
        "biomarker/fragmentomics_metrics/metrics/TFBS_entropy/{sample}_{type}_TFBS_entropy.txt",
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
        config.get("fragmentomics_metrics_calculate_TFBS_entropy", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates TFBS entropy for {wildcards.sample}_{wildcards.type}"
    shell:
        "Rscript workflow/scripts/fragmentomics_metrics/calculate_TFBS_entropy.R {input} {output} > {log} 2>&1"


rule fragmentomics_metrics_overlap_ATAC:
    input:
        sample_reads="biomarker/fragmentomics_metrics/data/{sample}_{type}.bed.gz_filtered",
    output:
        temp("biomarker/fragmentomics_metrics/metrics/ATAC_entropy/{sample}_{type}_ATAC_frag_count.txt.gz"),
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
        config.get("fragmentomics_metrics_overlap_ATAC", {}).get("container", config["default_container"])
    message:
        "{rule}: Overlaps ATAC peaks for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        (zcat {input.sample_reads} \
        | bedtools intersect -a stdin -b {params.atac} -wa -wb \
        | awk -F\"\\t\" '{{OFS=FS}}; {{split($12, arr, \"_\"); print arr[1], $3-$2}}' \
        | sort -k1,1 -k2,2n \
        | uniq -c \
        | gzip > {output}) > {log} 2>&1
        """


rule fragmentomics_metrics_calculate_ATAC_entropy:
    input:
        "biomarker/fragmentomics_metrics/metrics/ATAC_entropy/{sample}_{type}_ATAC_frag_count.txt.gz",
    output:
        "biomarker/fragmentomics_metrics/metrics/ATAC_entropy/{sample}_{type}_ATAC_entropy.txt",
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
        config.get("fragmentomics_metrics_calculate_ATAC_entropy", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates ATAC entropy for {wildcards.sample}_{wildcards.type}"
    shell:
        "Rscript workflow/scripts/fragmentomics_metrics/calculate_ATAC_entropy.R {input} {output} > {log} 2>&1"


rule fragmentomics_metrics_build_feature_tables:
    input:
        expand(
            "biomarker/fragmentomics_metrics/metrics/se/{sample}_{type}.SE.tsv",
            sample=get_samples(samples),
            type=get_unit_types(units, samples.index[0]),
        ),
        expand(
            "biomarker/fragmentomics_metrics/metrics/depth/{sample}_{type}.depth.tsv",
            sample=get_samples(samples),
            type=get_unit_types(units, samples.index[0]),
        ),
        expand(
            "biomarker/fragmentomics_metrics/metrics/frag_bins/{sample}_{type}.fragbins.tsv",
            sample=get_samples(samples),
            type=get_unit_types(units, samples.index[0]),
        ),
        expand(
            "biomarker/fragmentomics_metrics/metrics/small_frags/{sample}_{type}.smallfrag.tsv",
            sample=get_samples(samples),
            type=get_unit_types(units, samples.index[0]),
        ),
        expand(
            "biomarker/fragmentomics_metrics/metrics/full_gene_depth/{sample}_{type}.fullgenedepth.tsv",
            sample=get_samples(samples),
            type=get_unit_types(units, samples.index[0]),
        ),
        expand(
            "biomarker/fragmentomics_metrics/metrics/mds/{sample}_{type}_mds.txt",
            sample=get_samples(samples),
            type=get_unit_types(units, samples.index[0]),
        ),
        expand(
            "biomarker/fragmentomics_metrics/metrics/TFBS_entropy/{sample}_{type}_TFBS_entropy.txt",
            sample=get_samples(samples),
            type=get_unit_types(units, samples.index[0]),
        ),
        expand(
            "biomarker/fragmentomics_metrics/metrics/ATAC_entropy/{sample}_{type}_ATAC_entropy.txt",
            sample=get_samples(samples),
            type=get_unit_types(units, samples.index[0]),
        ),
    output:
        "biomarker/fragmentomics_metrics/feature_tables/se.rds",
        "biomarker/fragmentomics_metrics/feature_tables/depth.rds",
        "biomarker/fragmentomics_metrics/feature_tables/frag_bins.rds",
        "biomarker/fragmentomics_metrics/feature_tables/small_frags.rds",
        "biomarker/fragmentomics_metrics/feature_tables/full_gene_depth.rds",
        "biomarker/fragmentomics_metrics/feature_tables/mds.rds",
        "biomarker/fragmentomics_metrics/feature_tables/TFBS_entropy.rds",
        "biomarker/fragmentomics_metrics/feature_tables/ATAC_entropy.rds",
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
        config.get("fragmentomics_metrics_build_feature_tables", {}).get("container", config["default_container"])
    message:
        "{rule}: Builds feature tables for all samples"
    shell:
        "Rscript workflow/scripts/fragmentomics_metrics/build_feature_tables.R > {log} 2>&1"


rule fragmentomics_metrics_extract_first_exon:
    input:
        "biomarker/fragmentomics_metrics/feature_tables/se.rds",
        "biomarker/fragmentomics_metrics/feature_tables/depth.rds",
        "biomarker/fragmentomics_metrics/feature_tables/frag_bins.rds",
        "biomarker/fragmentomics_metrics/feature_tables/small_frags.rds",
        "biomarker/fragmentomics_metrics/feature_tables/full_gene_depth.rds",
        "biomarker/fragmentomics_metrics/feature_tables/mds.rds",
        "biomarker/fragmentomics_metrics/feature_tables/TFBS_entropy.rds",
        "biomarker/fragmentomics_metrics/feature_tables/ATAC_entropy.rds",
    output:
        "biomarker/fragmentomics_metrics/feature_tables/se_E1.rds",
        "biomarker/fragmentomics_metrics/feature_tables/depth_E1.rds",
        "biomarker/fragmentomics_metrics/feature_tables/small_frags_E1.rds",
        "biomarker/fragmentomics_metrics/feature_tables/mds_E1.rds",
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
        config.get("fragmentomics_metrics_extract_first_exon", {}).get("container", config["default_container"])
    message:
        "{rule}: Extracts first exons for feature tables"
    shell:
        "Rscript workflow/scripts/fragmentomics_metrics/filter_first_exons.R {params.strand_mapping} > {log} 2>&1"


rule fragmentomics_metrics_build_combined_ft:
    input:
        "biomarker/fragmentomics_metrics/feature_tables/se.rds",
        "biomarker/fragmentomics_metrics/feature_tables/depth.rds",
        "biomarker/fragmentomics_metrics/feature_tables/frag_bins.rds",
        "biomarker/fragmentomics_metrics/feature_tables/small_frags.rds",
        "biomarker/fragmentomics_metrics/feature_tables/full_gene_depth.rds",
        "biomarker/fragmentomics_metrics/feature_tables/mds.rds",
        "biomarker/fragmentomics_metrics/feature_tables/TFBS_entropy.rds",
        "biomarker/fragmentomics_metrics/feature_tables/ATAC_entropy.rds",
    output:
        "biomarker/fragmentomics_metrics/feature_tables/all_combined.rds",
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
        config.get("fragmentomics_metrics_build_combined_ft", {}).get("container", config["default_container"])
    message:
        "{rule}: Builds combined feature table for all samples"
    shell:
        "Rscript workflow/scripts/fragmentomics_metrics/generate_combined_ft.R > {log} 2>&1"
