rule all:
    input:
        "plots/read_depth_scatter.png"
rule download_data:
    output:
        "data/tu.r1.fq.gz",
        "data/tu.r2.fq.gz",
        "data/wt.r1.fq.gz",
        "data/wt.r2.fq.gz"
    log:
        "logs/download_data.log"
    params:
        memory="4"
    threads:
        1
    shell:
        """
        # download data to directorory data
        wget -O data/tu.r1.fq.gz https://gear.embl.de/data/.exercise/tu.r1.fq.gz &> {log}
        wget -O data/tu.r2.fq.gz https://gear.embl.de/data/.exercise/tu.r2.fq.gz &>> {log}
        wget -O data/wt.r1.fq.gz https://gear.embl.de/data/.exercise/wt.r1.fq.gz &>> {log}
        wget -O data/wt.r2.fq.gz https://gear.embl.de/data/.exercise/wt.r2.fq.gz &>> {log}
        """

rule download_reference:
    output:
        "data/reference.fa.gz"
    log:
        "logs/download_reference.log"
    params:
        memory="10"
    threads:
        1
    shell:
        """
        # download reference with filename reference.fa.gz to directory data
        wget -O data/reference.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz &> {log}
        """
rule gunzip_data:
    input:
        "data/{name}.{ext}.gz"
    output:
        "data/{name}.{ext}"
    log:
        "logs/gunzip_data_{name}_{ext}.log"
    params:
        memory="4"
    threads:
        1
    shell:
        """
        gunzip -c {input} > {output} 2> {log}
        """

rule index_reference:
  """
  Indexes the reference sequence.
  """
  input:
    reference = "data/reference.fa"
  output:
    "data/reference.fa.bwt"
  log:
    "logs/index_reference.log"
  params:
    memory="4"
  threads:
    1
  shell:
    """
    bwa index {input.reference} &> {log}
    """

rule align_reads:
    """
    Aligns the reads to the reference sequence.
    """
    input:
        reference = "data/reference.fa",
        reference_index = "data/reference.fa.bwt",
        r1_reads = "data/{sample}.r1.fq.gz",
        r2_reads = "data/{sample}.r2.fq.gz"
    output:
        "data/{sample}_align.bam"
    log:
        "logs/align_reads_{sample}.log"
    params:
        memory="4"
    threads:
        6
    shell:
        """
        bwa mem -t {threads} {input.reference} {input.r1_reads} {input.r2_reads} | samtools sort -@ {threads} -o {output} &> {log}
        """
        
rule index_align:
    """
    Indexes the alignment.
    """
    input:
        "data/{sample}_align.bam"
    output:
        "data/{sample}_align.bam.bai"
    log:
        "logs/index_align_{sample}.log"
    params:
        memory="4"
    threads:
        1
    shell:
        """
        samtools index {input} &> {log}
        """

rule subset_align:
    """
    Subset the BAM to the region of interest - chromosome X from 20Mbp to 40Mbp (GRCh37/hg19 coordinates)
    """
    input:
        alignment = "data/{sample}_align.bam",
        index="data/{sample}_align.bam.bai"
    output:
        "data/{sample}_align_subset.bam"
    log:
        "logs/subset_align_{sample}.log"
    params:
        memory="4"
    threads:
        1
    shell:
        """
        samtools view -b {input.alignment} chrX:20000000-40000000 > {output} 2> {log}
        """

rule get_depth:
    """
    Get the read depth of the subsetted BAM file.
    """
    input:
        alignment = "data/{sample}_align_subset.bam"
    output:
        "data/{sample}_read_depth.tsv"
    log:
        "logs/read_depth_{sample}.log"
    params:
        memory="4"
    threads:
        1
    shell:
        """
        samtools depth {input.alignment} > {output} 2> {log}
        """

rule plot:
    input:
        tu_depth = "data/tu_read_depth.tsv",
        wt_depth = "data/wt_read_depth.tsv"
    output:
        plot1 = "plots/read_depth_scatter.png",
        plot2 = "plots/log2_read_depth.png"
    log:
        "logs/read_depth_plot.log"
    params:
        memory="20"
    threads:
        1
    shell:
        """
        module add conda-modules-py37
        conda activate /storage/plzen1/home/calounovat/.conda/envs/cg_data_analysis_env
        scripts/plot.R {input.tu_depth} {input.wt_depth} {output.plot1} {output.plot2} 2> {log} 
        """
