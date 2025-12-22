
import pandas as pd

samples = pd.read_csv("samples.txt", header = None)[0]

rule all:
  input:
    "count_matrix/mmvd_dog.rc.csv"
    #expand("mapped/{sample}/{sample}ReadsPerGene.out.tab", sample = samples)

rule index_genome:
  input:
    fas = "ref/UU_Cfam_GSD_1.0_ROSY.fa"
  output:
    "/share/stern/mwvandew/mmvd/rna/ref/SAindex"
  params:
    ref_dir = "ref/",
    gtf = "ref/UU_Cfam_GSD_1.0_ROSY.refSeq.ensformat.gtf",
    overhang = 149
  conda: "env/star.yaml"
  threads: 12
  resources:
    time = 120,
    mem_mb = 40000
  shell:
    '''
      STAR \
        --runThreadN {threads} \
	--runMode genomeGenerate \
        --genomeFastaFiles {input.fas} \
        --genomeDir {params.ref_dir} \
        --sjdbOverhang {params.overhang} \
        --sjdbGTFfile {params.gtf}
    '''

rule trim:
  input:
    fq1 = "raw_fq/{sample}_R1_001.fastq.gz",
    fq2 = "raw_fq/{sample}_R2_001.fastq.gz"
  output:
    pe1 = "trimmed/{sample}_1.fq.gz",
    pe2 = "trimmed/{sample}_2.fq.gz",
    se1 = "trash/{sample}_1.fq.gz",
    se2 = "trash/{sample}_2.fq.gz"
  params:
    adapter = "/usr/local/usrapps/stern/mwvandew/conda/envs/star/share/trimmomatic-0.39-2/adapters/allAdapters.fa",
    length = 100
  conda: "env/star.yaml"
  threads: 12
  resources:  
    time = 200,
    mem_mb = 40000
  shell:
    '''
      trimmomatic PE \
        -threads {threads} \
        {input.fq1} {input.fq2} \
        {output.pe1} {output.se1} \
        {output.pe2} {output.se2} \
        ILLUMINACLIP:{params.adapter}:2:30:10 \
        SLIDINGWINDOW:5:20 \
        MINLEN:{params.length}
    '''

rule map:
  input:
    "/share/stern/mwvandew/mmvd/rna/ref/SAindex",
    fq1 = "trimmed/{sample}_1.fq.gz",
    fq2 = "trimmed/{sample}_2.fq.gz"
  output:
    "mapped/{sample}/{sample}ReadsPerGene.out.tab"
  params:
    ref_dir = "/share/stern/mwvandew/mmvd/rna/ref",
    gtf = "/share/stern/mwvandew/mmvd/rna/ref/UU_Cfam_GSD_1.0_ROSY.refSeq.ensformat.gtf",
    out_prefix = "mapped/{sample}/{sample}"
  conda: "env/star.yaml"
  threads: 6
  resources:
    time = 200,
    mem_mb = 40000
  shell:
    '''
      STAR \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --readFilesCommand zcat \
        --runThreadN {threads} \
        --genomeDir {params.ref_dir} \
        --sjdbGTFfile {params.gtf} \
        --readFilesIn {input.fq1} {input.fq2} \
        --outFileNamePrefix {params.out_prefix}
    '''

localrules: join_counts
rule join_counts:
  input:
    expand("mapped/{sample}/{sample}ReadsPerGene.out.tab", sample = samples)
  output:
    csv = "count_matrix/mmvd_dog.rc.csv"
  run:
    from functools import reduce
    import os
    
    dfs = []

    for f in input:
      sample = os.path.basename(os.path.dirname(f))
      df = pd.read_csv(f, sep = '\t', header = None, skiprows = 4)
      df = df[[0,1]]
      df.columns = ["GeneId", sample]
      dfs.append(df)

    merged = reduce(lambda left, right: pd.merge(left, right, on = "GeneId"), dfs)
    merged.to_csv(output.csv, index = False)
