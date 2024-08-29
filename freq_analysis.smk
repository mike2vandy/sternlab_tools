import re
import pandas as pd

singularity: '/panfs/jay/groups/0/fried255/shared/gatk4_workflow/rescuer/AFREQ/CalcFix/afreq.sif'

# MAKE THIS A PART OF THE AF CONTAINER - NO NEED TO REGEN EVERYTIME
chrmsDF = pd.read_table("/home/fried255/mvandewe/software/afs_in_vcf/dog.regions.txt")
regions = list(chrmsDF['chrom'])

rule all:
    input:
       #expand(
       #   #"results/{region}/{region}.split.vcf.gz",
       #    "results/{region}/{region}.csv",
       #    region=regions
       #   #region='chr27.31108325-46662488'
       #),
        "af_analysis.no_intergenic.csv"  

# ADD MORE STRUCTURE TO DIRS LIKE split/{region}/{region}.vcf.gz
rule select_variants_chrom:
    input:
        vcf = "joint_call.UU_Cfam_GSD_1.0_ROSY.20240718.vep.vcf.gz" 
    output:
        vcf = "output/{region}/{region}.split.vcf.gz"
    params:
        reg = lambda wildcards, output: re.sub(r'\.(?=[^.]*$)', ':', wildcards.region)
    #conda: "env/tools.yaml"
    threads: 12
    resources:
        time = 60,
        mem_mb = 40000
    shell:
        '''
            bcftools view \
                -o {output.vcf} \
                -r {params.reg} \
                --threads {threads} \
                -Oz \
                {input.vcf}
        '''

rule process_variants:
    input:
        "lists/tvds.txt",
        "lists/ctls.txt",
        "lists/oavrts.txt",
        vcf = "output/{region}/{region}.split.vcf.gz" 
    output:
        table = "output/{region}/{region}.csv"
    threads: 1
    #conda: "env/tools.yaml"
    resources:
        time = 1440,
        mem_mb = 20000
    shell:
        '''
            ./calcAFsFromGroups_noSift.py \
                {input.vcf} \
                {input[0]} \
                {input[1]} \
                {input[2]}  > {output.table}
        '''

rule concat_files:
    input:
        expand("output/{region}/{region}.csv", region=regions)
    output:
        "af_analysis.csv"  
    threads: 1
    resources:
        time = 600,
        mem_mb = 120000
    run:
        dfs = [pd.read_csv(f, low_memory=False) for f in input]
        concat_df = pd.concat(dfs)
        concat_df.to_csv(output[0], index = False)

localrules: remove_intergenic
rule remove_intergenic:
    input:
        "af_analysis.csv" 
    output:
        "af_analysis.no_intergenic.csv" 
    shell:
        '''
            awk 'NR == 1 || !/intergenic/' {input} > {output}
        '''

