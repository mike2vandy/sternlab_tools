import pandas as pd

singularity: "/home/fried255/mvandewe/software/megane/MEGAnE.sif"

rule all:
  input:
    "/scratch.global/friedlab_MV_MEIs/MEGAnE/merged/meh/meh_biallelic.vcf.gz"
    
rule dirList:
  input:
    confi = "meh.dogs.txt"
  output:
    conList = "meh.dirlist"
  params:
    time = 10,
    mem_mb = 20000
  shell:
    '''
      awk '{{print "MEGAnE/indCall/"$1}}' {input.confi} > {output.conList}
    '''
    
rule joinInsert:
  input:
    dirlist = "meh.dirlist" 
  output:
    vcf = "/scratch.global/friedlab_MV_MEIs/MEGAnE/merged/meh/meh_MEI_jointcall.vcf.gz"
  params:
    genome = "/home/fried255/mvandewe/universalDat/UU_Cfam_GSD_ROSY/UU_Cfam_GSD_1.0_ROSY.fa",
    lib = "/home/fried255/mvandewe/software/megane/UU_Cfam/library.fa",
    out_dir = "/scratch.global/friedlab_MV_MEIs/MEGAnE/merged/meh",
    cohort = "meh" 
  threads: 8
  resources:
    time = 600,
    mem_mb = 60000
  shell:
    '''
      joint_calling_hs \
        -merge_mei \
        -f {input.dirlist} \
        -fa {params.genome} \
        -rep {params.lib} \
        -outdir {params.out_dir} \
        -cohort_name {params.cohort} \
        -p {threads}
    '''

rule joinDels:
  input:
    dirlist = "meh.dirlist"
  output:
    vcf = "/scratch.global/friedlab_MV_MEIs/MEGAnE/merged/meh/meh_MEA_jointcall.vcf.gz"
  params:
    genome = "/home/fried255/mvandewe/universalDat/UU_Cfam_GSD_ROSY/UU_Cfam_GSD_1.0_ROSY.fa",
    lib = "/home/fried255/mvandewe/software/megane/UU_Cfam/library.fa",
    out_dir = "/scratch.global/friedlab_MV_MEIs/MEGAnE/merged/meh",
    cohort = "meh"
  threads: 8
  resources:
    time = 600,
    mem_mb = 60000
  shell:
    '''
      joint_calling_hs \
        -merge_absent_me \
        -f {input.dirlist} \
        -fa {params.genome} \
        -rep {params.lib} \
        -outdir {params.out_dir} \
        -cohort_name {params.cohort} \
        -p {threads}
    '''

rule joinAll:
  input:
    mei = "/scratch.global/friedlab_MV_MEIs/MEGAnE/merged/meh/meh_MEI_jointcall.vcf.gz",
    mea = "/scratch.global/friedlab_MV_MEIs/MEGAnE/merged/meh/meh_MEA_jointcall.vcf.gz"
  output:
    "/scratch.global/friedlab_MV_MEIs/MEGAnE/merged/meh/meh_biallelic.vcf.gz"
  params:
    out_dir = "/scratch.global/friedlab_MV_MEIs/MEGAnE/merged/meh",
    cohort = "meh"
  resources:
    time = 600,
    mem_mb = 60000
  shell:
    '''
      reshape_vcf \
        -i {input.mei} \
        -a {input.mea} \
        -outdir {params.out_dir} \
        -cohort_name {params.cohort}
    '''

