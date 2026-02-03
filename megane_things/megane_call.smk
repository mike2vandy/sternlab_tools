import pandas as pd

singularity: "/home/fried255/mvandewe/software/megane/MEGAnE.sif"

tmp = pd.read_csv("rnd4.dogs.list", header = None)

dogs = list(tmp[0])

rule all:
  input:
    expand("/scratch.global/friedlab_MV_MEIs/MEGAnE/indCall/{dogID}/MEI_final_percentile_genotyped.vcf", dogID = dogs)

rule callMEI:
  input:
    cram = "/scratch.global/friedlab_MV_MEIs/crams/{dogID}.UU_Cfam_GSD_1.0_ROSY.cram"
  output:
    vcf = "/scratch.global/friedlab_MV_MEIs/MEGAnE/indCall/{dogID}/MEI_final_percentile_genotyped.vcf" 
  params:
    genome = "/home/fried255/mvandewe/universalDat/UU_Cfam_GSD_ROSY/UU_Cfam_GSD_1.0_ROSY.fa",
    kmer = "/home/fried255/mvandewe/software/megane/UU_Cfam/kmer/UU_Cfam_GSD_1.0_ROSY.fa.mk",
    lib = "/home/fried255/mvandewe/software/megane/UU_Cfam/library.fa",
    rm_out = "/home/fried255/mvandewe/software/megane/UU_Cfam/dogLSine.out",
    norep = "/home/fried255/mvandewe/software/megane/UU_Cfam/norep.txt",
    polyA = "/home/fried255/mvandewe/software/megane/UU_Cfam/withPa.txt",
    chrlist = "/home/fried255/mvandewe/software/megane/UU_Cfam/chrlist.out",
    outdir = "/scratch.global/friedlab_MV_MEIs/MEGAnE/indCall/{dogID}",
    sample_name = "{dogID}"
  threads: 8
  resources:
    time = 1440,
    mem_mb = 90000
  shell:
    '''
      call_genotype \
        -i {input.cram} \
        -fa {params.genome} \
        -fadb {params.genome} \
        -mk {params.kmer} \
        -rep {params.lib} \
        -repout {params.rm_out} \
        -repremove {params.norep} \
        -pA_ME {params.polyA} \
        -mainchr {params.chrlist} \
        -outdir {params.outdir} \
        -sample_name {params.sample_name} \
        -p {threads}
    '''
