
singularity: config['sif']

rule all:
    input:
        "final_gather/colony_cats.vep.vcf.gz",

rule vep_scatter_intervals:
    output:
        acgt_ivals = "intervals/vep/acgt.interval_list"
    params:
        ref_fasta = config['ref_fasta'],
        contig_ns = config['nrun_length'],
    threads: 1
    resources:
         time   = 20,
         mem_mb = 8000
    shell:
        '''
            java -jar /opt/wags/src/picard.jar \
                ScatterIntervalsByNs \
                R={params.ref_fasta} \
                OT=ACGT \
                N={params.contig_ns} \
                O={output.acgt_ivals}

            sed -i '/^chrUn/d' {output.acgt_ivals}

	'''

checkpoint split_intervals:
    input:
        acgt_ivals = "intervals/vep/acgt.interval_list"
    output:
        directory("intervals/vep/scattered")
    params:
        ref_fasta    = config['ref_fasta'],
        scatter_size = config['scatter_size'],
    threads: 1
    resources:
         time   = 20,
         mem_mb = 8000
    shell:
        '''
            gatk SplitIntervals \
                -R {params.ref_fasta} \
                -L {input.acgt_ivals} \
                --scatter-count {params.scatter_size} \
                --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \
                -O {output}
        '''

rule vep_by_interval:
    input:
        final_vcf = "/share/stern/mwvandew/colony_cats/colony_cats.vcf.gz",
        final_tbi = "/share/stern/mwvandew/colony_cats/colony_cats.vcf.gz.tbi",
        interval  = "intervals/vep/scattered/{vep_interval}-scattered.interval_list"
    output:
        final_interval    = "final_gather/split/vep_{vep_interval}/joint_call.{vep_interval}.vcf.gz",
        interval_vep      = "final_gather/vep/vep_{vep_interval}/joint_call.{vep_interval}.vep.vcf.gz",
        interval_vep_tbi  = "final_gather/vep/vep_{vep_interval}/joint_call.{vep_interval}.vep.vcf.gz.tbi",
        interval_vep_html = "final_gather/vep/vep_{vep_interval}/joint_call.{vep_interval}.vep.vcf_summary.html",
    params:
        out_name = lambda wildcards, output: os.path.splitext(output.interval_vep)[0],
        ref_fasta = config["ref_fasta"],
        ref_gtf   = config["ref_gtf"]
    threads: 6
    resources:
         time   = 720,
         mem_mb = 60000
    shell:
        '''
            set -e

            source activate ensembl-vep

            gatk --java-options "-Xmx12g -Xms3g" \
                SelectVariants \
                -V {input.final_vcf} \
                -O {output.final_interval} \
                -L {input.interval}

            vep \
                -i {output.final_interval} \
                -o {params.out_name} \
                --gtf {params.ref_gtf} \
                --fasta {params.ref_fasta} \
                --fork {threads} \
                --everything \
                --force_overwrite \
                --vcf \
                --dont_skip

            bgzip --threads {threads} -c {params.out_name} > {output.interval_vep}
            tabix -p vcf {output.interval_vep}
        '''

def get_vep_vcfs(wildcards):
    # interval dir from split intervals
    ivals_dir = checkpoints.split_intervals.get(**wildcards).output[0]
    # variable number of intervals
    INTERVALS, = glob_wildcards(os.path.join(ivals_dir,"{vep_interval}-scattered.interval_list"))
    # return list of recal vcfs
    return sorted(expand(
        "final_gather/vep/vep_{vep_interval}/joint_call.{vep_interval}.vep.vcf.gz",
        vep_interval=INTERVALS
    ))

rule final_gather_veps:
    input:
        get_vep_vcfs
    output:
        vep_vcf = "final_gather/colony_cats.vep.vcf.gz",
        vep_tbi = "final_gather/colony_cats.vep.vcf.gz.tbi",
    params:
        vcf_tmp = "final_gather/joint_genotype.TMP.gz",
        veps    = lambda wildcards, input: " --input ".join(map(str,input)),
    threads: 12
    resources:
         time   = 770,
         mem_mb = 22000
    shell:
        '''
            set -e

            gatk --java-options "-Xmx18g -Xms6g" \
                GatherVcfsCloud \
                --ignore-safety-checks \
                --gather-type BLOCK \
                --input {params.veps} \
                --output {output.vep_vcf}

            tabix -p vcf {output.vep_vcf} 
        '''
