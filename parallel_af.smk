import pandas as pd

regions = list(pd.read_table("dog_regions.txt", header = None)[0])

rule all:
  input:
    "output/final/all_af_analysis.csv"
 
rule process_regions:
  input:
    list_dir = 'lists',
    vcf = '/path/to/vcf_file.vcf.gz'
  output:
    table = "output/{region}/{region}.csv"
  params:
    region = "{region}"
  threads: 1
  conda: "env/vcf_env.yaml"
  resources:
    time = 600,
    mem_mb = 20000
  shell:
    '''
      ./af_analysis.py \
        -v {input.vcf} \
        -d {input.list_dir} \
        -r {params.region} \
      > {output.table}
    '''

rule concat_files:
  input:
    expand("output/{region}/{region}.csv", region=regions)
  output:
    "output/final/all_af_analysis.csv"
  threads: 1
  resources:
    time = 200,
    mem_mb = 20000
  shell:
    '''
      cat {input} |awk 'NR==1 || !/chrm/' > {output}
    '''


