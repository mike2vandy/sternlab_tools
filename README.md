# sternlab_tools

## af_analysis.py  
This script calculates the number of genotypes and allele frequencies for predefined groups and lists the allele count for each sample provided in the format (0,1,2). 0 is homozgous allele 1, 1 is heterozygous, and 2 is homozygous for allele 2. It also outputs VEP categories, has options for SV and SNV files, and calculates Fisher Exact p-values for dominant and recessive MOIs. Output is `stdout` and `csv` delimited. 

#### Dependencies:  
af_analysis.py requires cyvcf2 and scipy.
```
conda create -n vcf_env -c bioconda -c conda-forge cyvcf2 scipy
conda activate vcf_env
```

#### Usage:
```
usage: af_analysis.py [-h] -v VCF -d DIR [-r REGION] [-c] [-f GROUP1 GROUP2] [-s] [-a]

counts the genotypes and allele frequency for defined groups

options:
  -h, --help            show this help message and exit
  -v VCF, --vcf VCF     Path to VCF file
  -d DIR, --dir DIR     Path to directory containing list files
  -r REGION, --region REGION
                        Optional: region to analyze. Format must be chr or chr:start-end
  -c, --consequences    Find and print all VEP CSQ attributes
  -f GROUP1 GROUP2, --fishers GROUP1 GROUP2
                        Perform Fisher's Exact Test between two groups. Group names must be identical to file names in -d <dir>
  -s, --sv              Use for SV vcf files
  -a, --alt             Genotypes and allele frequencies are oriented toward minor allele. Add this arguement to orient toward
                        alternate allele
```
The script expects files that contain single column lists of samples belonging to each group, without header, in the defined directory. Any number of groups can be defined.
```
lists/
|   hcm.list
|   control.list
|   shelties.list
|   beagles.list
```
Each file name will be used for header prefixes in the output:
```
hcm.count.homref,hcm.count.het,hcm.count.homvar,hcm.count.nocal,hcm.minor_af
```
`-c` is exploratory. It shows what VEP CSQ categories are in the VCF file (if present) and any missing categories the script is looking for.  
`-r` can be used exploring a specific range within the  vcf, i.e. a specific gene region. The format is chrm:start-stop  
Allele frequencies and genotypes are by default oriented toward the **minor** allele. Alternatively, `-a` can be used to orient allele frequencies and genotypes toward the alternative allele.  

If calculating allele frequencies on an SV vcf, use `-s`.  
`-f` Calls the Fisher's Exact test option. It accepts two group names, defined by titles of sample files, space delimited.  
```
af_analysis.py -v file.vcf.gz -d list/ -f hcm control
```

## parallel_af.smk

This snakemake wrapper is used to calcualte allele frequencies in a whole genome vcf among regions in parallel. A few preparations are necessary first.  
It's a little fragile so how I would run this is copy the `parallel_af.smk` and `af_analsis.py` to a new dirctory. 
Within that directory:  
create a directory called `env`
activate the `conda vcf_env` you created above and export that enviornment:  
```
conda activate vcf_env
conda env export > env/vcf_env.yaml
```
You'll need to create a regions file with `splitChr.py` in sternlab_tools like:  
```
splitChr.py /path/to/genome.fasta.fai > dog_regions.txt
```
So your directory structure should look like:

```
.
|-- af_anaylsis.py
|--parallel_af.smk
|--env/
   |--vcf_env.yaml
|lists/
|   hcm.list
|   control.list
|   shelties.list
|   beagles.list
```



