#! /bin/bash -l
#SBATCH -p msismall
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1
#SBATCH --mem=40g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mwvandew@ncsu.edu

cd $SLURM_SUBMIT_DIR

eval "$(conda shell.bash hook)"
conda activate hts

#intersectBed -a \
#	/home/fried255/fried255/working/pipeline/UU_Cfam_GSD_1.0_ROSY/20241122/joint_call.UU_Cfam_GSD_1.0_ROSY.20241122.vep_exact.vcf.gz \
#	-b UU_Cfam_GSD.exome.bed -header \
#	|bgzip -c > joint_call.UU_Cfam_GSD_1.0_ROSY.20241122.vep_exact.exome.vcf.gz

#bcftools annotate --threads 12 -x "INFO/CSQ"  joint_call.UU_Cfam_GSD_1.0_ROSY.20241122.vep_exact.exome.vcf.gz \
#       | bcftools norm --threads 12 -f /home/fried255/mvandewe/universalDat/UU_Cfam_GSD_ROSY/UU_Cfam_GSD_1.0_ROSY.fa -m - \
#       | bcftools norm --threads 12 -d exact \
#       | bcftools filter --threads 12 -e 'FORMAT/DP < 4 | FORMAT/GQ < 20' --set-GTs . \
#       | bcftools filter --threads 12 -s noPass -g 3 -G 10 -e 'MQ < 40 | QD < 2 | FS > 60 | MQRankSum < -12.4 | ReadPosRankSum < -8.0 |SOR > 3' \
#       |bcftools +fill-tags --threads 12 |bcftools view --threads 12 -i "AC>0" -Oz -o joint_call.UU_Cfam_GSD_1.0_ROSY.20241122.norm..filter.vcf.gz

#plink --vcf joint_call.UU_Cfam_GSD_1.0_ROSY.20241122.norm.filter.vcf.gz --double-id --maf 0.01 --vcf-filter --make-bed --out dogs --allow-extra-chr --chr-set 38

bcftools view --threads 12 -S ^dogToRemove.out joint_call.UU_Cfam_GSD_1.0_ROSY.20241122.norm.filter.vcf.gz \
	|bcftools view --threads 12 -i 'AC>0' \
	|bcftools +fill-tags --threads 12 -Oz -o joint_call.UU_Cfam_GSD_1.0_ROSY.20241122.norm.noBad.filter.vcf.gz

tabix joint_call.UU_Cfam_GSD_1.0_ROSY.20241122.norm.noBad.filter.vcf.gz
