#! /bin/bash -l
#SBATCH -p msismall
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mwvandew@ncsu.edu

cd $SLURM_SUBMIT_DIR

#conda activate vep 


#vep -i meh_pass.af.sorted.vcf.gz \
#	-o meh.vepped.vcf \
#	--gtf ~/universalDat/UU_Cfam_GSD_ROSY/UU_Cfam_GSD_1.0_ROSY.refSeq.ensformat.gtf.gz \
#	--fasta ~/universalDat/UU_Cfam_GSD_ROSY/UU_Cfam_GSD_1.0_ROSY.fa \
#	--vcf \
#	--fork 10 \
#	--dont_skip

conda activate commonTools

#./something.py meh.vepped.vcf > uniques.FriedDogs.tab

#./perDogMEI.py meh.vepped.vcf > dog.insCt.csv
#./getAF.py meh_pass.af.sorted.vcf.gz > mei.afs.csv
./gatherInfo.py meh.vepped.vcf.gz |uniq > dogs.insert.properties.csv
