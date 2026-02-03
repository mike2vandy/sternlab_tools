#! /bin/bash -l
#SBATCH -p msismall
#SBATCH --time=60:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=60g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mwvandew@ncsu.edu

cd $SLURM_SUBMIT_DIR

conda activate RM

~/software/RepeatMasker/RepeatMasker \
	-div 5 -a -inv -nolow -pa 30 \
	-species dog UU_Cfam_GSD_1.0_ROSY.fa -dir dogRM_Rep3

#conda activate commonTools

#bowtie-build --threads 10 UU_Cfam_GSD_1.0_ROSY.fa UU_Cfam_GSD_1.0_ROSY.fa
