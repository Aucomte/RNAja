#!/bin/bash

############      SLURM CONFIGURATION      ###################
#SBATCH --partition=highmemdell
#SBATCH --nodelist=node30
#SBATCH --job-name=RNAja
#SBATCH --output=RNAja-out
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --time=120:00:00

############################################################

# sbatch launch_RNAja.sh config.yaml

# module load

module load system/python/3.7.2
module load system/singularity/3.6.0
#module load system/Miniconda3/1.0

RNAja run_local -c ${CONFIG} -t 20
