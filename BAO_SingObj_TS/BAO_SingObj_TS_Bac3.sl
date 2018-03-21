#!/bin/bash
#SBATCH -J JobArray
#SBATCH --time=96:00:00     # Walltime
#SBATCH -A uoa00285         # Project Account
#SBATCH --ntasks=1          # number of tasks
#SBATCH --mem-per-cpu=8G  # memory/cpu (in MB)
#SBATCH --cpus-per-task=8   # 10 OpenMP Threads
ml Java/1.8.0_5
ml METIS/5.1.0-intel-2015a-shared
ml OpenBLAS/0.2.13-GCC-4.9.2-LAPACK-3.5.0
/projects/uoa00285/AMPL/bin/ampl_lic start
srun java -jar ./dist/BAO_SingObj_TS.jar $SLURM_JOBID inputFileBac3.txt
