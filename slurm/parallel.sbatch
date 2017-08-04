#!/bin/bash
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 120 # Runtime in minutes
#SBATCH -p serial_requeue # Partition to submit to
#SBATCH --mem=600 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o log_sa_parallel_%A_%a.out # Standard out goes to this file
#SBATCH -e log_sa_parallel_%A_%a.err # Standard err goes to this filehostname

module load python gcc/5.2.0-fasrc01 openmpi/2.0.1-fasrc01 fftw/3.3.5-fasrc01
source activate storm_analysis

cd /n/regal/zhuang_lab/hbabcock/2017_08_03/02_hbaic1_2017_06_09_02_biplane

rm ./m2_analysis/p_${SLURM_ARRAY_TASK_ID}_mlist.bin

python -u /n/home12/hbabcock/code/storm-analysis/storm_analysis/multi_plane/multi_plane.py --basename movie_0002_ --bin ./m2_analysis/p_${SLURM_ARRAY_TASK_ID}_mlist.bin --xml ./m2_analysis/job_${SLURM_ARRAY_TASK_ID}.xml
