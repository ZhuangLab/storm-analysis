#!/bin/bash
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 120 # Runtime in minutes
#SBATCH -p serial_requeue # Partition to submit to
#SBATCH --mem=1500 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o log_sa_parallel_%A_%a.out # Standard out goes to this file
#SBATCH -e log_sa_parallel_%A_%a.err # Standard err goes to this filehostname
 
module load python gcc/5.2.0-fasrc01 openmpi/2.0.1-fasrc01 fftw/3.3.5-fasrc01
source activate storm_analysis

cd /n/regal/zhuang_lab/hbabcock/2017_08_24/01_u2os_mc_cell1

MOVIE_ID="11"

rm ./m${MOVIE_ID}_analysis/p_${SLURM_ARRAY_TASK_ID}.hdf5

python -u /n/home12/hbabcock/code/storm-analysis/storm_analysis/multi_plane/multi_plane.py --basename movie_00${MOVIE_ID}_ --bin ./m${MOVIE_ID}_analysis/p_${SLURM_ARRAY_TASK_ID}.hdf5 --xml ./m${MOVIE_ID}_analysis/job_${SLURM_ARRAY_TASK_ID}.xml
