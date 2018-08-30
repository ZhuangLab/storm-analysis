
This is a collection of scripts for running storm-analysis, mostly 
object finding, in a distributed computing environment. The scripts
are configured for Harvard's Research Computing cluster, but with
some modification they should work on any cluster that uses SLURM 
for job management.

[slurm](https://slurm.schedmd.com/)

Analysis of a single STORM movie in parallel:
1. create a sub-directory for intermediate results (mkdir mXX_analysis).
2. run slurm/split_analysis_xml.py to create the necessary analysis.xml files.
   Specify the mXX_analysis directory as the working directory.
3. create a sub-directory for the slurm job files, etc. (mkdir mXX_run).
4. cd to this sub-directory.
5. copy parallel.sbatch to this directory and modify as necessary.
6. run analysis with sbatch (using job arrays).
7. check that all the analysis was done properly with slurm/check_analysis.py
8. merge individual localization files with slurm/merge_analysis.py
9. perform the rest of the analysis pipeline on the merged file with storm_analysis/sa_utilities/track_drift_correct.py
