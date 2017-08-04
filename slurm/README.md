
This is a collection of scripts for running storm-analysis, mostly 
object finding, in a distributed computing environment. The scripts
are configured for Harvard's Research Computing cluster, but with
some modification they should work on any cluster that uses SLURM 
for job management.

[slurm](https://slurm.schedmd.com/)

Analysis of a single STORM movie in parallel:
1. create a sub-directory for intermediate results.
2. cd to the sub-directory.
3. run slurm/split_analysis.xml to create the necessary analysis.xml files.
4. modify parallel.sbatch as necessary.
5. run analysis with sbatch (using job arrays).
6. check that all the analysis was done properly with slurm/check_analysis.py
7. merge individual localization files with slurm/merge_analysis.py
8. perform the rest of the analysis pipeline on the merged file with storm_analysis/sa_utilities/track_average_correct.py
