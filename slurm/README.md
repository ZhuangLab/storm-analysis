
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
4.