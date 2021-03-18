#! /bin/bash
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH --mem=50000  ## memory allotted for the job in MB

# This is to get e-mail notifications
# when the jobs start and end
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=

# Name of the job as it appears in squeue
#SBATCH --job-name=rand_summary

#! DELLA
# Load python
module load python 
module load anaconda

# The first argument is the total number of cores you
# want but we get it from a SLURM environment variable
python random_summary_files_TS.py