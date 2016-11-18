#!/bin/bash

#SBATCH --job-name="funbuts coverage"                # name of the job submitted
#SBATCH -p short                                        # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 40                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 02:00:00                                     # time allocated for this job hours:mins:seconds
#SBATCH --mail-user=julian.trachsel@ars.usda.gov        # enter your email address to receive emails
#SBATCH --mail-user=julestrachsel@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL                      # will receive an email when job starts, ends or fails
#SBATCH -o "stdout.%j.%N"                               # standard out %j adds job number to outputfile name and %N adds
the node name
#SBATCH -e "stderr.%j.%N"                               # optional but it prints our standard error

# ENTER COMMANDS HERE:

module load python_3
module load muscle
module load openblas

python3 butcoverage.py finalbuts.fasta >> log.txt

#End of file
