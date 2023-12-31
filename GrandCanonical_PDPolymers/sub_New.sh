#!/bin/bash
#-----------------------------------------------------------------
# Example SLURM job script to run serial applications on Mogon.
#
# This script requests one core (out of 64) on one node. The job
# will have access to all the memory in the node.  Note that this
# job will be charged as if all 64 cores were requested.
#-----------------------------------------------------------------
 
#SBATCH -J PD1.125sc0.8_9.0e-03           # Job name
#SBATCH -o PD1.125sc0.8_9.0e-03.out    # Specify stdout output file (%j expands to jobId)
#SBATCH -p long             # Queue name
#SBATCH -N 1                     # Total number of nodes requested (16 cores/node)
#SBATCH -n 1                    # Total number of tasks
#SBATCH -t 48:00:00              # Run time (hh:mm:ss) - 0.5 hours
 
#SBATCH -A komet331hpc               # Specify allocation to charge against
 
# Load all necessary modules if needed (these are examples)
# Loading modules in the script ensures a consistent environment.
# module load gcc/6.3.0
 
# Launch the executable
./out >&log_PD1.125sc0.8_9.0e_03
