#!/bin/bash -l
#SBATCH --job-name=sim_job
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=71:59:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH -o OUTPUTS/out/%j.%a.%N.%x.out
#SBATCH -e OUTPUTS/err/%j.%a.%N.%x.err
 
./run2_cylinders.sh
