#!/bin/bash

# Solve equilibrium problem for integrated h2-electricity system with ADMM

#SBATCH --job-name=H2_market
#SBATCH --partition=compute
#SBATCH --account=Education-EEMCS-MSc-SET
#SBATCH --time=09:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=9
#SBATCH --mem-per-cpu=2GB

module load 2022r2
srun ~/julia-1.9.0/bin/julia --threads=9 MAIN.jl
