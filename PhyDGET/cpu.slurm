#!/bin/bash
#SBATCH --job-name=phy
#SBATCH --partition=cpu
#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --time=99:20:00


python3 /share/home/zhou3lab/kangwen/01.Digit/05.PhyDeg/phydget-main/phydget.py emu.job
