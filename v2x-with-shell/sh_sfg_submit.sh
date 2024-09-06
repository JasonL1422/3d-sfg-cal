#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=64GB
#SBATCH --time=24:00:00
#SBATCH --account=open
#SBATCH --partition=open
#SBATCH --job-name=3d_sfg_calc
echo "==============================="
echo " 3d_sfg_calc "
echo "==============================="

python sfg-cals-v2x.py
