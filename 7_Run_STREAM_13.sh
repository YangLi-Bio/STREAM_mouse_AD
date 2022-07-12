#!/bin/bash
#SBATCH --job-name=Run_STREAM_13+_months
#SBATCH --time=34:20:59
#SBATCH --output=Run_STREAM_13+_months.out
#SBATCH --account=PCON0022
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=300GB
#SBATCH --gpus-per-node=1


set -e


cd /fs/ess/PCON0022/liyang/STREAM/Case_2_AD/Codes/


module load R/4.1.0-gnu9.1
Rscript 7_STREAM_13.R
