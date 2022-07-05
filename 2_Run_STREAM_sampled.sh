#!/bin/bash
#SBATCH --job-name=Run_STREAM_sampled
#SBATCH --time=03:20:59
#SBATCH --output=Run_STREAM_sampled.out
#SBATCH --account=PCON0022
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=100GB
#SBATCH --gpus-per-node=1


set -e


cd /fs/ess/PCON0022/liyang/STREAM/Case_2_AD/Codes/


module load R/4.1.0-gnu9.1
Rscript 2_Test_STREAM.R
