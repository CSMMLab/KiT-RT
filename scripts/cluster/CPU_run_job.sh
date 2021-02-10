#!/bin/bash
#BATCH --ntasks=24
#SBATCH --time=24:00:00
#SBATCH --mem=20gb
#SBATCH --partition=single
#SBATCH --job-name=EntropyTrainerCPU
#SBATCH --output=0_CPU_training_out_%j
#SBATCH --error=0_CPU_training_err_%j

python3 trainNNMK2.py
