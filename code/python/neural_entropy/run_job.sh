#!/bin/bash
#BATCH --ntasks=24
#SBATCH --partition=gpu_4
#SBATCH --time=24:00:00
#SBATCH --mem=20gb
#SBATCH --gres=gpu:1
#SBATCH --job-name=EntropyTrainerCPU

python3 trainNNMK2.py
