#!/bin/bash
#SBATCH --ntasks=20
#SBATCH --partition=gpu_4
#SBATCH --time=24:00:00
#SBATCH --mem=20gb
#SBATCH --gres=gpu:1
#SBATCH --job-name=EntropyTrainer 
#SBATCH --output=0_GPU_out_%j
#SBATCH --error=0_GPU_err_%j

python3 trainNNMK2.py
