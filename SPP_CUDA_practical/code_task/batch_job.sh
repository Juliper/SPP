#!/bin/bash
#SBATCH -J Cuda-Lab
#SBATCH -e err.%j
#SBATCH -o log.%j
#SBATCH -n 1                  	# Prozesse
#SBATCH --cpus-per-task=1 	# Kerne pro Prozess
#SBATCH --mem-per-cpu=1000    	# Hauptspeicher in MByte pro Rechenkern
#SBATCH -t 00:10:00           	# in Stunden und Minuten, oder '#SBATCH -t 10' - nur Minuten

#SBATCH -A kurs00062
#SBATCH -p kurs00062
#SBATCH --reservation=kurs00062
#SBATCH --gres=gpu

# -------------------------------
# Anschliessend schreiben Sie Ihre eigenen Befehle, wie z.B.
module purge
module load cuda
./CUDATracer 2 "CUDATracing/models/Venus.txt" 1
