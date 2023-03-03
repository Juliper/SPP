#!/bin/bash
#SBATCH -A kurs00062
#SBATCH -p kurs00062
#SBATCH --reservation=kurs00062
#SBATCH -J daxpy
#SBATCH --mail-type=ALL
#SBATCH -e /home/kurse/kurs00062/je55bela/Praktikum2/test.out.%j
#SBATCH -o /home/kurse/kurs00062/je55bela/Praktikum2/test.out.%j
#SBATCH -n 8
#SBATCH -c 1
#SBATCH --mem-per-cpu=3800
#SBATCH -t 00:10:00

module purge
module load gcc
module load openmpi
module load vampir
module load scorep
cd /home/kurse/kurs00062/je55bela/Praktikum2/code_task/Integral
export SCOREP_ENABLE_TRACING=true
export SCOREP_EXPERIMENT_DIRECTORY=$HPC_SCRATCH/Integral
make
mpirun -n 8 ./integral 100