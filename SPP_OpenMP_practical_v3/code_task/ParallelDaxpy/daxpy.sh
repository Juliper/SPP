#!/bin/bash
#SBATCH -A kurs00062
#SBATCH -p kurs00062
#SBATCH --reservation=kurs00062
#SBATCH -J daxpy
#SBATCH --mail-type=ALL
#SBATCH -e /home/kurse/kurs00062/je55bela/Daxpy/Benchmark/test.err.%j
#SBATCH -o /home/kurse/kurs00062/je55bela/Daxpy/Benchmark/test.out.%j
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem-per-cpu=3800
#SBATCH -t 00:10:00

module purge
module load gcc
cd /home/kurse/kurs00062/je55bela/Daxpy/daxpy_normal/
./daxpy 1
./daxpy 2
./daxpy 4
./daxpy 8
./daxpy 16
cd ..
cd daxpy_static/
./daxpy_static 1
./daxpy_static 2
./daxpy_static 4
./daxpy_static 8
./daxpy_static 16
