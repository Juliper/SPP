#!/bin/bash
#SBATCH -A kurs00062
#SBATCH -p kurs00062
#SBATCH --reservation=kurs00062
#SBATCH -J HeatedPlate
#SBATCH --mail-type=ALL
#SBATCH -e /home/kurse/kurs00062/ms24quty/SPP/P2/HeatedPlate/out/test.err.%j
#SBATCH -o /home/kurse/kurs00062/ms24quty/SPP/P2/HeatedPlate/out/test.out.%j
#SBATCH -n 16
#SBATCH -c 1
#SBATCH --mem-per-cpu=3800
#SBATCH -t 00:10:00
module purge
module load gcc
module load openmpi
cd /home/kurse/kurs00062/ms24quty/SPP/P2/HeatedPlate/
#mpirun -n 1 heated-plate heated_plate_ref_500.txt
mpirun -n 2 heated-plate-parallel-synchronous #heated_plate_ref_500.txt
mpirun -n 2 heated-plate-parallel-asynchronous #heated_plate_ref_500.txt
mpirun -n 4 heated-plate-parallel-synchronous #heated_plate_ref_500.txt
mpirun -n 4 heated-plate-parallel-asynchronous #heated_plate_ref_500.txt
mpirun -n 8 heated-plate-parallel-synchronous #heated_plate_ref_500.txt
mpirun -n 8 heated-plate-parallel-asynchronous #heated_plate_ref_500.txt
mpirun -n 16 heated-plate-parallel-synchronous #heated_plate_ref_500.txt
mpirun -n 16 heated-plate-parallel-asynchronous #heated_plate_ref_500.txt