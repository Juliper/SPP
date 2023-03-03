#!/bin/bash
#SBATCH -A kurs00062
#SBATCH -p kurs00062
#SBATCH --reservation=kurs00062
#SBATCH -J hello
#SBATCH --mail-type=ALL
#SBATCH -e /home/kurse/kurs00062/ah93degi/Batch/test.err.%j
#SBATCH -o /home/kurse/kurs00062/ah93degi/Batch/test.out.%j
#SBATCH -n 1
#SBATCH --mem-per-cpu=3800
#SBATCH -t 00:30:00
#SBATCH -c 16

module purge 
module load gcc/10
cd /home/kurse/kurs00062/ah93degi/ParallelNbody

# ---- Test für 500 Körper (und 50.000 Zeitschritte) ----
./bin/spp-nbody-omp /home/kurse/kurs00062/ah93degi/ParallelNbody/scenarios/random500.txt

# ---- Test für 2.000 Körper (und 10.000 Zeitschritte) ----
./bin/spp-nbody-omp /home/kurse/kurs00062/ah93degi/ParallelNbody/scenarios/random2000.txt

# ---- Test für 5.000 Körper (und 1.000 Zeitschritte) ----
./bin/spp-nbody-omp /home/kurse/kurs00062/ah93degi/ParallelNbody/scenarios/random5000.txt





