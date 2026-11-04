#!/bin/bash
#SBATCH -J fdtd_uniform
#SBATCH -t 00:10:00
#SBATCH -o logs/perf_%j.out
#SBATCH -e logs/perf_%j.err 

source /etc/profile.d/modules.sh
module load gcc-9.5.0
module load mpi/2021.8.0

srun --cpu_bind=cores ../build/unix_makefiles/src/MPI_PerformanceTests/mpi_perf_tests --grid-longitudinal 128 --grid-transverse 128 --split-mode uniform