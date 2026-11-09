#!/bin/bash
#SBATCH -J fdtd_uniform
#SBATCH -t 00:10:00
#SBATCH -o logs/perf_%j.out
#SBATCH -e logs/perf_%j.err 
#SBATCH --cpus-per-task=4

export OMP_NUM_THREADS=4

export I_MPI_THREAD_RUNTIME=openmp
export I_MPI_THREAD_MAX=4
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
srun --cpu_bind=cores ../build/unix_makefiles/src/MPI_PerformanceTests/mpi_perf_tests --grid-longitudinal 128 --grid-transverse 128 --split-mode uniform