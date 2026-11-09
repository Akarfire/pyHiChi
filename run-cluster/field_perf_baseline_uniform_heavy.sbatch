#!/bin/bash
#SBATCH -J fdtd_uniform
#SBATCH -t 02:00:00
#SBATCH -o logs/perf_%j.out
#SBATCH -e logs/perf_%j.err 
#SBATCH --cpus-per-task=4

export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

export OMP_NUM_THREADS=1
srun --cpu_bind=cores ../build/unix_makefiles/src/MPI_PerformanceTests/reference_perf_test --grid-longitudinal 512 --grid-transverse 512 --split-mode uniform
mv BASELINE_perf_log_L512_T512_uniform_sec1x1x1_rank0.log log_omp1_baseline_${SLURM_JOB_ID}.log

export OMP_NUM_THREADS=4
srun --cpu_bind=cores ../build/unix_makefiles/src/MPI_PerformanceTests/reference_perf_test --grid-longitudinal 512 --grid-transverse 512 --split-mode uniform
mv BASELINE_perf_log_L512_T512_uniform_sec1x1x1_rank0.log log_omp4_baseline_${SLURM_JOB_ID}.log