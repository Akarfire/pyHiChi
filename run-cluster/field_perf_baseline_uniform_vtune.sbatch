#!/bin/bash
#SBATCH -J fdtd_vtune
#SBATCH -t 02:00:00
#SBATCH -o logs/vtune_%j.out
#SBATCH -e logs/vtune_%j.err 

thread_num=${1:-1}

export OMP_NUM_THREADS=$thread_num
numactl --physcpubind=0-$(($thread_num - 1)) \
    vtune -collect hpc-performance -result-dir vtune_hpcperformance_baseline_omp${thread_num} -- \
    ../build/unix_makefiles/src/MPI_PerformanceTests/reference_perf_test \
    --grid-longitudinal 512 --grid-transverse 512 --log-file "log_omp${thread_num}_baseline_${SLURM_JOB_ID}.log"
unset OMP_NUM_THREADS