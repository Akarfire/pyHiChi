#!/bin/bash
#SBATCH -J fdtd_uniform
#SBATCH -t 02:00:00
#SBATCH -o logs/perf_%j.out
#SBATCH -e logs/perf_%j.err 

thread_num=${1:-1}

export OMP_NUM_THREADS=$thread_num
numactl --physcpubind=0-$(($thread_num - 1)) ../build/unix_makefiles/src/MPI_PerformanceTests/reference_perf_test --grid-longitudinal 512 --grid-transverse 512
mv BASELINE_perf_log_L512_T512_uniform_sec1x1x1_rank0.log log_omp${thread_num}_baseline_${SLURM_JOB_ID}.log
unset OMP_NUM_THREADS