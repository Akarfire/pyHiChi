#!/bin/bash
#SBATCH -J fdtd_uniform
#SBATCH -t 02:00:00
#SBATCH -o logs/perf_%j.out
#SBATCH -e logs/perf_%j.err 

#SBATCH -N 1 
#SBATCH -n 1 

thread_num=1
export OMP_NUM_THREADS=$thread_num
numactl --physcpubind=0-$(($thread_num - 1)) ../build/unix_makefiles/src/MPI_PerformanceTests/reference_perf_test \
    --grid-longitudinal 512 --grid-transverse 512 --log-file "log_explicit_omp${thread_num}_baseline_${SLURM_JOB_ID}.log"
unset OMP_NUM_THREADS

thread_num=2
export OMP_NUM_THREADS=$thread_num
numactl --physcpubind=0-$(($thread_num - 1)) ../build/unix_makefiles/src/MPI_PerformanceTests/reference_perf_test \
    --grid-longitudinal 512 --grid-transverse 512 --log-file "log_explicit_omp${thread_num}_baseline_${SLURM_JOB_ID}.log"
unset OMP_NUM_THREADS

thread_num=4
export OMP_NUM_THREADS=$thread_num
numactl --physcpubind=0-$(($thread_num - 1)) ../build/unix_makefiles/src/MPI_PerformanceTests/reference_perf_test \
    --grid-longitudinal 512 --grid-transverse 512 --log-file "log_explicit_omp${thread_num}_baseline_${SLURM_JOB_ID}.log"
unset OMP_NUM_THREADS

thread_num=5
export OMP_NUM_THREADS=$thread_num
numactl --physcpubind=0-$(($thread_num - 1)) ../build/unix_makefiles/src/MPI_PerformanceTests/reference_perf_test \
    --grid-longitudinal 512 --grid-transverse 512 --log-file "log_explicit_omp${thread_num}_baseline_${SLURM_JOB_ID}.log"
unset OMP_NUM_THREADS

thread_num=8
export OMP_NUM_THREADS=$thread_num
numactl --physcpubind=0-$(($thread_num - 1)) ../build/unix_makefiles/src/MPI_PerformanceTests/reference_perf_test \
    --grid-longitudinal 512 --grid-transverse 512 --log-file "log_explicit_omp${thread_num}_baseline_${SLURM_JOB_ID}.log"
unset OMP_NUM_THREADS

thread_num=10
export OMP_NUM_THREADS=$thread_num
numactl --physcpubind=0-$(($thread_num - 1)) ../build/unix_makefiles/src/MPI_PerformanceTests/reference_perf_test \
    --grid-longitudinal 512 --grid-transverse 512 --log-file "log_explicit_omp${thread_num}_baseline_${SLURM_JOB_ID}.log"
unset OMP_NUM_THREADS

thread_num=16
export OMP_NUM_THREADS=$thread_num
numactl --physcpubind=0-$(($thread_num - 1)) ../build/unix_makefiles/src/MPI_PerformanceTests/reference_perf_test \
    --grid-longitudinal 512 --grid-transverse 512 --log-file "log_explicit_omp${thread_num}_baseline_${SLURM_JOB_ID}.log"
unset OMP_NUM_THREADS

thread_num=20
export OMP_NUM_THREADS=$thread_num
numactl --physcpubind=0-$(($thread_num - 1)) ../build/unix_makefiles/src/MPI_PerformanceTests/reference_perf_test \
    --grid-longitudinal 512 --grid-transverse 512 --log-file "log_explicit_omp${thread_num}_baseline_${SLURM_JOB_ID}.log"
unset OMP_NUM_THREADS