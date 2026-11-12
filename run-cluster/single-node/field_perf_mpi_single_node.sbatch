#!/bin/bash
#SBATCH -J fdtd_uniform_single_node
#SBATCH -t 02:00:00
#SBATCH -o logs/perf_%j.out
#SBATCH -e logs/perf_%j.err 
#SBATCH -N 1

proc_num=${1:-1}
thread_num=${2:-2}
subdir=${3:-}

if [ -n "$subdir" ]; then
    log_dir="$subdir"
else
    log_dir="."
fi

mkdir -p "$log_dir"
log_file="${log_dir}/log_${proc_num}p_${thread_num}t_${SLURM_JOB_ID}.log"

export OMP_NUM_THREADS=$thread_num

srun --ntasks=$proc_num \
     --cpus-per-task=$thread_num \
     --cpu_bind=cores \
     numactl --localalloc \
     ../../build/unix_makefiles/src/MPI_PerformanceTests/mpi_perf_tests \
     --grid-longitudinal 512 \
     --grid-transverse 512 \
     --split-mode uniform \
     --log-file "$log_file"

unset OMP_NUM_THREADS