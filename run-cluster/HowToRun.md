### Run command

~~~
cd <put your path here>/pyHiChi/run-cluster
mkdir -p logs
sbatch -N 1 -n 1 field_perf_mpi_line.sbatch
sbatch -N 2 -n 2 field_perf_mpi_line.sbatch
sbatch -N 3 -n 3 field_perf_mpi_line.sbatch
sbatch -N 1 -n 1 field_perf_mpi_uniform.sbatch
sbatch -N 2 -n 2 field_perf_mpi_uniform.sbatch
sbatch -N 3 -n 3 field_perf_mpi_uniform.sbatch
sbatch field_perf_mpi_baseline.sbatch
~~~