### Run command

~~~
cd <put your path here>/pyHiChi/run-cluster
mkdir -p logs
module load gcc-9.5.0
module load mpi/2021.8.0
sbatch --export=ALL -N 1 -n 1 field_perf_mpi_line.sbatch
sbatch --export=ALL -N 2 -n 2 field_perf_mpi_line.sbatch
sbatch --export=ALL -N 3 -n 3 field_perf_mpi_line.sbatch
sbatch --export=ALL -N 1 -n 1 field_perf_mpi_uniform.sbatch
sbatch --export=ALL -N 2 -n 2 field_perf_mpi_uniform.sbatch
sbatch --export=ALL -N 3 -n 3 field_perf_mpi_uniform.sbatch
sbatch --export=ALL field_perf_mpi_baseline.sbatch
~~~