import pyHiChi as hichi

from mpi4py import MPI
import pyHiChi as pfc

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

topology = pfc.Topology(pfc.IntVector3d(2, 2, 1), pfc.MPI_LoopType.LoopXY, size)

print(f"[{rank}/{size}] sections={topology.sections}"
      f"size={topology.sections} valid={topology.is_valid_on_this_rank()}")

# Run
# From:
#   pyHiChi\example-tests> 
# Command:
#   mpiexec -n 6 path\to\pyHiChi\.venv\Scripts\python.exe ./example_MPI.py