#include <mpi.h>

#include "Constants.h"
#include "FieldSolver.h"
#include "Grid.h"
#include "PmlFdtd.h"
#include "FieldBoundaryConditionFdtd.h"
#include "FieldGeneratorFdtd.h"
#include "Vectors.h"

using GridType = pfc::YeeGrid;
namespace pfc {

 FP3 eFunc(FP3 coords, FP t) 
 {
    FP3 e;
    for (int i = 0; i < 3; i++)
        e[i] = sin((coords.x + coords.y) * coords.z + i*i);
    return e;
}

FP3 bFunc(FP3 coords, FP t) 
{
    FP3 e;
    for (int i = 0; i < 3; i++)
        e[i] = sin((coords.x + coords.y) * coords.z + i*i);
    return e;
}

void initializeGrid(std::unique_ptr<GridType>& grid, const FP3& coordOffset) 
{
    Int3 begin = Int3(0, 0, 0);
    Int3 end = grid->numCells;

    for (int i = begin.x; i < end.x; i++)
        for (int j = begin.y; j < end.y; j++)
            for (int k = begin.z; k < end.z; k++) {
                grid->Ex(i, j, k) = eFunc(grid->ExPosition(i, j, k) + coordOffset, 0).x;
                grid->Ey(i, j, k) = eFunc(grid->EyPosition(i, j, k) + coordOffset, 0).y;
                grid->Ez(i, j, k) = eFunc(grid->EzPosition(i, j, k) + coordOffset, 0).z;

                grid->Bx(i, j, k) = bFunc(grid->BxPosition(i, j, k) + coordOffset, 0).x;
                grid->By(i, j, k) = bFunc(grid->ByPosition(i, j, k) + coordOffset, 0).y;
                grid->Bz(i, j, k) = bFunc(grid->BzPosition(i, j, k) + coordOffset, 0).z;
            }
}

}

int main(int argc, char** argv)
{
    using FP3 = pfc::FP3;
    using Int3 = pfc::Int3;

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Int3 gridSize = Int3(8, 8, 8);
    FP3 minCoords = FP3(0, 0, 0);
    FP3 maxCoords = FP3(1, 1, 1);
    FP3 gridStep = (maxCoords - minCoords) / (FP3)gridSize;

    std::unique_ptr<GridType> grid;
    grid.reset(new GridType(gridSize, minCoords, gridStep, gridSize));

    FP3 coordOffset = FP3(0.0, 0.0, 0.0);

    if (rank == 1) coordOffset += FP3(1.0, 0.0, 0.0);
    else if (rank == 2) coordOffset += FP3(0.0, 1.0, 0.0);
    else if (rank == 3) coordOffset += FP3(1.0, 1.0, 0.0);

    pfc::initializeGrid(grid, coordOffset);

    std::cout << std::endl << "Hello Grid! From rank:" << rank << " : " << grid->getB(pfc::FP3(0.5, 0.5, 0.5)) << std::endl;

    MPI_Finalize();
    return 0;
}


// int main(int argc, char** argv) {
//     MPI_Init(&argc, &argv);

//     int rank, size;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &size);
    
//     std::cout << "Hello Grid from rank " << rank << " of " << size << std::endl;
    
//     MPI_Finalize();
//     return 0;
// }