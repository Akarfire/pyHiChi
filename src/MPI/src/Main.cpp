#include <mpi.h>

#include <functional>
#include <string>

#include "Constants.h"
#include "FieldSolver.h"
#include "Grid.h"
#include "PmlFdtd.h"
#include "FieldBoundaryConditionFdtd.h"
#include "FieldGeneratorFdtd.h"
#include "Vectors.h"

#include "MPI_Utilities.h"

using GridType = pfc::YeeGrid;
using FP = pfc::FP;
using FP3 = pfc::FP3;
using Int3 = pfc::Int3;

#define testFieldComponent Ex

FP3 eFunc_0(FP3 coords, FP t) 
{
    return coords;
    // FP3 e;
    // for (int i = 0; i < 3; i++)
    //     e[i] = sin((coords.x + coords.y) * coords.z + i*i);
    // return e;
}

FP3 bFunc_0(FP3 coords, FP t) 
{
    return coords;
    // FP3 e;
    // for (int i = 0; i < 3; i++)
    //     e[i] = sin((coords.x + coords.y) * coords.z + i*i);
    // return e;
}

FP3 eFunc_1(FP3 coords, FP t) 
{
    return coords + FP3(2, 2, 2); //FP3(0, 0, 0);
    // FP3 e;
    // for (int i = 0; i < 3; i++)
    //     e[i] = sin((coords.x + coords.y) * coords.z + i*i);
    // return e;
}

FP3 bFunc_1(FP3 coords, FP t) 
{
    return coords + FP3(2, 2, 2);
    // FP3 e;
    // for (int i = 0; i < 3; i++)
    //     e[i] = sin((coords.x + coords.y) * coords.z + i*i);
    // return e;
}

void initializeGrid(std::unique_ptr<GridType>& grid, int rank) 
{
    Int3 begin = Int3(0, 0, 0);
    Int3 end = grid->numCells;

    std::function<FP3(FP3, FP)> eFunc = eFunc_0;
    std::function<FP3(FP3, FP)> bFunc = bFunc_0;

    if (rank == 1)
    {
        eFunc = eFunc_1;
        bFunc = bFunc_1;
    }

    for (int i = begin.x; i < end.x; i++)
        for (int j = begin.y; j < end.y; j++)
            for (int k = begin.z; k < end.z; k++) {
                grid->Ex(i, j, k) = eFunc(grid->ExPosition(i, j, k), 0).x;
                grid->Ey(i, j, k) = eFunc(grid->EyPosition(i, j, k), 0).y;
                grid->Ez(i, j, k) = eFunc(grid->EzPosition(i, j, k), 0).z;

                grid->Bx(i, j, k) = bFunc(grid->BxPosition(i, j, k), 0).x;
                grid->By(i, j, k) = bFunc(grid->ByPosition(i, j, k), 0).y;
                grid->Bz(i, j, k) = bFunc(grid->BzPosition(i, j, k), 0).z;
            }
}

void debugPrintGrid_xyplane(std::unique_ptr<GridType>& grid, int z)
{
    Int3 size = grid->numCells;

    for (int i = 0; i < size.x; i++)
    {
        std::string line = "";
        for (int j = 0; j < size.y; j++)
            line += std::to_string(grid->testFieldComponent(i, j, z)) + ", ";
        std::cout << line << std::endl;
    }

    std::cout << std::endl << std::endl;
}

void debugPrintGrid_xzplane(std::unique_ptr<GridType>& grid, int y)
{
    Int3 size = grid->numCells;

    for (int i = 0; i < size.z; i++)
    {
        std::string line = "";
        for (int j = 0; j < size.x; j++)
            line += std::to_string(grid->testFieldComponent(j, y, i)) + ", ";
        std::cout << line << std::endl;
    }

    std::cout << std::endl << std::endl;
}

// int main()
// {
//     mpi::Topology topology(mpi::Topology::Type::lineX, mpi::Topology::LoopType::loopXYZ, 3);

//     for (int rank = 0; rank < 3; rank++)
//     {
//         std::cout << "Rank " << rank << std::endl;
//         for (int n = 0; n < 6; n++)
//         {
//             std::cout << n << " : " << topology.getNeighbor(rank, static_cast<mpi::Direction>(n)) << std::endl;
//         }
//     }
// }

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Int3 gridSize = Int3(4, 4, 4);
    FP3 minCoords = FP3(0, 0, 0);
    FP3 maxCoords = FP3(1, 1, 1);
    FP3 gridStep = (maxCoords - minCoords) / (FP3)gridSize;

    std::unique_ptr<GridType> grid;
    grid.reset(new GridType(gridSize, minCoords, gridStep, gridSize));

    initializeGrid(grid, rank);

    // Debug print
    for (int r = 0; r < size; r++)
    {
        if (rank == r)
            debugPrintGrid_xyplane(grid, 0);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Exchange logic
    mpi::Topology topology(mpi::Topology::Type::lineX, mpi::Topology::LoopType::loopX, size);
    mpi::FieldExchanger exchanger(grid->numCells);
    exchanger.PerformExchangeSequence(grid->testFieldComponent.getData(), topology, rank, MPI_COMM_WORLD);

    // Debug print
    for (int r = 0; r < size; r++)
    {
        if (rank == r)
            debugPrintGrid_xyplane(grid, 0);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}