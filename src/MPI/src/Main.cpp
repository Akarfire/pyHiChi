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

#include "MpiTransmissionUtilities.h"

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

    if (rank == 0)
    {    
        debugPrintGrid_xyplane(grid, 0);
        //debugPrintGrid_xzplane(grid, 0);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 1)
    {    
        debugPrintGrid_xyplane(grid, 0);
        //debugPrintGrid_xzplane(grid, 0);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 2)
    {    
        debugPrintGrid_xyplane(grid, 0);
        //debugPrintGrid_xplane(grid, 0);
    }

    // Exchange logic

    std::vector<MPI_Datatype> mpi_type_cache;

    auto direction_1 = MpiTransmissionUtilities::Direction::positiveX;
    auto direction_2 = MpiTransmissionUtilities::Direction::negativeX;

    if (rank == 0)
    {
        // Wave 1
        int send_offset;
        MPI_Datatype send_type;
        MpiTransmissionUtilities::defineMpiTransmission_Send(send_offset, send_type, grid->numCells, direction_1, mpi_type_cache);

        MPI_Send(grid->testFieldComponent.getData() + send_offset, 1, send_type, 1, 0, MPI_COMM_WORLD);
        
        // Wave 2
        int recv_offset;
        MPI_Datatype recv_type;
        MpiTransmissionUtilities::defineMpiTransmission_Recv(recv_offset, recv_type, grid->numCells, direction_2, mpi_type_cache);

        MPI_Status status;
        MPI_Recv(grid->testFieldComponent.getData() + recv_offset, 1, recv_type, 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    }

    if (rank == 1)
    {
        // Wave 1
        int recv_offset;
        MPI_Datatype recv_type;
        MpiTransmissionUtilities::defineMpiTransmission_Recv(recv_offset, recv_type, grid->numCells, direction_1, mpi_type_cache);

        MPI_Status status;
        MPI_Recv(grid->testFieldComponent.getData() + recv_offset, 1, recv_type, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        int send_offset_to_2;
        MPI_Datatype send_type_to_2;
        MpiTransmissionUtilities::defineMpiTransmission_Send(send_offset_to_2, send_type_to_2, grid->numCells, direction_1, mpi_type_cache);

        MPI_Send(grid->testFieldComponent.getData() + send_offset_to_2, 1, send_type_to_2, 2, 0, MPI_COMM_WORLD);

        // Wave 2
        int send_offset;
        MPI_Datatype send_type;
        MpiTransmissionUtilities::defineMpiTransmission_Send(send_offset, send_type, grid->numCells, direction_2, mpi_type_cache);

        MPI_Send(grid->testFieldComponent.getData() + send_offset, 1, send_type, 0, 0, MPI_COMM_WORLD);

        int recv_offset_from_2;
        MPI_Datatype recv_type_from_2;
        MpiTransmissionUtilities::defineMpiTransmission_Recv(recv_offset_from_2, recv_type_from_2, grid->numCells, direction_2, mpi_type_cache);

        MPI_Recv(grid->testFieldComponent.getData() + recv_offset_from_2, 1, recv_type_from_2, 2, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    }

    if (rank == 2)
    {
        // Wave 1
        int recv_offset;
        MPI_Datatype recv_type;
        MpiTransmissionUtilities::defineMpiTransmission_Recv(recv_offset, recv_type, grid->numCells, direction_1, mpi_type_cache);

        MPI_Status status;
        MPI_Recv(grid->testFieldComponent.getData() + recv_offset, 1, recv_type, 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        // Wave 2
        int send_offset;
        MPI_Datatype send_type;
        MpiTransmissionUtilities::defineMpiTransmission_Send(send_offset, send_type, grid->numCells, direction_2, mpi_type_cache);

        MPI_Send(grid->testFieldComponent.getData() + send_offset, 1, send_type, 1, 0, MPI_COMM_WORLD);
    }

    MpiTransmissionUtilities::cleanUpMpiTypes(mpi_type_cache);

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
    {    
        debugPrintGrid_xyplane(grid, 0);
        //debugPrintGrid_xzplane(grid, 0);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 1)
    {    
        debugPrintGrid_xyplane(grid, 0);
        //debugPrintGrid_xzplane(grid, 0);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 2)
    {    
        debugPrintGrid_xyplane(grid, 0);
        //debugPrintGrid_xzplane(grid, 0);
    }


    MPI_Finalize();
    return 0;
}