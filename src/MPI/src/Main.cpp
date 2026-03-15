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

using GridType = pfc::YeeGrid;
using FP = pfc::FP;
using FP3 = pfc::FP3;
using Int3 = pfc::Int3;

#define testFieldComponent Ey

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
    return coords + FP3(1, 1, 1); //FP3(0, 0, 0);
    // FP3 e;
    // for (int i = 0; i < 3; i++)
    //     e[i] = sin((coords.x + coords.y) * coords.z + i*i);
    // return e;
}

FP3 bFunc_1(FP3 coords, FP t) 
{
    return coords + FP3(1, 1, 1);
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

void define_line_type(MPI_Datatype &out_type, int count, int stride, int block_size = 1)
{
    MPI_Type_vector(count, block_size, stride, MPI_DOUBLE, &out_type);
    MPI_Type_commit(&out_type);
}

void define_plane_type(MPI_Datatype &out_type, const MPI_Datatype &line_type, int count, int stride)
{
    // Calculating stride in bytes
    MPI_Aint byte_stride;
    MPI_Type_extent(MPI_DOUBLE, &byte_stride);
    byte_stride *= stride;
    
    // Creating a vector type with a byte stride
    MPI_Type_create_hvector(count, 1, byte_stride, line_type, &out_type);
    MPI_Type_commit(&out_type);
}

int main(int argc, char** argv)
{
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

    initializeGrid(grid, rank);

    if (rank == 1)
    {    
        debugPrintGrid_xyplane(grid, 0);
        debugPrintGrid_xyplane(grid, 1);
        debugPrintGrid_xyplane(grid, 2);
    }

    // MPI_Datatype X_lineType;
    // define_line_type(X_lineType, grid->numCells.x, grid->numCells.y * grid->numCells.z, 1);

    // MPI_Datatype XZ_planeType;
    // define_plane_type(XZ_planeType, X_lineType, grid->numCells.z, 1);

    // if (rank == 0)
    //     MPI_Send(grid->testFieldComponent.getData(), 1, XZ_planeType, 1, 0, MPI_COMM_WORLD);
    // if (rank == 1)
    // {
    //     MPI_Status status;
    //     MPI_Recv(grid->testFieldComponent.getData(), 1, XZ_planeType, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    //     std::cout << "MPI SEND ---------------------------------------------- MPI RECV" << std::endl << std::endl << std::endl;

    //     debugPrintGrid_xyplane(grid, 0);
    //     debugPrintGrid_xyplane(grid, 1);
    //     debugPrintGrid_xyplane(grid, 2);
    // }

    // MPI_Barrier(MPI_COMM_WORLD);

    MPI_Datatype Y_lineType;
    define_line_type(Y_lineType, grid->numCells.y, grid->numCells.z, 1);

    MPI_Datatype YZ_planeType;
    define_plane_type(YZ_planeType, Y_lineType, grid->numCells.z, 1);

    if (rank == 0)
        MPI_Send(grid->testFieldComponent.getData(), 1, YZ_planeType, 1, 0, MPI_COMM_WORLD);
    if (rank == 1)
    {
        MPI_Status status;
        int offset = grid->numCells.x * grid->numCells.y * grid->numCells.z - grid->numCells.y * grid->numCells.z;
        MPI_Recv(grid->testFieldComponent.getData() + offset, 1, YZ_planeType, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        std::cout << "MPI SEND ---------------------------------------------- MPI RECV" << std::endl << std::endl << std::endl;

        debugPrintGrid_xyplane(grid, 0);
        debugPrintGrid_xyplane(grid, 1);
        debugPrintGrid_xyplane(grid, 2);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // MPI_Datatype XY_planeType;
    // define_plane_type(XY_planeType, X_lineType, grid->numCells.y, grid->numCells.z);

    // if (rank == 0)
    //     MPI_Send(grid->testFieldComponent.getData(), 1, XY_planeType, 1, 0, MPI_COMM_WORLD);
    // if (rank == 1)
    // {
    //     MPI_Status status;
    //     MPI_Recv(grid->testFieldComponent.getData(), 1, XY_planeType, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    //     std::cout << "MPI SEND ---------------------------------------------- MPI RECV" << std::endl << std::endl << std::endl;

    //     debugPrintGrid_xyplane(grid, 0);
    //     debugPrintGrid_xyplane(grid, 1);
    //     debugPrintGrid_xyplane(grid, 2);
    // }

    // // Clean up
    // MPI_Type_free(&X_lineType);
    // MPI_Type_free(&XZ_planeType);
    // MPI_Type_free(&Y_lineType);
    // MPI_Type_free(&YZ_planeType);

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