#include "MPI_Utilities.h"

namespace mpi
{

// Calls MPI_Type_free on all types stored in mpi_type_cache
void cleanUpMpiTypes(std::vector<MPI_Datatype>& mpi_type_cache)
{
    for (auto& type : mpi_type_cache)
        MPI_Type_free(&type);
}

// Defines an mpi data type for transmitting sub arrays of the grid
void MPI_FieldUtils::defineSubArrayType(MPI_Datatype &out_type, int sizes[3], int sub_sizes[3], int starts[3], std::vector<MPI_Datatype>& mpi_type_cache)
{
    MPI_Type_create_subarray(3, sizes, sub_sizes, starts, MPI_ORDER_C, MPI_DOUBLE, &out_type);
    MPI_Type_commit(&out_type);
    mpi_type_cache.push_back(out_type);
}

// Resolving "send" operation parameters based on direction
void MPI_FieldUtils::resolveSendParameters(int out_sizes[3], int out_sub_sizes[3], int out_starts[3], 
        const pfc::Int3& num_cells, Direction direction)
{
    out_sizes[0] = num_cells.x;
    out_sizes[1] = num_cells.y;
    out_sizes[2] = num_cells.z;

    // General outputs (to avoid code repettition)
    out_sub_sizes[0] = num_cells.x;
    out_sub_sizes[1] = num_cells.y;
    out_sub_sizes[2] = num_cells.z;
    out_starts[0] = 0;
    out_starts[1] = 0;
    out_starts[2] = 0;

    // Specifying outputs based on direction
    switch (direction)
    {
    case Direction::positiveX:
        // Second to last YZ plane
        out_sub_sizes[0] = 1;
        out_starts[0] = num_cells.x - 2;
        break;

    case Direction::negativeX:
        // Second YZ plane
        out_sub_sizes[0] = 1;
        out_starts[0] = 1;
        break;

    case Direction::positiveY:
        // Second to last XZ plane
        out_sub_sizes[1] = 1;
        out_starts[1] = num_cells.y - 2;
        break;

    case Direction::negativeY:
        // Second XZ plane
        out_sub_sizes[1] = 1;
        out_starts[1] = 1;
        break;

    case Direction::positiveZ:
        // Second to last XY plane
        out_sub_sizes[2] = 1;
        out_starts[2] = num_cells.z - 2;
        break;

    case Direction::negativeZ:
        // Second XY plane
        out_sub_sizes[2] = 1;
        out_starts[2] = 1;
        break;

    default: break;
    }
}

// Resolving "recv" operation parameters based on direction
void MPI_FieldUtils::resolveRecvParameters(int out_sizes[3], int out_sub_sizes[3], int out_starts[3], 
        const pfc::Int3& num_cells, Direction direction)
{
    out_sizes[0] = num_cells.x;
    out_sizes[1] = num_cells.y;
    out_sizes[2] = num_cells.z;

    // General outputs (to avoid code repettition)
    out_sub_sizes[0] = num_cells.x;
    out_sub_sizes[1] = num_cells.y;
    out_sub_sizes[2] = num_cells.z;
    out_starts[0] = 0;
    out_starts[1] = 0;
    out_starts[2] = 0;

    // Specifying outputs based on direction
    // IMPORTANT: Direction is specified relative to the SENDER
    switch (direction)
    {
    case Direction::positiveX:
        // First YZ plane
        out_sub_sizes[0] = 1;
        out_starts[0] = 0;
        break;

    case Direction::negativeX:
        // Last YZ plane
        out_sub_sizes[0] = 1;
        out_starts[0] = num_cells.x - 1;
        break;

    case Direction::positiveY:
        // First XZ plane
        out_sub_sizes[1] = 1;
        out_starts[1] = 0;
        break;

    case Direction::negativeY:
        // Last XZ plane
        out_sub_sizes[1] = 1;
        out_starts[1] = num_cells.y - 1;
        break;

    case Direction::positiveZ:
        // First XY plane
        out_sub_sizes[2] = 1;
        out_starts[2] = 0;
        break;

    case Direction::negativeZ:
        // Last XY plane
        out_sub_sizes[2] = 1;
        out_starts[2] = num_cells.z - 1;
        break;

    default: break;
    }
}


// Determines mpi data type and required offset to send field data in the required direction
// All created types are stored into mpi_type_cache, for it to be cleaned up later (using cleanUpMpiTypes)
void MPI_FieldUtils::defineTransmission_Send(MPI_Datatype& out_type, pfc::Int3 grid_num_cells, Direction direction,
    std::vector<MPI_Datatype>& mpi_type_cache)
{
    // Parameters
    int sizes[3];
    int sub_sizes[3];
    int starts[3];

    // Resolving parameters
    resolveSendParameters(sizes, sub_sizes, starts, grid_num_cells, direction);

    // Creating type
    defineSubArrayType(out_type, sizes, sub_sizes, starts, mpi_type_cache);
}

// Determines mpi data type and required offset to received field data from the required direction (direction relative to the sender)
// All created types are stored into mpi_type_cache, for it to be cleaned up later (using cleanUpMpiTypes)
void MPI_FieldUtils::defineTransmission_Recv(MPI_Datatype& out_type, pfc::Int3 grid_num_cells, Direction direction,
    std::vector<MPI_Datatype>& mpi_type_cache)
{
    // Parameters
    int sizes[3];
    int sub_sizes[3];
    int starts[3];

    // Resolving parameters
    resolveRecvParameters(sizes, sub_sizes, starts, grid_num_cells, direction);

    // Creating types
    defineSubArrayType(out_type, sizes, sub_sizes, starts, mpi_type_cache);
}

}