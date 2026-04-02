#include "MPI_FieldUtils.h"

namespace mpi
{

// Defines an mpi data type for transmitting sub arrays of the grid
void FieldUtils::defineSubArrayType(MPI_Datatype &out_type, int sizes[3], int sub_sizes[3], int starts[3])
{
    MPI_Type_create_subarray(3, sizes, sub_sizes, starts, MPI_ORDER_C, MPI_DOUBLE, &out_type);
    MPI_Type_commit(&out_type);
}

// Resolving "send" operation parameters based on direction
void FieldUtils::resolveSendParameters(int out_sizes[3], int out_sub_sizes[3], int out_starts[3], 
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
    case Direction::pX:
        // Second to last YZ plane
        out_sub_sizes[0] = 1;
        out_starts[0] = num_cells.x - 2;
        break;

    case Direction::nX:
        // Second YZ plane
        out_sub_sizes[0] = 1;
        out_starts[0] = 1;
        break;

    case Direction::pY:
        // Second to last XZ plane
        out_sub_sizes[1] = 1;
        out_starts[1] = num_cells.y - 2;
        break;

    case Direction::nY:
        // Second XZ plane
        out_sub_sizes[1] = 1;
        out_starts[1] = 1;
        break;

    case Direction::pZ:
        // Second to last XY plane
        out_sub_sizes[2] = 1;
        out_starts[2] = num_cells.z - 2;
        break;

    case Direction::nZ:
        // Second XY plane
        out_sub_sizes[2] = 1;
        out_starts[2] = 1;
        break;

    default: break;
    }
}

// Resolving "recv" operation parameters based on direction
void FieldUtils::resolveRecvParameters(int out_sizes[3], int out_sub_sizes[3], int out_starts[3], 
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
    case Direction::pX:
        // First YZ plane
        out_sub_sizes[0] = 1;
        out_starts[0] = 0;
        break;

    case Direction::nX:
        // Last YZ plane
        out_sub_sizes[0] = 1;
        out_starts[0] = num_cells.x - 1;
        break;

    case Direction::pY:
        // First XZ plane
        out_sub_sizes[1] = 1;
        out_starts[1] = 0;
        break;

    case Direction::nY:
        // Last XZ plane
        out_sub_sizes[1] = 1;
        out_starts[1] = num_cells.y - 1;
        break;

    case Direction::pZ:
        // First XY plane
        out_sub_sizes[2] = 1;
        out_starts[2] = 0;
        break;

    case Direction::nZ:
        // Last XY plane
        out_sub_sizes[2] = 1;
        out_starts[2] = num_cells.z - 1;
        break;

    default: break;
    }
}


// Determines mpi data type and required offset to send field data in the required direction
void FieldUtils::defineTransmission_Send(MPI_Datatype& out_type, pfc::Int3 grid_num_cells, Direction direction)
{
    // Parameters
    int sizes[3];
    int sub_sizes[3];
    int starts[3];

    // Resolving parameters
    resolveSendParameters(sizes, sub_sizes, starts, grid_num_cells, direction);

    // Creating type
    defineSubArrayType(out_type, sizes, sub_sizes, starts);
}

// Determines mpi data type and required offset to received field data from the required direction (direction relative to the sender)
void FieldUtils::defineTransmission_Recv(MPI_Datatype& out_type, pfc::Int3 grid_num_cells, Direction direction)
{
    // Parameters
    int sizes[3];
    int sub_sizes[3];
    int starts[3];

    // Resolving parameters
    resolveRecvParameters(sizes, sub_sizes, starts, grid_num_cells, direction);

    // Creating types
    defineSubArrayType(out_type, sizes, sub_sizes, starts);
}

}