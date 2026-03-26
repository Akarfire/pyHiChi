#pragma once

#include <mpi.h>

#include "Grid.h"

namespace mpi
{
// Defines the direction, in which the transmission will be performed (From sender to receiver)
enum class Direction
{
    positiveX,
    negativeX,
    positiveY,
    negativeY,
    positiveZ,
    negativeZ
};


// Calls MPI_Type_free on all types stored in mpi_type_cache
void cleanUpMpiTypes(std::vector<MPI_Datatype>& mpi_type_cache);

// Wrapper class for mpi utility functions, CAN NOT BE INSTANCED!
// All method are static
class MPI_FieldUtils final
{
private:
    // Private constructor to prevent instancing
    MPI_FieldUtils() {}

private:
    // Defines an mpi data type for transmitting sub arrays of the grid
    static void defineSubArrayType(MPI_Datatype &out_type, int sizes[3], int sub_sizes[3], int starts[3], std::vector<MPI_Datatype>& mpi_type_cache);

    // Resolving "send" operation parameters based on direction
    static void resolveSendParameters(int out_sizes[3], int out_sub_sizes[3], int out_starts[3], 
        const pfc::Int3& num_cells, Direction direction);

    // Resolving "recv" operation parameters based on direction
    static void resolveRecvParameters(int out_sizes[3], int out_sub_sizes[3], int out_starts[3], 
        const pfc::Int3& num_cells, Direction direction);

public:

    // Determines mpi data type required to send field data in the required direction
    // All created types are stored into mpi_type_cache, for it to be cleaned up later (using cleanUpMpiTypes)
    static void defineTransmission_Send(MPI_Datatype& out_type, pfc::Int3 grid_num_cells, Direction direction,
        std::vector<MPI_Datatype>& mpi_type_cache);

     // Determines mpi data type required to received field data from the required direction (direction relative to the sender)
    // All created types are stored into mpi_type_cache, for it to be cleaned up later (using cleanUpMpiTypes)
    static void defineTransmission_Recv(MPI_Datatype& out_type, pfc::Int3 grid_num_cells, Direction direction,
        std::vector<MPI_Datatype>& mpi_type_cache);
};

}