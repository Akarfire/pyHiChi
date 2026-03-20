#pragma once

#include <mpi.h>

#include "Grid.h"

// Wrapper class for mpi utility functions, CAN NOT BE INSTANCED!
// All method are static
class MpiTransmissionUtilities final
{
private:
    // Private constructor to prevent instancing
    MpiTransmissionUtilities() {}

public:
    // Defines the direction, in which the transmission will be performed (From sender to receiver)
    enum class MpiTransmissionDirection
    {
        positiveX,
        negativeX,
        positiveY,
        negativeY,
        positiveZ,
        negativeZ
    };

private:
    // Defines an mpi data type for transmitting a line of elements (alligned with coordinate axes)
    static void defineLineType(MPI_Datatype &out_type, int count, int stride);

    // Defines an mpi data type for transmitting a plane of elements (alligned with coordinate axes)
    static void definePlaneType(MPI_Datatype &out_type, const MPI_Datatype &line_type, int count, int stride);

    // Resolving "send" operation parameters based on direction
    static void resolveSendParameters(int &line_count, int &line_stride, int &plane_count, int &plane_stride, int &offset, 
        const pfc::Int3& num_cells, MpiTransmissionDirection direction);

    // Resolving "recv" operation parameters based on direction
    static void resolveRecvParameters(int &line_count, int &line_stride, int &plane_count, int &plane_stride, int &offset, 
        const pfc::Int3& num_cells, MpiTransmissionDirection direction);

public:

    // Determines mpi data type and required offset to send field data in the required direction
    // All created types are stored into mpi_type_cache, for it to be cleaned up later (using cleanUpMpiTypes)
    static void defineMpiTransmission_Send(int& out_offset, MPI_Datatype& out_type, pfc::Int3 grid_num_cells, MpiTransmissionDirection direction,
        std::vector<MPI_Datatype>& mpi_type_cache);

    // Determines mpi data type and required offset to received field data from the required direction (direction relative to the sender)
    // All created types are stored into mpi_type_cache, for it to be cleaned up later (using cleanUpMpiTypes)
    static void defineMpiTransmission_Recv(int& out_offset, MPI_Datatype& out_type, pfc::Int3 grid_num_cells, MpiTransmissionDirection direction,
        std::vector<MPI_Datatype>& mpi_type_cache);

    // Calls MPI_Type_free on all types stored in mpi_type_cache
    static void cleanUpMpiTypes(std::vector<MPI_Datatype>& mpi_type_cache);
};