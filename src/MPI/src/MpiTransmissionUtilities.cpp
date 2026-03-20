#include "MpiTransmissionUtilities.h"

// Defines an mpi data type for transmitting a line of elements (alligned with coordinate axes)
void MpiTransmissionUtilities::defineLineType(MPI_Datatype &out_type, int count, int stride)
{
    MPI_Type_vector(count, 1, stride, MPI_DOUBLE, &out_type);
    MPI_Type_commit(&out_type);
}

// Defines an mpi data type for transmitting a plane of elements (alligned with coordinate axes)
void MpiTransmissionUtilities::definePlaneType(MPI_Datatype &out_type, const MPI_Datatype &line_type, int count, int stride)
{
    // Calculating stride in bytes
    MPI_Aint byte_stride;
    MPI_Type_extent(MPI_DOUBLE, &byte_stride);
    byte_stride *= stride;
    
    // Creating a vector type with a byte stride
    MPI_Type_create_hvector(count, 1, byte_stride, line_type, &out_type);
    MPI_Type_commit(&out_type);
}

// Resolving "send" operation parameters based on direction
void MpiTransmissionUtilities::resolveSendParameters(int &line_count, int &line_stride, int &plane_count, int &plane_stride, int &offset, 
        const pfc::Int3& num_cells, MpiTransmissionDirection direction)
{
    switch (direction)
    {
    case MpiTransmissionDirection::positiveX:
        // YZ plane
        line_count = num_cells.y;
        line_stride = num_cells.z;
        plane_count = num_cells.z;
        plane_stride = 1;
        // Far end of X axis (second from last)
        offset = num_cells.x * num_cells.y * num_cells.z - 2 * line_count * plane_count;
        break;

    case MpiTransmissionDirection::negativeX:
        // YZ plane
        line_count = num_cells.y;
        line_stride = num_cells.z;
        plane_count = num_cells.z;
        plane_stride = 1;
        // Second YZ plane
        offset = line_count * plane_count;
        break;

    case MpiTransmissionDirection::positiveY:
        // XZ plane
        line_count = num_cells.x;
        line_stride = num_cells.y * num_cells.z;
        plane_count = num_cells.z;
        plane_stride = 1;
        // Far end of Y axis (second from last)
        offset = num_cells.z * (num_cells.y - 2);
        break;

    case MpiTransmissionDirection::negativeY:
        // XZ plane
        line_count = num_cells.x;
        line_stride = num_cells.y * num_cells.z;
        plane_count = num_cells.z;
        plane_stride = 1;
        // Second XZ plane
        offset = num_cells.z;
        break;

    case MpiTransmissionDirection::positiveZ:
        // XY plane
        line_count = num_cells.x;
        line_stride = num_cells.y * num_cells.z;
        plane_count = num_cells.y;
        plane_stride = num_cells.z;
        // Far end of Z axis (second from last)
        offset = num_cells.z - 2;
        break;

    case MpiTransmissionDirection::negativeZ:
        // XY plane
        line_count = num_cells.x;
        line_stride = num_cells.y * num_cells.z;
        plane_count = num_cells.y;
        plane_stride = num_cells.z;
        // Second XY plane
        offset = 1;
        break;

    default: break;
    }
}

// Resolving "recv" operation parameters based on direction
void MpiTransmissionUtilities::resolveRecvParameters(int &line_count, int &line_stride, int &plane_count, int &plane_stride, int &offset, 
        const pfc::Int3& num_cells, MpiTransmissionDirection direction)
{
    // IMPORTANT: Direction is specified relative to the SENDER
    switch (direction)
    {
    case MpiTransmissionDirection::positiveX:
        // YZ plane
        line_count = num_cells.y;
        line_stride = num_cells.z;
        plane_count = num_cells.z;
        plane_stride = 1;
        // First YZ plane
        offset = 0;
        break;

    case MpiTransmissionDirection::negativeX:
        // YZ plane
        line_count = num_cells.y;
        line_stride = num_cells.z;
        plane_count = num_cells.z;
        plane_stride = 1;
        // Last YZ plane
        offset = num_cells.x * num_cells.y * num_cells.z - line_count * plane_count;
        break;

    case MpiTransmissionDirection::positiveY:
        // XZ plane
        line_count = num_cells.x;
        line_stride = num_cells.y * num_cells.z;
        plane_count = num_cells.z;
        plane_stride = 1;
        // First XZ plane
        offset = 0;
        break;

    case MpiTransmissionDirection::negativeY:
        // XZ plane
        line_count = num_cells.x;
        line_stride = num_cells.y * num_cells.z;
        plane_count = num_cells.z;
        plane_stride = 1;
        // Last XZ plane
        offset = num_cells.z * (num_cells.y - 1);
        break;

    case MpiTransmissionDirection::positiveZ:
        // XY plane
        line_count = num_cells.x;
        line_stride = num_cells.y * num_cells.z;
        plane_count = num_cells.y;
        plane_stride = num_cells.z;
        // First XY plane
        offset = 0;
        break;

    case MpiTransmissionDirection::negativeZ:
        // XY plane
        line_count = num_cells.x;
        line_stride = num_cells.y * num_cells.z;
        plane_count = num_cells.y;
        plane_stride = num_cells.z;
        // Last XY plane
        offset = num_cells.z - 1;
        break;

    default: break;
    }
}


// Determines mpi data type and required offset to send field data in the required direction
// All created types are stored into mpi_type_cache, for it to be cleaned up later (using cleanUpMpiTypes)
void MpiTransmissionUtilities::defineMpiTransmission_Send(int& out_offset, MPI_Datatype& out_type, pfc::Int3 grid_num_cells, MpiTransmissionDirection direction,
    std::vector<MPI_Datatype>& mpi_type_cache)
{
    // Parameters
    int line_count;
    int line_stride;
    int plane_count;
    int plane_stride;

    int offset;

    // Resolving parameters
    resolveSendParameters(  line_count, line_stride, 
                            plane_count, plane_stride, 
                            offset, grid_num_cells, direction);

    // Creating types
    MPI_Datatype lineType;
    mpi_type_cache.push_back(lineType);
    defineLineType(mpi_type_cache[mpi_type_cache.size() - 1], line_count, line_stride);

    MPI_Datatype planeType;
    mpi_type_cache.push_back(planeType);
    definePlaneType(mpi_type_cache[mpi_type_cache.size() - 1], lineType, plane_count, plane_stride);

    // Outputting
    out_offset = offset;
    out_type = mpi_type_cache[mpi_type_cache.size() - 1];
}

// Determines mpi data type and required offset to received field data from the required direction (direction relative to the sender)
// All created types are stored into mpi_type_cache, for it to be cleaned up later (using cleanUpMpiTypes)
void MpiTransmissionUtilities::defineMpiTransmission_Recv(int& out_offset, MPI_Datatype& out_type, pfc::Int3 grid_num_cells, MpiTransmissionDirection direction,
    std::vector<MPI_Datatype>& mpi_type_cache)
{
    // Parameters
    int line_count;
    int line_stride;
    int plane_count;
    int plane_stride;

    int offset;

    // Resolving parameters
    resolveRecvParameters(  line_count, line_stride, 
                            plane_count, plane_stride, 
                            offset, grid_num_cells, direction);

    // Creating types
    MPI_Datatype lineType;
    mpi_type_cache.push_back(lineType);
    defineLineType(mpi_type_cache[mpi_type_cache.size() - 1], line_count, line_stride);

    MPI_Datatype planeType;
    mpi_type_cache.push_back(planeType);
    definePlaneType(mpi_type_cache[mpi_type_cache.size() - 1], lineType, plane_count, plane_stride);

    // Outputting
    out_offset = offset;
    out_type = mpi_type_cache[mpi_type_cache.size() - 1];
}

// Calls MPI_Type_free on all types stored in mpi_type_cache
void MpiTransmissionUtilities::cleanUpMpiTypes(std::vector<MPI_Datatype>& mpi_type_cache)
{
    for (auto& type : mpi_type_cache)
        MPI_Type_free(&type);
}