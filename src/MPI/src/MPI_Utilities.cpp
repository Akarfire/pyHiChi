#include "MPI_Utilities.h"

namespace mpi
{

// Calls MPI_Type_free on all types stored in mpi_type_cache
void cleanUpMpiTypes(std::vector<MPI_Datatype>& mpi_type_cache)
{
    for (auto& type : mpi_type_cache)
        MPI_Type_free(&type);
}


// Defines an mpi data type for transmitting a line of elements (alligned with coordinate axes)
void MPI_FieldUtils::defineLineType(MPI_Datatype &out_type, int count, int stride)
{
    MPI_Type_vector(count, 1, stride, MPI_DOUBLE, &out_type);
    MPI_Type_commit(&out_type);
}

// Defines an mpi data type for transmitting a plane of elements (alligned with coordinate axes)
void MPI_FieldUtils::definePlaneType(MPI_Datatype &out_type, const MPI_Datatype &line_type, int count, int stride)
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
void MPI_FieldUtils::resolveSendParameters(int &line_count, int &line_stride, int &plane_count, int &plane_stride, int &offset, 
        const pfc::Int3& num_cells, Direction direction)
{
    switch (direction)
    {
    case Direction::positiveX:
        // YZ plane
        line_count = num_cells.y;
        line_stride = num_cells.z;
        plane_count = num_cells.z;
        plane_stride = 1;
        // Far end of X axis (second from last)
        offset = num_cells.x * num_cells.y * num_cells.z - 2 * line_count * plane_count;
        break;

    case Direction::negativeX:
        // YZ plane
        line_count = num_cells.y;
        line_stride = num_cells.z;
        plane_count = num_cells.z;
        plane_stride = 1;
        // Second YZ plane
        offset = line_count * plane_count;
        break;

    case Direction::positiveY:
        // XZ plane
        line_count = num_cells.x;
        line_stride = num_cells.y * num_cells.z;
        plane_count = num_cells.z;
        plane_stride = 1;
        // Far end of Y axis (second from last)
        offset = num_cells.z * (num_cells.y - 2);
        break;

    case Direction::negativeY:
        // XZ plane
        line_count = num_cells.x;
        line_stride = num_cells.y * num_cells.z;
        plane_count = num_cells.z;
        plane_stride = 1;
        // Second XZ plane
        offset = num_cells.z;
        break;

    case Direction::positiveZ:
        // XY plane
        line_count = num_cells.x;
        line_stride = num_cells.y * num_cells.z;
        plane_count = num_cells.y;
        plane_stride = num_cells.z;
        // Far end of Z axis (second from last)
        offset = num_cells.z - 2;
        break;

    case Direction::negativeZ:
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
void MPI_FieldUtils::resolveRecvParameters(int &line_count, int &line_stride, int &plane_count, int &plane_stride, int &offset, 
        const pfc::Int3& num_cells, Direction direction)
{
    // IMPORTANT: Direction is specified relative to the SENDER
    switch (direction)
    {
    case Direction::positiveX:
        // YZ plane
        line_count = num_cells.y;
        line_stride = num_cells.z;
        plane_count = num_cells.z;
        plane_stride = 1;
        // First YZ plane
        offset = 0;
        break;

    case Direction::negativeX:
        // YZ plane
        line_count = num_cells.y;
        line_stride = num_cells.z;
        plane_count = num_cells.z;
        plane_stride = 1;
        // Last YZ plane
        offset = num_cells.x * num_cells.y * num_cells.z - line_count * plane_count;
        break;

    case Direction::positiveY:
        // XZ plane
        line_count = num_cells.x;
        line_stride = num_cells.y * num_cells.z;
        plane_count = num_cells.z;
        plane_stride = 1;
        // First XZ plane
        offset = 0;
        break;

    case Direction::negativeY:
        // XZ plane
        line_count = num_cells.x;
        line_stride = num_cells.y * num_cells.z;
        plane_count = num_cells.z;
        plane_stride = 1;
        // Last XZ plane
        offset = num_cells.z * (num_cells.y - 1);
        break;

    case Direction::positiveZ:
        // XY plane
        line_count = num_cells.x;
        line_stride = num_cells.y * num_cells.z;
        plane_count = num_cells.y;
        plane_stride = num_cells.z;
        // First XY plane
        offset = 0;
        break;

    case Direction::negativeZ:
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
void MPI_FieldUtils::defineMpiTransmission_Send(int& out_offset, MPI_Datatype& out_type, pfc::Int3 grid_num_cells, Direction direction,
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
    defineLineType(lineType, line_count, line_stride);
    mpi_type_cache.push_back(lineType);

    MPI_Datatype planeType;
    definePlaneType(planeType, lineType, plane_count, plane_stride);
    mpi_type_cache.push_back(planeType);

    // Outputting
    out_offset = offset;
    out_type = planeType;
}

// Determines mpi data type and required offset to received field data from the required direction (direction relative to the sender)
// All created types are stored into mpi_type_cache, for it to be cleaned up later (using cleanUpMpiTypes)
void MPI_FieldUtils::defineMpiTransmission_Recv(int& out_offset, MPI_Datatype& out_type, pfc::Int3 grid_num_cells, Direction direction,
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
    defineLineType(lineType, line_count, line_stride);
    mpi_type_cache.push_back(lineType);

    MPI_Datatype planeType;
    definePlaneType(planeType, lineType, plane_count, plane_stride);
    mpi_type_cache.push_back(planeType);

    // Outputting
    out_offset = offset;
    out_type = planeType;
}
}