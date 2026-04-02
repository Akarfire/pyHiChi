#include "MPI_Topology.h"

namespace mpi
{
    
// Constructing an abstracted MPI topology based on the number of nodes (ranks) and the desired topology type
Topology::Topology(const pfc::Int3& sections, LoopType loop_type, int node_count) : sections_(sections), loopType(loop_type)
{
    // Checking if topology is possible
    if (node_count < sections.x * sections.y * sections.z)
        throw std::runtime_error("MPI TOPOLOGY : Impossible topology for the specified number of nodes");

    size = sections.x * sections.y * sections.z;
}

// Converts a 3D index into a 1D index
// Returns MPI_INVALID_RANK if at least one of coordinates is MPI_INVALID_RANK
int Topology::convertIndexToRank(int x, int y, int z) const
{
    if (x == MPI_INVALID_RANK || y == MPI_INVALID_RANK || z == MPI_INVALID_RANK)
        return MPI_INVALID_RANK;

    return x * sections_.y * sections_.z + y * sections_.z + z;
}

 // Converts a rank index into a 3D index
pfc::Int3 Topology::convertRankToIndex(int rank) const
{
    if (rank == MPI_INVALID_RANK)
        return pfc::Int3(MPI_INVALID_RANK, MPI_INVALID_RANK, MPI_INVALID_RANK);
    
    pfc::Int3 index;

    index.x = rank / (sections_.y * sections_.z);
    index.y = (rank / sections_.z) % sections_.y;
    index.z = rank % sections_.z;
    
    return index;
}

bool Topology::doesLoopOverX(const LoopType& loop_type)
{
    return loop_type == LoopType::loopX 
        || loop_type == LoopType::loopXY 
        || loop_type == LoopType::loopXZ 
        || loop_type == LoopType::loopXYZ;
}
bool Topology::doesLoopOverY(const LoopType& loop_type)
{
    return loop_type == LoopType::loopY
        || loop_type == LoopType::loopXY 
        || loop_type == LoopType::loopYZ 
        || loop_type == LoopType::loopXYZ;
}
bool Topology::doesLoopOverZ(const LoopType& loop_type)
{
    return loop_type == LoopType::loopZ
        || loop_type == LoopType::loopXZ 
        || loop_type == LoopType::loopYZ 
        || loop_type == LoopType::loopXYZ;
}


// Returns rank of the node's neighbor in the specified direction
// MPI_INVALID_RANK if no neighbor exists
int Topology::getNeighbor(int node_rank, Direction direction) const
{
    return getNeighbor(node_rank, mpi::DirectionToOffset[static_cast<int>(direction)]);
}

// Returns rank of the node's neighbor with the specified offset
// MPI_INVALID_RANK if no neighbor exists
int Topology::getNeighbor(int node_rank, pfc::Int3 offset) const
{
    if (node_rank == MPI_INVALID_RANK)
        return MPI_INVALID_RANK;

    pfc::Int3 index = convertRankToIndex(node_rank);
    pfc::Int3 neighbor_index = index + offset;
    
    // Validating and looping neighbor_index
    neighbor_index.x = (neighbor_index.x >= 0) ? (neighbor_index.x) : ( doesLoopOverX(loopType) ? (sections_.x + neighbor_index.x) : MPI_INVALID_RANK );
    neighbor_index.x = (neighbor_index.x < sections_.x) ? (neighbor_index.x) : ( doesLoopOverX(loopType) ? (neighbor_index.x - sections_.x) : MPI_INVALID_RANK );

    neighbor_index.y = (neighbor_index.y >= 0) ? (neighbor_index.y) : ( doesLoopOverX(loopType) ? (sections_.y + neighbor_index.y) : MPI_INVALID_RANK );
    neighbor_index.y = (neighbor_index.y < sections_.y) ? (neighbor_index.y) : ( doesLoopOverX(loopType) ? (neighbor_index.y - sections_.y) : MPI_INVALID_RANK );

    neighbor_index.z = (neighbor_index.z >= 0) ? (neighbor_index.z) : ( doesLoopOverX(loopType) ? (sections_.z + neighbor_index.z) : MPI_INVALID_RANK );
    neighbor_index.z = (neighbor_index.z < sections_.z) ? (neighbor_index.z) : ( doesLoopOverX(loopType) ? (neighbor_index.z - sections_.z) : MPI_INVALID_RANK );

    // Converting neighbor_index to rank
    return convertIndexToRank(neighbor_index);
}

}