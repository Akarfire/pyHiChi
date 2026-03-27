#include "MPI_Topology.h"

namespace mpi
{
    
// Constructing an abstracted MPI topology based on the number of nodes (ranks) and the desired topology type
Topology::Topology(const pfc::Int3& sections, LoopType loop_type, int node_count) : sections_(sections), loopType(loop_type)
{
    // Checking if topology is possible
    if (node_count < sections.x * sections.y * sections.z)
        throw std::runtime_error("MPI TOPOLOGY : Impossible topology for the specified number of nodes");

    // Constructing topology
    for (int x = 0; x < sections.x; x++)
        for (int y = 0; y < sections.y; y++)
            for (int z = 0; z < sections.z; z++)
            {
                NodeData nodeData = NodeData();

                int p_x = convertRankIndex((x + 1 < sections.x ? (x + 1) : (doesLoopOverX(loop_type) ? 0 : MPI_INVALID_RANK)), y, z);
                int n_x = convertRankIndex((x - 1 >= 0) ? (x - 1) : (doesLoopOverX(loop_type) ? (sections.x - 1) : MPI_INVALID_RANK), y, z);
                
                int p_y = convertRankIndex(x, (y + 1 < sections.y ? (y + 1) : (doesLoopOverY(loop_type) ? 0 : MPI_INVALID_RANK)), z);
                int n_y = convertRankIndex(x, (y - 1 >= 0) ? (y - 1) : (doesLoopOverY(loop_type) ? (sections.y - 1) : MPI_INVALID_RANK), z);

                int p_z = convertRankIndex(x, y, (z + 1 < sections.z ? (z + 1) : (doesLoopOverZ(loop_type) ? 0 : MPI_INVALID_RANK)));
                int n_z = convertRankIndex(x, y, (z - 1 >= 0) ? (z - 1) : (doesLoopOverZ(loop_type) ? (sections.z - 1) : MPI_INVALID_RANK));

                nodeData.neighbors.resize(6);
                nodeData.neighbors[static_cast<int>(Direction::positiveX)] = p_x;
                nodeData.neighbors[static_cast<int>(Direction::negativeX)] = n_x;
                nodeData.neighbors[static_cast<int>(Direction::positiveY)] = p_y;
                nodeData.neighbors[static_cast<int>(Direction::negativeY)] = n_y;
                nodeData.neighbors[static_cast<int>(Direction::positiveZ)] = p_z;
                nodeData.neighbors[static_cast<int>(Direction::negativeZ)] = n_z;

                data.push_back(nodeData);
            }
}

// Converts a 3D index into a 1D index (used only in constructor)
// Returns MPI_INVALID_RANK if at least one of coordinates is MPI_INVALID_RANK
int Topology::convertRankIndex(int x, int y, int z)
{
    if (x == MPI_INVALID_RANK || y == MPI_INVALID_RANK || z == MPI_INVALID_RANK)
        return MPI_INVALID_RANK;

    return x * sections_.y * sections_.z + y * sections_.z + z;
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
    const NodeData& data = getNodeData(node_rank);
    return data.neighbors[static_cast<int>(direction)];
}

}