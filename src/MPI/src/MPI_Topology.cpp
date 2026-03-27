#include "MPI_Topology.h"

namespace mpi
{
    
// Constructing an abstracted MPI topology based on the number of nodes (ranks) and the desired topology type
Topology::Topology(Type topology_type, LoopType loop_type, int node_count) : type(topology_type), loopType(loop_type)
{
    // Line topology
    if (topology_type == Type::lineX || topology_type == Type::lineY || topology_type == Type::lineZ)
        initLineTopology(data, topology_type, loop_type, node_count);

    // 2D Grid topology
    else if (topology_type == Type::gridXY || topology_type == Type::gridYZ || topology_type == Type::gridXZ)
        init2DGridTopology(data, topology_type, loop_type, node_count);
    
    // 3D Grid topology
    init3DGridTopology(data, loop_type, node_count);
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


// Initializes line-type (lineX, lineY or lineZ) type topologies (used in constructor)
void Topology::initLineTopology(std::vector<NodeData>& topology_data, Type type, LoopType loop_type, int node_count)
{
    topology_data.clear();

    for (int rank = 0; rank < node_count; rank++)
    {
        NodeData nodeData = NodeData();
        
        for (int i = 0; i < 6; i++)
            nodeData.neighbors.push_back(MPI_INVALID_RANK);

        if (node_count > 1)
        {           
            switch(type)
            {
            case Type::lineX:
                nodeData.neighbors[static_cast<int>(Direction::positiveX)] = 
                    (rank + 1 < node_count) ? (rank + 1) : (doesLoopOverX(loop_type) ? 0 : MPI_INVALID_RANK);

                nodeData.neighbors[static_cast<int>(Direction::negativeX)] = 
                    (rank - 1 >= 0) ? (rank - 1) : (doesLoopOverX(loop_type) ? (node_count - 1) : MPI_INVALID_RANK);

                if (doesLoopOverY(loop_type))
                {
                    nodeData.neighbors[static_cast<int>(Direction::positiveY)] = rank;
                    nodeData.neighbors[static_cast<int>(Direction::negativeY)] = rank;
                }
                if (doesLoopOverZ(loop_type))
                {
                    nodeData.neighbors[static_cast<int>(Direction::positiveZ)] = rank;
                    nodeData.neighbors[static_cast<int>(Direction::negativeZ)] = rank;
                }
                break;

            case Type::lineY:
                nodeData.neighbors[static_cast<int>(Direction::positiveY)] = 
                    (rank + 1 < node_count) ? (rank + 1) : (doesLoopOverY(loop_type) ? 0 : MPI_INVALID_RANK);

                nodeData.neighbors[static_cast<int>(Direction::negativeY)] = 
                    (rank - 1 >= 0) ? (rank - 1) : (doesLoopOverY(loop_type) ? (node_count - 1) : MPI_INVALID_RANK);

                if (doesLoopOverX(loop_type))
                {
                    nodeData.neighbors[static_cast<int>(Direction::positiveX)] = rank;
                    nodeData.neighbors[static_cast<int>(Direction::negativeX)] = rank;
                }
                if (doesLoopOverZ(loop_type))
                {
                    nodeData.neighbors[static_cast<int>(Direction::positiveZ)] = rank;
                    nodeData.neighbors[static_cast<int>(Direction::negativeZ)] = rank;
                }
                break;

            case Type::lineZ:
                nodeData.neighbors[static_cast<int>(Direction::positiveZ)] = 
                    (rank + 1 < node_count) ? (rank + 1) : (doesLoopOverZ(loop_type) ? 0 : MPI_INVALID_RANK);

                nodeData.neighbors[static_cast<int>(Direction::negativeZ)] = 
                    (rank - 1 >= 0) ? (rank - 1) : (doesLoopOverZ(loop_type) ? (node_count - 1) : MPI_INVALID_RANK);

                if (doesLoopOverX(loop_type))
                {
                    nodeData.neighbors[static_cast<int>(Direction::positiveX)] = rank;
                    nodeData.neighbors[static_cast<int>(Direction::negativeX)] = rank;
                }
                if (doesLoopOverY(loop_type))
                {
                    nodeData.neighbors[static_cast<int>(Direction::positiveY)] = rank;
                    nodeData.neighbors[static_cast<int>(Direction::negativeY)] = rank;
                }
                break;

            default: break;
            }
        }

        topology_data.push_back(nodeData);
    }
}

// Initializes 2D-Grid-type (gridXY, gridYZ or gridXZ) type topologies (used in constructor)
void Topology::init2DGridTopology(std::vector<NodeData>& topology_data, Type type, LoopType loop_type, int node_count)
{
    // TO DO : Implement this
}

// Initializes 3D-Grid-type (gridXYZ) type topologies (used in constructor)
void Topology::init3DGridTopology(std::vector<NodeData>& topology_data, LoopType loop_type, int node_count)
{
    // TO DO : Implement this
}

// Returns rank of the node's neighbor in the specified direction
// MPI_INVALID_RANK if no neighbor exists
int Topology::getNeighbor(int node_rank, Direction direction) const
{
    const NodeData& data = getNodeData(node_rank);
    return data.neighbors[static_cast<int>(direction)];
}

}