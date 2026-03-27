#pragma once

#include <vector>

#include "MPI_Utilities.h"

namespace mpi
{

// A class that describes abstracted mpi topology
class Topology final
{
public:
    // MPI topology type
    enum class Type
    {
        lineX, lineY, lineZ,
        gridXY, gridYZ, gridXZ,
        gridXYZ
    };

    // MPI topology looping type
    enum class LoopType
    {
        None,
        loopX, loopY, loopZ,
        loopXY, loopYZ, loopXZ,
        loopXYZ
    };

    // MPI-related parameters associated with a specifid node (rank)
    struct NodeData
    {
        // List of neighbouring ranks
        // Access using `mpi::Direction direction`
        std::vector<int> neighbors;
    };

public:
    // Constructing an abstracted MPI topology based on the number of nodes (ranks) and the desired topology type
    Topology(Type topology_type, LoopType loop_type, int node_count);

private:

    // Topology type
    Type type;
    // Topology looping type
    LoopType loopType;

    // NodeData for each node (rank), that exists in this topology
    // Access using `int rank`
    std::vector<NodeData> data;

    bool doesLoopOverX(const LoopType& loop_type);
    bool doesLoopOverY(const LoopType& loop_type);
    bool doesLoopOverZ(const LoopType& loop_type);

    // Initializes line-type (lineX, lineY or lineZ) type topologies (used in constructor)
    void initLineTopology(std::vector<NodeData>& topology_data, Type type, LoopType loop_type, int node_count);

    // Initializes 2D-Grid-type (gridXY, gridYZ or gridXZ) type topologies (used in constructor)
    void init2DGridTopology(std::vector<NodeData>& topology_data, Type type, LoopType loop_type, int node_count);

    // Initializes 3D-Grid-type (gridXYZ) type topologies (used in constructor)
    void init3DGridTopology(std::vector<NodeData>& topology_data, LoopType loop_type, int node_count);

public:

    // Returns the type of the topology
    Type getType() { return type; }

    // Returns the type of looping used for this topology
    LoopType getLoopType() { return loopType; }

    // Returns rank of the node's neighbor in the specified direction
    // MPI_INVALID_RANK if no neighbor exists
    int getNeighbor(int node_rank, Direction direction) const;

    // Returns data about the specified node (rank)
    const NodeData& getNodeData(int node_rank) const { return data[node_rank]; }
};

}