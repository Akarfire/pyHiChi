#pragma once

#include <vector>

#include "Vectors.h"
#include "MPI_Utilities.h"

namespace mpi
{

// A class that describes abstracted mpi topology
class Topology final
{
public:

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
    // `sections` - Number of sections for each axis
    Topology(const pfc::Int3& sections, LoopType loop_type, int node_count);

private:

    // Number of sections for each axis
    pfc::Int3 sections_;
    // Topology looping type
    LoopType loopType;

    // NodeData for each node (rank), that exists in this topology
    // Access using `int rank`
    std::vector<NodeData> data;

    // Converts a 3D index into a 1D index (used only in constructor)
    // Returns MPI_INVALID_RANK if at least one of coordinates is MPI_INVALID_RANK
    int convertRankIndex(int x, int y, int z);

    bool doesLoopOverX(const LoopType& loop_type);
    bool doesLoopOverY(const LoopType& loop_type);
    bool doesLoopOverZ(const LoopType& loop_type);

public:

    // Returns sections
    const pfc::Int3& getSections() { return sections_; }

    // Returns the type of looping used for this topology
    LoopType getLoopType() { return loopType; }

    // Returns rank of the node's neighbor in the specified direction
    // MPI_INVALID_RANK if no neighbor exists
    int getNeighbor(int node_rank, Direction direction) const;

    // Returns data about the specified node (rank)
    const NodeData& getNodeData(int node_rank) const { return data[node_rank]; }
};

}