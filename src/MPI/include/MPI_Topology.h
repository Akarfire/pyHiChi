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

public:
    // Constructing an abstracted MPI topology based on the number of nodes (ranks) and the desired topology type
    // `sections` - Number of sections for each axis
    Topology(const pfc::Int3& sections, LoopType loop_type, int node_count);

private:

    // Number of sections for each axis
    pfc::Int3 sections_;
    // Topology looping type
    LoopType loopType;
    // Number of nodes in this topology
    int size;

    // Converts a 3D index into a 1D index
    // Returns MPI_INVALID_RANK if at least one of coordinates is MPI_INVALID_RANK
    int convertIndexToRank(int x, int y, int z) const;
    // Converts a 3D index into a 1D index
    // Returns MPI_INVALID_RANK if at least one of coordinates is MPI_INVALID_RANK
    int convertIndexToRank(pfc::Int3 index) const { return convertIndexToRank(index.x, index.y, index.z); }

    // Converts a rank index into a 3D index
    pfc::Int3 convertRankToIndex(int rank) const;

    static bool doesLoopOverX(const LoopType& loop_type);
    static bool doesLoopOverY(const LoopType& loop_type);
    static bool doesLoopOverZ(const LoopType& loop_type);

public:

    // Returns sections
    const pfc::Int3& getSections() const { return sections_; }

    // Returns the type of looping used for this topology
    LoopType getLoopType() const { return loopType; }

    // Returns the number of nodes in this topology
    int getSize() const { return size; }

    // Returns rank of the node's neighbor in the specified direction
    // MPI_INVALID_RANK if no neighbor exists
    int getNeighbor(int node_rank, Direction direction) const;

    // Returns rank of the node's neighbor with the specified offset
    // MPI_INVALID_RANK if no neighbor exists
    int getNeighbor(int node_rank, pfc::Int3 offset) const;
};

}