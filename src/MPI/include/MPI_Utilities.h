#pragma once

#include <mpi.h>

#include "Grid.h"

#define MPI_INVALID_RANK -1

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

        // Node color, used to arrange transmissions in the correct order
        // Neighboring nodes MUST have differing colors 
        unsigned int color;
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
    int getNeighbor(int node_rank, Direction direction);

    // Returns data about the specified node (rank)
    const NodeData& getNodeData(int node_rank) { return data[node_rank]; }
};


// Calls MPI_Type_free on all types stored in mpi_type_cache
void cleanUpMpiTypes(std::vector<MPI_Datatype>& mpi_type_cache);


// Wrapper class for mpi field utility functions, CAN NOT BE INSTANCED!
// All methods are static
class FieldUtils final
{
private:
    // Private constructor to prevent instancing
    FieldUtils() {}

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