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

// Returns the inverse of the specified direction
Direction invertDirection(const Direction& direction);

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


// Wrapper class for mpi field utility functions, CAN NOT BE INSTANCED!
// All methods are static
class FieldUtils final
{
private:
    // Private constructor to prevent instancing
    FieldUtils() {}

private:
    // Defines an mpi data type for transmitting sub arrays of the grid
    static void defineSubArrayType(MPI_Datatype &out_type, int sizes[3], int sub_sizes[3], int starts[3]);

    // Resolving "send" operation parameters based on direction
    static void resolveSendParameters(int out_sizes[3], int out_sub_sizes[3], int out_starts[3], 
        const pfc::Int3& num_cells, Direction direction);

    // Resolving "recv" operation parameters based on direction
    static void resolveRecvParameters(int out_sizes[3], int out_sub_sizes[3], int out_starts[3], 
        const pfc::Int3& num_cells, Direction direction);

public:

    // Determines mpi data type required to send field data in the required direction
    static void defineTransmission_Send(MPI_Datatype& out_type, pfc::Int3 grid_num_cells, Direction direction);

     // Determines mpi data type required to received field data from the required direction (direction relative to the sender)
    static void defineTransmission_Recv(MPI_Datatype& out_type, pfc::Int3 grid_num_cells, Direction direction);
};


// Performs field exchange sequence
class FieldExchanger final
{
private:
    // MPI types used for sending field data
    // Access using `mpi::Direction direction`
    std::vector<MPI_Datatype> sendTypes;

    // MPI types used for receiving field data
    // Access using `mpi::Direction direction`
    std::vector<MPI_Datatype> recvTypes;

public:

    // Constructor that will create types for sending and recieving 
    FieldExchanger(const pfc::Int3& grid_num_cells);
    // Frees created mpi types
    ~FieldExchanger();

    // Performs echange sequence iteration for node "rank" in the specified topology
    // MUST BE CALLED BY EVERY PROCESS IN THE TOPOLOGY
    void PerformExchangeSequence(pfc::FP* data, const Topology& topology, int rank, MPI_Comm communicator = MPI_COMM_WORLD);
};

}