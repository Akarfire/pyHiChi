#include "MPI_Utilities.h"

namespace mpi
{

// Returns the inverse of the specified direction
Direction invertDirection(const Direction& direction)
{
    switch (direction)
    {
    case Direction::positiveX: return Direction::negativeX;
    case Direction::negativeX: return Direction::positiveX;
    case Direction::positiveY: return Direction::negativeY;
    case Direction::negativeY: return Direction::positiveY;
    case Direction::positiveZ: return Direction::negativeZ;
    case Direction::negativeZ: return Direction::positiveZ;

    default: return Direction::positiveX;
    }
}

// TOPOLOGY

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


// FIELD UTILS

// Defines an mpi data type for transmitting sub arrays of the grid
void FieldUtils::defineSubArrayType(MPI_Datatype &out_type, int sizes[3], int sub_sizes[3], int starts[3])
{
    MPI_Type_create_subarray(3, sizes, sub_sizes, starts, MPI_ORDER_C, MPI_DOUBLE, &out_type);
    MPI_Type_commit(&out_type);
}

// Resolving "send" operation parameters based on direction
void FieldUtils::resolveSendParameters(int out_sizes[3], int out_sub_sizes[3], int out_starts[3], 
        const pfc::Int3& num_cells, Direction direction)
{
    out_sizes[0] = num_cells.x;
    out_sizes[1] = num_cells.y;
    out_sizes[2] = num_cells.z;

    // General outputs (to avoid code repettition)
    out_sub_sizes[0] = num_cells.x;
    out_sub_sizes[1] = num_cells.y;
    out_sub_sizes[2] = num_cells.z;
    out_starts[0] = 0;
    out_starts[1] = 0;
    out_starts[2] = 0;

    // Specifying outputs based on direction
    switch (direction)
    {
    case Direction::positiveX:
        // Second to last YZ plane
        out_sub_sizes[0] = 1;
        out_starts[0] = num_cells.x - 2;
        break;

    case Direction::negativeX:
        // Second YZ plane
        out_sub_sizes[0] = 1;
        out_starts[0] = 1;
        break;

    case Direction::positiveY:
        // Second to last XZ plane
        out_sub_sizes[1] = 1;
        out_starts[1] = num_cells.y - 2;
        break;

    case Direction::negativeY:
        // Second XZ plane
        out_sub_sizes[1] = 1;
        out_starts[1] = 1;
        break;

    case Direction::positiveZ:
        // Second to last XY plane
        out_sub_sizes[2] = 1;
        out_starts[2] = num_cells.z - 2;
        break;

    case Direction::negativeZ:
        // Second XY plane
        out_sub_sizes[2] = 1;
        out_starts[2] = 1;
        break;

    default: break;
    }
}

// Resolving "recv" operation parameters based on direction
void FieldUtils::resolveRecvParameters(int out_sizes[3], int out_sub_sizes[3], int out_starts[3], 
        const pfc::Int3& num_cells, Direction direction)
{
    out_sizes[0] = num_cells.x;
    out_sizes[1] = num_cells.y;
    out_sizes[2] = num_cells.z;

    // General outputs (to avoid code repettition)
    out_sub_sizes[0] = num_cells.x;
    out_sub_sizes[1] = num_cells.y;
    out_sub_sizes[2] = num_cells.z;
    out_starts[0] = 0;
    out_starts[1] = 0;
    out_starts[2] = 0;

    // Specifying outputs based on direction
    // IMPORTANT: Direction is specified relative to the SENDER
    switch (direction)
    {
    case Direction::positiveX:
        // First YZ plane
        out_sub_sizes[0] = 1;
        out_starts[0] = 0;
        break;

    case Direction::negativeX:
        // Last YZ plane
        out_sub_sizes[0] = 1;
        out_starts[0] = num_cells.x - 1;
        break;

    case Direction::positiveY:
        // First XZ plane
        out_sub_sizes[1] = 1;
        out_starts[1] = 0;
        break;

    case Direction::negativeY:
        // Last XZ plane
        out_sub_sizes[1] = 1;
        out_starts[1] = num_cells.y - 1;
        break;

    case Direction::positiveZ:
        // First XY plane
        out_sub_sizes[2] = 1;
        out_starts[2] = 0;
        break;

    case Direction::negativeZ:
        // Last XY plane
        out_sub_sizes[2] = 1;
        out_starts[2] = num_cells.z - 1;
        break;

    default: break;
    }
}


// Determines mpi data type and required offset to send field data in the required direction
void FieldUtils::defineTransmission_Send(MPI_Datatype& out_type, pfc::Int3 grid_num_cells, Direction direction)
{
    // Parameters
    int sizes[3];
    int sub_sizes[3];
    int starts[3];

    // Resolving parameters
    resolveSendParameters(sizes, sub_sizes, starts, grid_num_cells, direction);

    // Creating type
    defineSubArrayType(out_type, sizes, sub_sizes, starts);
}

// Determines mpi data type and required offset to received field data from the required direction (direction relative to the sender)
void FieldUtils::defineTransmission_Recv(MPI_Datatype& out_type, pfc::Int3 grid_num_cells, Direction direction)
{
    // Parameters
    int sizes[3];
    int sub_sizes[3];
    int starts[3];

    // Resolving parameters
    resolveRecvParameters(sizes, sub_sizes, starts, grid_num_cells, direction);

    // Creating types
    defineSubArrayType(out_type, sizes, sub_sizes, starts);
}


// FIELD EXCHANGER

// Constructor that will create types for sending and recieving 
FieldExchanger::FieldExchanger(const pfc::Int3& grid_num_cells)
{
    // Creating types
    for (int d = 0; d < 6; d++)
    {
        Direction direction = static_cast<Direction>(d);

        // Send
        MPI_Datatype send_type;
        FieldUtils::defineTransmission_Send(send_type, grid_num_cells, direction);
        sendTypes.push_back(send_type);

        // Recv
        MPI_Datatype recv_type;
        FieldUtils::defineTransmission_Recv(recv_type, grid_num_cells, direction);
        recvTypes.push_back(recv_type);
    }
}

// Frees created mpi types
FieldExchanger::~FieldExchanger()
{
    int mpiInit;
    MPI_Initialized(&mpiInit);

    int mpiFinalized;
    MPI_Finalized(&mpiFinalized);

    if (mpiInit && !mpiFinalized)
        for (int i = 0; i < 6; i++)
        {
            MPI_Type_free(&sendTypes[i]);
            MPI_Type_free(&recvTypes[i]);
        }
}

// Performs echange sequence iteration for node "rank" in the specified topology
// MUST BE CALLED BY EVERY PROCESS IN THE TOPOLOGY AT THE SAME TIME!
void FieldExchanger::PerformExchangeSequence(pfc::FP* data, const Topology& topology, int rank, MPI_Comm communicator)
{
    for (int dir = 0; dir < 6; dir++)
    {
        Direction direction = static_cast<Direction>(dir);

        // Neighbors
        int send_neighbor = topology.getNeighbor(rank, direction);
        int recv_neighbor = topology.getNeighbor(rank, invertDirection(direction));

        // Sync nodes
        MPI_Barrier(communicator);

        MPI_Request request;
        MPI_Status status;

        // Sending
        if (send_neighbor != MPI_INVALID_RANK)
            MPI_Isend(data, 1, sendTypes[dir], send_neighbor, 0, communicator, &request);
        
        // Receiving
        if (recv_neighbor != MPI_INVALID_RANK)
            MPI_Recv(data, 1, recvTypes[dir], recv_neighbor, 0, communicator, &status);

        // Waiting for send operation to complete
        if (send_neighbor != MPI_INVALID_RANK)
            MPI_Wait(&request, &status);
    }

    // Sync nodes
    MPI_Barrier(communicator);
}

}