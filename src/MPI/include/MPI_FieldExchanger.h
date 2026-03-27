#pragma once

#include <mpi.h>
#include <vector>

#include "FP.h"
#include "Vectors.h"


namespace mpi
{
    
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
    // MUST BE CALLED BY EVERY PROCESS IN THE TOPOLOGY AT THE SAME TIME!
    void PerformExchangeSequence(pfc::FP* data, const class Topology& topology, int rank, MPI_Comm communicator = MPI_COMM_WORLD);
};
}