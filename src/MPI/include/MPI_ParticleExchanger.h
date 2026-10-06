#pragma once

#include <mpi.h>
#include <vector>

#include "FP.h"
#include "Vectors.h"
#include "Particle.h"
#include "Ensemble.h"

#include "MPI_Utilities.h"

namespace mpi
{

// Class responsible for exchanging particles between MPI ranks
class ParticleExchanger final
{
public:
    // Constructs necessary mpi types (if needed, not sure yet)
    ParticleExchanger() {}

private:
    // ... If needed mpi types will be here ...

private:
    // Helper functions

    // Returns an offset corresponding to the neighbouring rank, that a particle 
    // with the specified position shall be sent to
    // If resulting offset is (0, 0, 0) : out_send = false, otherwise: out_send = true
    pfc::Int3 getResponsibleNeighbour(  bool& out_send,
                                        const pfc::FP3& position, 
                                        const pfc::FP3& lower_bound, 
                                        const pfc::FP3& upper_bound);

public:

    // Exchanges out-of-bound particles between MPI ranks
    void exchangeParticles( int mpi_rank, std::unique_ptr<Ensemble3d>& ensemble, 
                            const pfc::FP3& lower_bound, const pfc::FP3& upper_bound);

};

}