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

class ParticleExchanger
{
protected:
    static const int PARTICLE_HEADER_TAG_BASE = 100;
    static const int PARTICLE_PAYLOAD_TAG_BASE = 200;
};

// Class responsible for exchanging particles between MPI ranks
// All methods are static
class ParticleExchanger_Uniform_SoA
{
private:
    // Private constructor to prevent instancing
    ParticleExchanger_Uniform_SoA() {}

private:
    // Helper functions


    // Prepares SendData for every direction, by extracting out-of-bound particles from the ensemble.
    // Particles are deleted from the ensemble
    static void prepareSendData(ParticleData3D out_data[27], std::unique_ptr<Ensemble3d>& ensemble, 
                                const pfc::FP3& lower_bound, const pfc::FP3& upper_bound);

public:

    // Exchanges out-of-bound particles between MPI ranks
    static void exchangeParticles(  int mpi_rank, std::unique_ptr<Ensemble3d>& ensemble, 
                                    const pfc::FP3& lower_bound, const pfc::FP3& upper_bound,
                                    std::shared_ptr<class Topology> topology, MPI_Comm communicator);

};

}