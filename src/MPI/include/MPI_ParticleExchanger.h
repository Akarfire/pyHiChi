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
    // Describes data, that is exchanged between ranks
    template<int positionDimension, int momentumDimension>
    struct SendData
    {
        // Header
        int count;
        int type_counts[pfc::sizeParticleTypes];

        // Data
        std::vector<double> positions[positionDimension];
        std::vector<double> momentums[momentumDimension];
        std::vector<double> weights;
        std::vector<double> gammas;
    };

    using SendData3D = SendData<3, 3>;

private:
    // Helper functions

    // Returns an offset corresponding to the neighboring rank, that a particle 
    // with the specified position shall be sent to
    pfc::Int3 getResponsibleNeighbor(  const pfc::FP3& position, 
                                        const pfc::FP3& lower_bound, 
                                        const pfc::FP3& upper_bound);

public:

    // Exchanges out-of-bound particles between MPI ranks
    void exchangeParticles( int mpi_rank, std::unique_ptr<Ensemble3d>& ensemble, 
                            const pfc::FP3& lower_bound, const pfc::FP3& upper_bound,
                            std::shared_ptr<class Topology> topology);

};

}