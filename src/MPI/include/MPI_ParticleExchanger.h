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
// All methods are static
class ParticleExchanger final
{
private:
    // Private constructor to prevent instancing
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

    // Returns the cube corner id corresponding to the neighboring rank, that a particle 
    // with the specified position shall be sent to
    static int getResponsibleNeighborID(FP x, FP y, FP z, 
                                        const pfc::FP3& lower_bound, 
                                        const pfc::FP3& upper_bound);

    // Returns an offset corresponding to the neighboring rank, that a particle 
    // with the specified position shall be sent to
    static pfc::Int3 getResponsibleNeighbor(const pfc::FP3& position, 
                                            const pfc::FP3& lower_bound, 
                                            const pfc::FP3& upper_bound);


    // Prepares SendData for every direction, by extracting out-of-bound particles from the ensemble.
    // Particles are deleted from the ensemble
    static void prepareSendData(SendData3D out_data[27], std::unique_ptr<Ensemble3d>& ensemble, 
                                const pfc::FP3& lower_bound, const pfc::FP3& upper_bound);

public:

    // Exchanges out-of-bound particles between MPI ranks
    static void exchangeParticles(  int mpi_rank, std::unique_ptr<Ensemble3d>& ensemble, 
                                    const pfc::FP3& lower_bound, const pfc::FP3& upper_bound,
                                    std::shared_ptr<class Topology> topology);

};

}