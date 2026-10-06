#include "MPI_ParticleExchanger.h"

#include "ParticleTypes.h"

namespace mpi
{

// Returns an offset corresponding to the neighbouring rank, that a particle 
// with the specified position shall be sent to
// If resulting offset is (0, 0, 0) : out_send = false, otherwise: out_send = true
pfc::Int3 ParticleExchanger::getResponsibleNeighbour(   bool& out_send,
                                                        const pfc::FP3& position, 
                                                        const pfc::FP3& lower_bound, 
                                                        const pfc::FP3& upper_bound)
{    
    // Formula without branching (to make it faster, since it runs per particle)

    // If particle is below the lower bound (0 or 1)
    int lower_x = position.x < lower_bound.x;
    int lower_y = position.y < lower_bound.y;
    int lower_z = position.z < lower_bound.z;

    // If particle is above the upper bound (0 or 1)
    int higher_x = position.x > upper_bound.x;
    int higher_y = position.y > upper_bound.y;
    int higher_z = position.z > upper_bound.z;

    // higher_c and lower_c are NEVER equal to 1 at the same time

    // Constructing offsets
    int x = higher_x - lower_x;
    int y = higher_y - lower_y;
    int z = higher_z - lower_z;

    out_send = (x | y | z) != 0;
    return pfc::Int3(x, y, z);
}

// Exchanges out-of-bound particles between MPI ranks
void ParticleExchanger::exchangeParticles(  int mpi_rank, std::unique_ptr<Ensemble3d>& ensemble, 
                                            const pfc::FP3& lower_bound, const pfc::FP3& upper_bound)
{
    // WIP
    // Loop over the ensemble and pack out-of-bound particles
    for (int type = 0; type < pfc::sizeParticleTypes; type++)
        for (int i = 0; i < (*ensemble)[type].size(); i++)
        {
            bool send;
            pfc::ParticleProxy3d particleProxy = (*ensemble)[type][i];
            getResponsibleNeighbour(send, particleProxy.getPosition(), lower_bound, upper_bound);
            if (send)
            {
                // Pack the particle
            }
        }
}

}