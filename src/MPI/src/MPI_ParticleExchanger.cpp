#include "MPI_ParticleExchanger.h"

#include "ParticleTypes.h"

namespace mpi
{

// Returns an offset corresponding to the neighbouring rank, that a particle 
// with the specified position shall be sent to
// If resulting offset is (0, 0, 0) : out_send = false, otherwise: out_send = true
pfc::Int3 ParticleExchanger::getResponsibleNeighbor(   const pfc::FP3& position, 
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

    return pfc::Int3(x, y, z);
}

// Exchanges out-of-bound particles between MPI ranks
void ParticleExchanger::exchangeParticles(  int mpi_rank, std::unique_ptr<Ensemble3d>& ensemble, 
                                            const pfc::FP3& lower_bound, const pfc::FP3& upper_bound,
                                            std::shared_ptr<class Topology> topology)
{
    // Send indices for every direction
    // (0, 0, 0) offset is mapped to 13, so it will always be emtpy
    std::vector<int> send_indices[27][pfc::sizeParticleTypes];

    // Loop over the ensemble and select the out-of-bound particles
    for (int type = 0; type < pfc::sizeParticleTypes; type++)
        for (int i = 0; i < (*ensemble)[type].size(); i++)
        {
            pfc::ParticleProxy3d particleProxy = (*ensemble)[type][i];
            pfc::Int3 neighbor_offset = getResponsibleNeighbor(particleProxy.getPosition(), lower_bound, upper_bound);
            if (neighbor_offset != MPI_SAME_RANK)
            {
                send_indices[offsetToCubeCornerID(neighbor_offset)][type].push_back(i);
            }
        }

    // Packaging loop
    SendData3D data[27];
    for (int dir = 0; dir < 27; dir++)
    {
        if (dir == 13) continue; // Skipping the (0, 0, 0) offset (it is not used)

        // Filling out the header
        data[dir].count = 0;
        for (int type = 0; type < pfc::sizeParticleTypes; type++)
        {
            int type_count = static_cast<int>(send_indices[dir][type].size());
            data[dir].count += type_count;
            data[dir].type_counts[type] = type_count;
        }

        // Package positions
        for (int type = 0; type < pfc::sizeParticleTypes; type++)
        {
            std::vector<int>& indices = send_indices[dir][type];
            for (int i = 0; i < indices.size(); i++)
            {
                int index = indices[i];
                pfc::ParticleProxy3d particleProxy = (*ensemble)[type][index];

                data[dir].positions[0].push_back(particleProxy.getPosition().x);
                data[dir].positions[1].push_back(particleProxy.getPosition().y);
                data[dir].positions[2].push_back(particleProxy.getPosition().z);
            }
        }

        // ...
    }

    // Deletion loop
    // ...

    // Communication loop
    for (int dir = 0; dir < 27; dir++)
    {
        if (dir == 13) continue; // Skipping the (0, 0, 0) offset (it is not used)
        pfc::Int3 offset = cubeCornerIdToOffset(dir);

        // ... MPI_Send and stuff here ...
    }
}

}