#pragma once

#include <mpi.h>
#include <vector>

#include "FP.h"
#include "Vectors.h"
#include "Particle.h"
#include "Ensemble.h"

#include "MPI_Utilities.h"
#include "MPI_ParticleUtils.h"
#include "MPI_ParticleData.h"

namespace mpi
{

// 1) Extracts out-of-bound particles from the particle array and prepares them to be sent to another rank
// 2) Inserts particles received from other ranks into the particle array
// Base class with common helper functions
class ParticleStorageInterface
{
protected:
    
    // Helper functions

    // Returns the cube corner id corresponding to the neighboring rank, that a particle 
    // with the specified position shall be sent to
    static int getResponsibleNeighborID(FP x, FP y, FP z, 
                                        const pfc::FP3& lower_bound, 
                                        const pfc::FP3& upper_bound)
    {
        // Formula without branching (to make it faster, since it runs per particle)

        // If particle is below the lower bound (0 or 1)
        int lower_x = x < lower_bound.x;
        int lower_y = y < lower_bound.y;
        int lower_z = z < lower_bound.z;

        // If particle is above the upper bound (0 or 1)
        int higher_x = x > upper_bound.x;
        int higher_y = y > upper_bound.y;
        int higher_z = z > upper_bound.z;

        // higher_c and lower_c are NEVER equal to 1 at the same time

        // Constructing offsets
        int ix = higher_x - lower_x;
        int iy = higher_y - lower_y;
        int iz = higher_z - lower_z;

        return offsetToCubeCornerID(ix, iy, iz);
    }

// Implement in derived classes
/*
public:
    static void collectParticles(   ParticleArrayType* particleArray, 
                                    ParticleDataType outParticleData[27],
                                    const pfc::FP3& lower_bound, const pfc::FP3& upper_bound)

    static void insertParticles(    pfc::ParticleArray3d* particleArray, 
                                    const ParticleData3D_Uniform_SoA inParticleData[27],
                                    const pfc::FP3& lower_bound, const pfc::FP3& upper_bound)
*/
};


// 1) Extracts out-of-bound particles from the particle array and prepares them to be sent to another rank
// 2) Inserts particles received from other ranks into the particle array
// Version for a 3D uniform particle array with SoA data type
class ParticleStorageInterface3D_Uniform_SoA : public ParticleStorageInterface
{
public:
    static void collectParticles(   pfc::ParticleArray3d* particleArray, 
                                    ParticleData3D_Uniform_SoA outParticleData[27],
                                    const pfc::FP3& lower_bound, const pfc::FP3& upper_bound)
    {
        // Send indices for every direction
        // (0, 0, 0) offset is mapped to 13, so it will always be emtpy
        std::vector<int> send_indices[27]; // TO DO : Add a heuristic reservation mechanism
        std::vector<int> deletion_indices;

        // Loop over the ensemble and select the out-of-bound particles
        auto& pos_x = particleArray->getPositions(0);
        auto& pos_y = particleArray->getPositions(1);
        auto& pos_z = particleArray->getPositions(2);

        for (int i = 0; i < particleArray->size(); i++)
        {
            int neighbor_id = getResponsibleNeighborID(pos_x[i], pos_y[i], pos_z[i], lower_bound, upper_bound);
            if (neighbor_id != NO_DIRECTION)
            {
                send_indices[neighbor_id].push_back(i);
                deletion_indices.push_back(i);
            }
        }

        auto& momentums_x = particleArray->getMomentums(0);
        auto& momentums_y = particleArray->getMomentums(1);
        auto& momentums_z = particleArray->getMomentums(2);

        // Packaging loop
        for (int dir = 0; dir < 27; dir++)
        {
            if (dir == NO_DIRECTION) continue; // Skipping the (0, 0, 0) offset (it is not used)

            // Filling out the header
            outParticleData[dir].setCount(send_indices[dir].size(), true);
            outParticleData[dir].setType(particleArray->getType());

            // Package positions
            std::vector<int>& indices = send_indices[dir];

            // TEMPORARY, WILL BE REPLACED WITH COPYING THE END OF ARRAYS
            // (Once the swapping functionality is implemented)

            // Positions
            for (int i = 0; i < indices.size(); i++)
            {
                int index = indices[i];
                outParticleData[dir].payload.push_back(pos_x[index]);
            }

            for (int i = 0; i < indices.size(); i++)
            {
                int index = indices[i];
                outParticleData[dir].payload.push_back(pos_y[index]);
            }

            for (int i = 0; i < indices.size(); i++)
            {
                int index = indices[i];
                outParticleData[dir].payload.push_back(pos_z[index]);
            }

            // Momentums
            for (int i = 0; i < indices.size(); i++)
            {
                int index = indices[i];
                outParticleData[dir].payload.push_back(momentums_x[index]);
            }

            for (int i = 0; i < indices.size(); i++)
            {
                int index = indices[i];
                outParticleData[dir].payload.push_back(momentums_y[index]);
            }

            for (int i = 0; i < indices.size(); i++)
            {
                int index = indices[i];
                outParticleData[dir].payload.push_back(momentums_z[index]);
            }
            
            // Weights
            for (int i = 0; i < indices.size(); i++)
            {
                int index = indices[i];
                outParticleData[dir].payload.push_back(particleArray->getWeights()[index]);
            }

            // Gammas
            for (int i = 0; i < indices.size(); i++)
            {
                int index = indices[i];
                outParticleData[dir].payload.push_back(particleArray->getGammas()[index]);
            }
        }

        // Deletion loop (ALSO TEMPORARY)
        for (int i = static_cast<int>(deletion_indices.size()) - 1; i >= 0; i--)
        {
            particleArray->deleteParticle(deletion_indices[i]);
        }
    }

    static void insertParticles(    pfc::ParticleArray3d* particleArray, 
                                    const ParticleData3D_Uniform_SoA inParticleData[27],
                                    const pfc::FP3& lower_bound, const pfc::FP3& upper_bound)
    {
        // Type check
        if (particleArray->getType() != inParticleData[0].getType())
            throw std::runtime_error("MPI PARTICLE STORAGE INTEFACE : Attempted to insert particles with incorrect type!");

        // Cache initial write index before resizing the buffers
        int writeIndex = particleArray->size();

        // Reserving space for particles
        int totalCount = 0;
        for (int dir = 0; dir < 27; dir++)
        {
            if (dir == NO_DIRECTION) continue; // Skipping the (0, 0, 0) offset (it is not used)
            totalCount += inParticleData[dir].getCount();
        }

        int reserveSize = totalCount + particleArray->size();

        particleArray->getPositions(0).resize(reserveSize);
        particleArray->getPositions(1).resize(reserveSize);
        particleArray->getPositions(2).resize(reserveSize);

        particleArray->getMomentums(0).resize(reserveSize);
        particleArray->getMomentums(1).resize(reserveSize);
        particleArray->getMomentums(2).resize(reserveSize);

        particleArray->getWeights().resize(reserveSize);
        particleArray->getGammas().resize(reserveSize);

        // Unpack loop
        for (int dir = 0; dir < 27; dir++)
        {
            if (dir == NO_DIRECTION) continue;

            const ParticleData3D_Uniform_SoA &src = inParticleData[dir];
            int count = src.getCount();

            if (count == 0) continue;

            // Positions
            std::memcpy(particleArray->getPositions(0).data() + writeIndex,
                        src.payload.data() + 0 * count, count * sizeof(pfc::FP));
            std::memcpy(particleArray->getPositions(1).data() + writeIndex,
                        src.payload.data() + 1 * count, count * sizeof(pfc::FP));
            std::memcpy(particleArray->getPositions(2).data() + writeIndex,
                        src.payload.data() + 2 * count, count * sizeof(pfc::FP));

            // Momentums
            std::memcpy(particleArray->getMomentums(0).data() + writeIndex,
                        src.payload.data() + 3 * count, count * sizeof(pfc::FP));
            std::memcpy(particleArray->getMomentums(1).data() + writeIndex,
                        src.payload.data() + 4 * count, count * sizeof(pfc::FP));
            std::memcpy(particleArray->getMomentums(2).data() + writeIndex,
                        src.payload.data() + 5 * count, count * sizeof(pfc::FP));

            // Weight
            std::memcpy(particleArray->getWeights().data() + writeIndex,
                        src.payload.data() + 6 * count, count * sizeof(pfc::FP));

            // Gamma
            std::memcpy(particleArray->getGammas().data() + writeIndex,
                        src.payload.data() + 7 * count, count * sizeof(pfc::FP));

            writeIndex += count;
        }
    }
};

}