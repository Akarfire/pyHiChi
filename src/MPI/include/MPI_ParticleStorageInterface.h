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
#include "MPI_Topology.h"

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
    static pfc::Int3 getResponsibleNeighbor(FP x, FP y, FP z, 
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

        return pfc::Int3(ix, iy, iz);
    }

// Implement in derived classes
/*
public:
    static void collectParticles(   ParticleArrayType* particleArray, 
                                    ParticleDataType outParticleData[27],
                                    const pfc::FP3& lower_bound, const pfc::FP3& upper_bound,
                                    int mpi_rank, std::shared_ptr<class Topology> topology,
                                    int expected_num_particles[27])

    static void insertParticles(    pfc::ParticleArray3d* particleArray, 
                                    const ParticleData3D_Uniform_SoA inParticleData[27],
                                    const pfc::FP3& lower_bound, const pfc::FP3& upper_bound,
                                    int mpi_rank, std::shared_ptr<class Topology> topology)
*/
};


// 1) Extracts out-of-bound particles from the particle array and prepares them to be sent to another rank
// 2) Inserts particles received from other ranks into the particle array
// Version for a 3D uniform particle array with SoA data type
class ParticleStorageInterface3D_Uniform_SoA : public ParticleStorageInterface
{
    // Swaps particle at particleArray[index] with the particle before the tail
    // Returns the new index of the particle
    static int swapParticle(pfc::ParticleArray3d* particleArray, int index, int tail)
    {        
        int size = particleArray->size();
        int new_index = size - tail - 1;

        if (new_index == index) return index;

        for (int d = 0; d < particleArray->positionDimension; d++)
        {
            auto& positions = particleArray->getPositions(d);
            std::swap(positions[index], positions[new_index]);
        }
        for (int d = 0; d < particleArray->momentumDimension; d++)
        {
            auto& momentums = particleArray->getMomentums(d);
            std::swap(momentums[index], momentums[new_index]);
        }

        auto& weights = particleArray->getWeights();
        std::swap(weights[index], weights[new_index]);
        auto& gammas = particleArray->getGammas();
        std::swap(gammas[index], gammas[new_index]);

        return new_index;
    }

public:
    static void collectParticles(   pfc::ParticleArray3d* particleArray, 
                                    ParticleData3D_Uniform_SoA outParticleData[27],
                                    const pfc::FP3& lower_bound, const pfc::FP3& upper_bound,
                                    int mpi_rank, std::shared_ptr<Topology> topology,
                                    int expected_num_particles[27])
    {
        // Send indices for every direction
        // (0, 0, 0) offset is mapped to 13, so it will always be emtpy
        std::vector<int> send_indices[27];
        bool validDirection[27] = {};

        for (int dir = 0; dir < 27; dir++)
        {
            if (dir == NO_DIRECTION) continue;
            send_indices[dir].reserve(expected_num_particles[dir]);
            validDirection[dir] = topology->getNeighbor(mpi_rank, cubeCornerIdToOffset(dir)) != MPI_INVALID_RANK;
        }

        int size = particleArray->size();

        // Loop over the ensemble and select the out-of-bound particles
        auto& pos_x = particleArray->getPositions(0);
        auto& pos_y = particleArray->getPositions(1);
        auto& pos_z = particleArray->getPositions(2);

        int num_transfered = 0;
        for (int i = 0; i < size - num_transfered; i++)
        {
            pfc::Int3 neighbor_offset = getResponsibleNeighbor(pos_x[i], pos_y[i], pos_z[i], lower_bound, upper_bound);
            int neighbor_id = offsetToCubeCornerID(neighbor_offset);
            if (neighbor_id != NO_DIRECTION && validDirection[neighbor_id])
            {
                int new_index = swapParticle(particleArray, i, num_transfered);
                send_indices[neighbor_id].push_back(new_index);
                num_transfered++;

                i--;
            }
        }

        auto& momentums_x = particleArray->getMomentums(0);
        auto& momentums_y = particleArray->getMomentums(1);
        auto& momentums_z = particleArray->getMomentums(2);

        auto& weights = particleArray->getWeights();
        auto& gammas = particleArray->getGammas();

        // Packaging loop
        for (int dir = 0; dir < 27; dir++)
        {
            if (dir == NO_DIRECTION) continue; // Skipping the (0, 0, 0) offset (it is not used)

            // Filling out the header
            outParticleData[dir].setCount(send_indices[dir].size());
            outParticleData[dir].setType(particleArray->getType());

            // Resizing payload
            outParticleData[dir].resizePayload(send_indices[dir].size());

            // Package positions
            std::vector<FP>& payload = outParticleData[dir].payload;
            std::vector<int>& indices = send_indices[dir];
            int n = static_cast<int>(indices.size());
    
            // Positions
            for (int i = 0; i < n; i++) payload[n * 0 + i] = pos_x[indices[i]];
            for (int i = 0; i < n; i++) payload[n * 1 + i] = pos_y[indices[i]];
            for (int i = 0; i < n; i++) payload[n * 2 + i] = pos_z[indices[i]];

            // Momentums
            for (int i = 0; i < n; i++) payload[n * 3 + i] = momentums_x[indices[i]];
            for (int i = 0; i < n; i++) payload[n * 4 + i] = momentums_y[indices[i]];
            for (int i = 0; i < n; i++) payload[n * 5 + i] = momentums_z[indices[i]];
            
            // Weights
            for (int i = 0; i < n; i++) payload[n * 6 + i] = weights[indices[i]];

            // Gammas
            for (int i = 0; i < n; i++) payload[n * 7 + i] = gammas[indices[i]];
        }

        // Deletion
        particleArray->getPositions(0).resize(size - num_transfered);
        particleArray->getPositions(1).resize(size - num_transfered);
        particleArray->getPositions(2).resize(size - num_transfered);

        particleArray->getMomentums(0).resize(size - num_transfered);
        particleArray->getMomentums(1).resize(size - num_transfered);
        particleArray->getMomentums(2).resize(size - num_transfered);

        particleArray->getWeights().resize(size - num_transfered);
        particleArray->getGammas().resize(size - num_transfered);

        // Handling periodical conditions - Pre send
        for (int dir = 0; dir < 27; dir++)
        {
            if (dir == NO_DIRECTION) continue;

            pfc::Int3 offset = cubeCornerIdToOffset(dir);
            bool periodical = topology->isNeighborOverLoop(mpi_rank, offset);

            if (!periodical) continue;

            pfc::FP* data = outParticleData[dir].payload.data();
            int particle_count = outParticleData[dir].getCount();

            pfc::FP* dpos_x = data + 0 * particle_count;
            pfc::FP* dpos_y = data + 1 * particle_count;
            pfc::FP* dpos_z = data + 2 * particle_count;

            // X
            if (offset.x < 0)
                for (int p = 0; p < particle_count; p++)
                    dpos_x[p] = lower_bound.x - dpos_x[p];

            else if (offset.x > 0)
                for (int p = 0; p < particle_count; p++)
                    dpos_x[p] = dpos_x[p] - upper_bound.x;

            // Y
            if (offset.y < 0)
                for (int p = 0; p < particle_count; p++)
                    dpos_y[p] = lower_bound.y - dpos_y[p];

            else if (offset.y > 0)
                for (int p = 0; p < particle_count; p++)
                    dpos_y[p] = dpos_y[p] - upper_bound.y;

            // Z
            if (offset.z < 0)
                for (int p = 0; p < particle_count; p++)
                    dpos_z[p] = lower_bound.z - dpos_z[p];

            else if (offset.z > 0)
                for (int p = 0; p < particle_count; p++)
                    dpos_z[p] = dpos_z[p] - upper_bound.z;
        }
    }

    static void insertParticles(    pfc::ParticleArray3d* particleArray, 
                                    ParticleData3D_Uniform_SoA inParticleData[27],
                                    const pfc::FP3& lower_bound, const pfc::FP3& upper_bound,
                                    int mpi_rank, std::shared_ptr<Topology> topology)
    {
        // Type check (since ParticleData3D_Uniform_SoA is single type)
        if (particleArray->getType() != inParticleData[0].getType())
            throw std::runtime_error("MPI PARTICLE STORAGE INTEFACE : Attempted to insert particles with incorrect type!");


        // Handling periodical conditions - Post send
        for (int dir = 0; dir < 27; dir++)
        {
            if (dir == NO_DIRECTION) continue;

            pfc::Int3 offset = cubeCornerIdToOffset(dir);
            bool periodical = topology->isNeighborOverLoop(mpi_rank, -1 * offset); // Offset is inversed

            if (!periodical) continue;

            pfc::FP* data = inParticleData[dir].payload.data();
            int particle_count = inParticleData[dir].getCount();

            pfc::FP* dpos_x = data + 0 * particle_count;
            pfc::FP* dpos_y = data + 1 * particle_count;
            pfc::FP* dpos_z = data + 2 * particle_count;

            // X
            if (offset.x < 0)
                for (int p = 0; p < particle_count; p++)
                    dpos_x[p] = upper_bound.x - dpos_x[p];

            else if (offset.x > 0)
                for (int p = 0; p < particle_count; p++)
                    dpos_x[p] = dpos_x[p] + lower_bound.x;

            // Y
            if (offset.y < 0)
                for (int p = 0; p < particle_count; p++)
                    dpos_y[p] = upper_bound.y - dpos_y[p];

            else if (offset.y > 0)
                for (int p = 0; p < particle_count; p++)
                    dpos_y[p] = dpos_y[p] + lower_bound.y;

            // Z
            if (offset.z < 0)
                for (int p = 0; p < particle_count; p++)
                    dpos_z[p] = upper_bound.z - dpos_z[p];

            else if (offset.z > 0)
                for (int p = 0; p < particle_count; p++)
                    dpos_z[p] = dpos_z[p] + lower_bound.z;
        }


        // Cache initial write index before resizing the buffers
        int writeIndex = particleArray->size();

        // Reserving space for particles
        int totalCount = 0;
        for (int dir = 0; dir < 27; dir++)
        {
            if (dir == NO_DIRECTION) continue; // Skipping the (0, 0, 0) offset (it is not used)
            totalCount += inParticleData[dir].getTypeCount(particleArray->getType());
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
            int count = src.getTypeCount(particleArray->getType());
            const pfc::FP *data = src.getTypePayload(particleArray->getType());

            if (count == 0) continue;

            // Positions
            std::memcpy(particleArray->getPositions(0).data() + writeIndex,
                        data + 0 * count, count * sizeof(pfc::FP));
            std::memcpy(particleArray->getPositions(1).data() + writeIndex,
                        data + 1 * count, count * sizeof(pfc::FP));
            std::memcpy(particleArray->getPositions(2).data() + writeIndex,
                        data + 2 * count, count * sizeof(pfc::FP));

            // Momentums
            std::memcpy(particleArray->getMomentums(0).data() + writeIndex,
                        data + 3 * count, count * sizeof(pfc::FP));
            std::memcpy(particleArray->getMomentums(1).data() + writeIndex,
                        data + 4 * count, count * sizeof(pfc::FP));
            std::memcpy(particleArray->getMomentums(2).data() + writeIndex,
                        data + 5 * count, count * sizeof(pfc::FP));

            // Weight
            std::memcpy(particleArray->getWeights().data() + writeIndex,
                        data + 6 * count, count * sizeof(pfc::FP));

            // Gamma
            std::memcpy(particleArray->getGammas().data() + writeIndex,
                        data + 7 * count, count * sizeof(pfc::FP));

            writeIndex += count;
        }
    }
};

}