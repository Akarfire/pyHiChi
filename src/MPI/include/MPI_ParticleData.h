#pragma once

#include <mpi.h>
#include <vector>

#include "FP.h"
#include "Dimension.h"

namespace mpi
{
// Data Structures
// Uniform - particles are stored in a single Particle Array (std::vector)
// Grid - particles are stored in an array of arrays

// Describes data, that is exchanged between ranks
// ALL TYPES ONE MESSAGE
/*
struct ParticleData_Uniform_SoA
{
    // Header
    // [0] -> count
    // [1, 2, ..., pfc::sizeParticleTypes] -> type counts
    int header[pfc::sizeParticleTypes + 1] = {}

    // Data
    // positions_x, positions_y, positions_z
    // momentums_x, momentums_y, momentums_z
    // weights, gammas
    std::vector<double> payload;

    // Accessors

    // Count of all stored particles
    int getCount() { return header[0]; }
    void setCount(int count, bool update_payload_size = true) 
    { 
        header[0] = count; 
        payload.resize(count * 8); // 8 particle data components
    }

    // C-style array of particle counts for each particle type (size = pfc::sizeParticleTypes)
    int* getTypeCounts() { return &header[1]; }
};
*/

// Describes data, that is exchanged between ranks
template<pfc::Dimension dimension>
struct ParticleData_Uniform_SoA
{
    // Header
    // [0] -> count
    // [1] -> type
    int header[2] = {};

    // Data
    // positions_x, positions_y, positions_z
    // momentums_x, momentums_y, momentums_z
    // weights, gammas
    std::vector<pfc::FP> payload;

    // Number of double values per each particle
    constexpr int perParticleDataSize = 
          static_cast<int>(dimension) // Position components
        + 3 // Momentum always has 3 components
        + 2; // Weight and Gamma

    // Accessors

    // Count of all stored particles
    int getCount() const { return header[0]; }
    void setCount(int count, bool update_payload_size = true) 
    { 
        header[0] = count; 
        payload.reserve(count * perParticleDataSize);
    }

    // Type of particles stored in the data
    int getType() const { return header[1]; }
    void setType(int type) { header[1] = type; } 

    // MPI Communication data

    static const MPI_Datatype headerMpiDataType = MPI_INT;
    static const int headerMpiSize = 2;

    static const MPI_Datatype payloadMpiDataType = MPI_DOUBLE;
    // Only valid after header is written
    int getPayloadMpiSize() const { return header[0] * perParticleDataSize; }
};


// Defines
typedef typename ParticleData_Uniform_SoA<Three> ParticleData3D_Uniform_SoA;
}