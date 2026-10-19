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

/* Example (mandatory only)

struct ParticleData
{
    // Header
    int header[...] = {};

    // Data
    std::vector<PayloadDataType> payload;

    // Accessors

    // Count of all stored particles
    // Total count of stored particles of all types
    int getCount() const { return <resolve based on header>; }

    // Per type data access
    int getTypeCount(int type) const { <resolve based on header>; }
    PayloadDataType* getTypePayload(int type) { return <pointer to payload section for this type>; } 

    // MPI Communication data
    static const MPI_Datatype headerMpiDataType = MPI_INT;
    static const int headerMpiSize = 2;

    static const MPI_Datatype payloadMpiDataType = MPI_DOUBLE;
    // Only valid after header is written
    int getPayloadMpiSize() const { return <resolve based on header>; }

    // Custom mpi types support
    static const bool useCustomHeaderMpiDatatype = <true / false>;
    static const bool useCustomPayloadMpiDatatype = <true / false>;
    static const int headerCustomMpiDatatype = <index (same as in createCustomMpiDatatypes) or -1 if not used>;
    static const int payloadCustomMpiDatatype = <index (same as in createCustomMpiDatatypes) or -1 if not used>;
    static void createCustomMpiDatatypes(std::vector<MPI_Datatype>& outCustomTypeBuffer) 
    { 
        // Create types and put them into outCustomTypeBuffer onto desired indices
    }
};

*/

// Describes data, that is exchanged between ranks
// Single type
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
    static inline constexpr int perParticleDataSize()
    { 
        return  static_cast<int>(dimension) // Position components
                + 3 // Momentum always has 3 components
                + 2; // Weight and Gamma
    }

    // Accessors

    // Count of all stored particles
    int getCount() const { return header[0]; }
    void setCount(int count, bool update_payload_size = true) 
    { 
        header[0] = count; 
        payload.reserve(count * perParticleDataSize());
    }
    // Type (for message validation only)
    int getType() const { return header[1]; }
    void setType(int type) { header[1] = type; }

    // Per type data access (any type in this case)
    int getTypeCount(int type) const { return header[0]; }
    const pfc::FP* getTypePayload(int type) const { return payload.data(); } 

    // MPI Communication data
    static const MPI_Datatype headerMpiDataType = MPI_INT;
    static const int headerMpiSize = 2;

    static const MPI_Datatype payloadMpiDataType = MPI_DOUBLE;
    // Only valid after header is written
    int getPayloadMpiSize() const { return header[0] * perParticleDataSize(); }

    // Custom mpi types support
    static const bool useCustomHeaderMpiDatatype = false;
    static const bool useCustomPayloadMpiDatatype = false;
    static const int headerCustomMpiDatatypeIndex = -1;
    static const int payloadCustomMpiDatatypeIndex = -1;
    static void createCustomMpiDatatypes(std::vector<MPI_Datatype>& outCustomTypeBuffer) {}
};


// Defines
typedef typename ParticleData_Uniform_SoA<Three> ParticleData3D_Uniform_SoA;
}