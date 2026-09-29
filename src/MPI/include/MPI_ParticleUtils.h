#pragma once

#include <mpi.h>
#include <vector>

#include "FP.h"
#include "Vectors.h"
#include "Particle.h"

#include "MPI_Utilities.h"

namespace mpi
{
// Wrapper class for mpi particle (pfc::Particle<dimension>) utility functions, CAN NOT BE INSTANCED!
// All methods are static
template<pfc::Dimension dimension>
class ParticleUtils final
{
private:
    // Private constructor to prevent instancing
    ParticleUtils() {}

private:
    // Layout mirror to avoid editing Particle.h
    struct ParticleLayoutMirror
    {
        typename pfc::Particle<dimension>::PositionType position;
        typename pfc::Particle<dimension>::MomentumType p;
        typename pfc::Particle<dimension>::WeightType weight;
        typename pfc::Particle<dimension>::GammaType gamma;
        typename pfc::Particle<dimension>::TypeIndexType typeIndex;
    };
    // Number of fields in pfc::Particle<dimension> (same as in ParticleLayoutMirror):
    static constexpr int ParticleMemberCount = 5;

    // Validates whether pfc::Particle
    static void compileTimeValidation()
    {
        static_assert(std::is_trivially_copyable<pfc::Particle<dimension>>::value, "MPI PARTICLE UTILITIES: pfc::Particle<dimension> type must be trivially copyable!");
        static_assert(sizeof(ParticleLayoutMirror) == sizeof(pfc::Particle<dimension>), "MPI PARTICLE UTILITIES: Layout mirror does not match pfc::Particle<dimension> size!");
        static_assert(std::is_same<FP, double>::value, "MPI PARTICLE UTILITIES: FP is assumed to be double while it is not!");
        static_assert(sizeof(pfc::ParticleTypes) == sizeof(int), "MPI PARTICLE UTILITIES: ParticleTypes must be int-sized!");
    }

public:

    // Creates an MPI type for pfc::Particle<dimension>
    static MPI_Datatype defineParticleType()
    {
        compileTimeValidation();

        // Defining member types

        // Position Type
        MPI_Datatype position_type;
        MPI_Type_contiguous(dimension, MPI_DOUBLE, &position_type);
        MPI_Type_commit(&position_type);

        // Momentum Type
        MPI_Datatype momentum_type;
        MPI_Type_contiguous(3, MPI_DOUBLE, &momentum_type);
        MPI_Type_commit(&momentum_type);

        // The rest of the types are trivial:
            // WeightType type is MPI_DOUBLE
            // GammaType type is MPI_DOUBLE
            // TypeIndexType is MPI_Int

        // Calculating member offsets
        ParticleLayoutMirror mirror;

        MPI_Aint base;
        MPI_Get_address(&mirror, &base);
        
        MPI_Aint offsets[ParticleMemberCount];
        MPI_Get_address(&mirror.position,   &offsets[0]);
        MPI_Get_address(&mirror.p,          &offsets[1]);
        MPI_Get_address(&mirror.weight,     &offsets[2]);
        MPI_Get_address(&mirror.gamma,      &offsets[3]);
        MPI_Get_address(&mirror.typeIndex,  &offsets[4]);
        for (int i = 0; i < ParticleMemberCount; i++)
            offsets[i] -= base;

        // Creating particle MPI type
        int block_lengths[ParticleMemberCount] = {1, 1, 1, 1, 1};
        MPI_Datatype types[ParticleMemberCount] = {
            position_type,
            momentum_type,
            MPI_DOUBLE,
            MPI_DOUBLE,
            MPI_INT
        };

        MPI_Datatype particle_type;
        MPI_Type_create_struct(ParticleMemberCount, block_lengths, offsets, types, &particle_type);
        MPI_Type_commit(&particle_type);
        
        // Resizing to account for tail padding
        MPI_Aint lower_bound, extent;
        MPI_Type_get_extent(particle_type, &lower_bound, &extent);

        MPI_Datatype resized_particle_type;
        MPI_Type_create_resized(particle_type, lower_bound, sizeof(pfc::Particle<dimension>), &resized_particle_type);
        MPI_Type_commit(&resized_particle_type);

        // Freeing member types (not used after struct creation and resizing)
        MPI_Type_free(&particle_type);
        MPI_Type_free(&position_type);
        MPI_Type_free(&momentum_type);

        // Returning type
        return resized_particle_type;
    }
};
}