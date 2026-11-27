#pragma once

#include "MPI_ParticleData.h"
#include "MPI_ParticleStorageInterface.h"

namespace mpi
{

// Helper classes for setting up particle transfering 

enum class ParticleTransferType
{
    ArrayUniform3DSoA,
    EnsembleUniform3DSoA,
};

template<ParticleTransferType TType>
struct ParticleTransfer {};

template<>
struct ParticleTransfer<ParticleTransferType::ArrayUniform3DSoA> 
{
    typedef ParticleData3D_Uniform_SoA ParticleData;
    typedef ParticleStorageInterface3D_Uniform_SoA StorageInterface;
};

template<>
struct ParticleTransfer<ParticleTransferType::EnsembleUniform3DSoA> 
{
    typedef ParticleData3D_Uniform_SoA_Ensemble ParticleData;
    typedef ParticleStorageInterface3D_Uniform_SoA_Ensemble StorageInterface;
};

} // namepsace mpi