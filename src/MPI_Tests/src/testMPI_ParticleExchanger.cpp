#include "gtest/gtest.h"
#include <mpi.h>

#include <random>
#include <iostream>

#include "MPI_Utilities.h"
#include "MPI_ParticleStorageInterface.h"
#include "MPI_ParticleExchanger.h"

using namespace mpi;

using ParticleWriteType = pfc::Particle<pfc::Three>;
using ParticleReadType = pfc::ParticleProxy<pfc::Three>;

void setupParticles()
{
    pfc::ParticleInfo::typesVector = { { pfc::constants::electronMass, pfc::constants::electronCharge },
        { pfc::constants::electronMass, -pfc::constants::electronCharge },
        { pfc::constants::protonMass, -pfc::constants::electronCharge } };
        pfc::ParticleInfo::types = &pfc::ParticleInfo::typesVector[0];
        pfc::ParticleInfo::numTypes = sizeParticleTypes;
}

pfc::FP3 limitPosition(const pfc::FP3& pos, const pfc::FP3& lower_lim, const pfc::FP3& upper_lim)
{
    pfc::FP3 new_pos = pos;
    if (pos.x < lower_lim.x)
        new_pos.x = lower_lim.x;
    if (pos.x > upper_lim.x)
        new_pos.x = upper_lim.x;

    if (pos.y < lower_lim.y)
        new_pos.y = lower_lim.y;
    if (pos.y > upper_lim.y)
        new_pos.y = upper_lim.y;

    if (pos.z < lower_lim.z)
        new_pos.z = lower_lim.z;
    if (pos.z > upper_lim.z)
        new_pos.z = upper_lim.z;

    return new_pos;
}

TEST(MPI_ParticleExchanger, Exhange_Non_Periodical_One_To_Many_3_3_3)
{
    // MPI Data
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    // Topology
    pfc::Int3 sections(3, 3, 3);

    int required_size = sections.x * sections.y * sections.z;
    if (mpi_size < required_size)
        GTEST_SKIP() << "Not enough processes " << mpi_size << "/" << required_size;

    std::shared_ptr<mpi::Topology> topology = std::make_shared<mpi::Topology>(sections, mpi::Topology::LoopType::None, mpi_size);
    int center_rank = 13;

    // Particle exchanger
    mpi::ParticleExchanger<mpi::ParticleData3D_Uniform_SoA> particleExchanger;

    // Creating particle type & array
    MPI_Datatype particle_type = mpi::ParticleUtils<pfc::Three>::defineParticleType();
    pfc::ParticleArray3d particle_array;    

    // Config
    size_t num_particles = 50;
    pfc::FP3 lower_bound(0, 0, 0);
    pfc::FP3 upper_bound(1, 1, 1);

    pfc::FP3 center_lower_bound = pfc::FP3(0.25, 0.25, 0.25);
    pfc::FP3 center_upper_bound = pfc::FP3(0.75, 0.75, 0.75);

    if (mpi_rank == center_rank)
    {
        lower_bound = center_lower_bound;
        upper_bound = center_upper_bound;
    }

    else
    {
        pfc::Int3 rank_offset = cubeCornerIdToOffset(mpi_rank);

        // x
        if (rank_offset.x < 0)
        {
            lower_bound.x = 0;
            upper_bound.x = center_lower_bound.x;
        }
        if (rank_offset.x > 0)
        {
            lower_bound.x = center_upper_bound.x;
            upper_bound.x = 1.0;
        }

        // y
        if (rank_offset.y < 0)
        {
            lower_bound.y = 0;
            upper_bound.y = center_lower_bound.y;
        }
        if (rank_offset.y > 0)
        {
            lower_bound.y = center_upper_bound.y;
            upper_bound.y = 1.0;
        }

        // z
        if (rank_offset.z < 0)
        {
            lower_bound.z = 0;
            upper_bound.z = center_lower_bound.z;
        }
        if (rank_offset.z > 0)
        {
            lower_bound.z = center_upper_bound.z;
            upper_bound.z = 1.0;
        }
    }

    // Fill with particles
    if (mpi_rank == center_rank)
        for (size_t i = 0; i < num_particles; i++)
        {
            ParticleWriteType p;
            p.setVelocity(pfc::FP3(mpi_rank + 1, mpi_rank + 1, mpi_rank + 1));
            p.setPosition(pfc::FP3( (float)(rand()) / RAND_MAX, 
                                    (float)(rand()) / RAND_MAX, 
                                    (float)(rand()) / RAND_MAX));
            p.setWeight((mpi_rank + 1) * 5);
            particle_array.pushBack(p);
        }

    // Communication
    mpi::ParticleData3D_Uniform_SoA send_data[27];
    mpi::ParticleData3D_Uniform_SoA recv_data[27];

    int expected[27] = {};
    mpi::ParticleStorageInterface3D_Uniform_SoA::collectParticles(&particle_array, send_data, lower_bound, upper_bound, mpi_rank, topology, expected);

    particleExchanger.exchangeParticles(send_data, recv_data, mpi_rank, topology);

    mpi::ParticleStorageInterface3D_Uniform_SoA::insertParticles(&particle_array, recv_data, lower_bound, upper_bound, mpi_rank, topology);

    // Validation

    // Total number of particles
    int final_num_particles = 0;
    
    int local_count = particle_array.size();
    MPI_Reduce(&local_count, &final_num_particles, 1, MPI_INT, MPI_SUM, center_rank, topology->getTopologyCommunicator());
    if (mpi_rank == center_rank)
    {
        ASSERT_EQ(num_particles, final_num_particles);
    }

    // All particles withing their rank's domain
    for (int i = 0; i < particle_array.size(); i++)
    {
        ParticleReadType p = particle_array[i];

        ASSERT_TRUE(lower_bound.x <= p.getPosition().x && p.getPosition().x <= upper_bound.x);
        ASSERT_TRUE(lower_bound.y <= p.getPosition().y && p.getPosition().y <= upper_bound.y);
        ASSERT_TRUE(lower_bound.z <= p.getPosition().z && p.getPosition().z <= upper_bound.z);
    }
}

TEST(MPI_ParticleExchanger, Exhange_Non_Periodical_All_3_3_3)
{
    // MPI Data
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    // Topology
    pfc::Int3 sections(3, 3, 3);

    int required_size = sections.x * sections.y * sections.z;
    if (mpi_size < required_size)
        GTEST_SKIP() << "Not enough processes " << mpi_size << "/" << required_size;

    std::shared_ptr<mpi::Topology> topology = std::make_shared<mpi::Topology>(sections, mpi::Topology::LoopType::None, mpi_size);
    int center_rank = 13;

    // Particle exchanger
    mpi::ParticleExchanger<mpi::ParticleData3D_Uniform_SoA> particleExchanger;

    // Creating particle type & array
    MPI_Datatype particle_type = mpi::ParticleUtils<pfc::Three>::defineParticleType();
    pfc::ParticleArray3d particle_array;    

    // Config
    size_t num_particles = 50;

    pfc::FP3 generator_lower_bound(-1, -1, -1);
    pfc::FP3 generator_upper_bound(2, 2, 2);
    pfc::FP generator_out_of_bounds_bias = 1.25;

    pfc::FP3 lower_bound(0, 0, 0);
    pfc::FP3 upper_bound(1, 1, 1);

    pfc::Int3 rank_offset = cubeCornerIdToOffset(mpi_rank);
    pfc::FP3 rank_offset_fp(rank_offset.x, rank_offset.y, rank_offset.z);
    lower_bound += rank_offset_fp;
    upper_bound += rank_offset_fp;

    // Fill with particles
    for (size_t i = 0; i < num_particles; i++)
    {
        ParticleWriteType p;
        p.setVelocity(pfc::FP3(mpi_rank + 1, mpi_rank + 1, mpi_rank + 1));
        pfc::FP3 pos = lower_bound + generator_out_of_bounds_bias * (upper_bound - lower_bound) 
            * pfc::FP3( (float)(rand()) / RAND_MAX, 
                        (float)(rand()) / RAND_MAX, 
                        (float)(rand()) / RAND_MAX);
        p.setPosition(limitPosition(pos, generator_lower_bound, generator_upper_bound));
        
        p.setWeight(mpi_rank);
        particle_array.pushBack(p);
    }

    // Communication
    mpi::ParticleData3D_Uniform_SoA send_data[27];
    mpi::ParticleData3D_Uniform_SoA recv_data[27];

    int expected[27] = {};
    mpi::ParticleStorageInterface3D_Uniform_SoA::collectParticles(&particle_array, send_data, lower_bound, upper_bound, mpi_rank, topology, expected);

    particleExchanger.exchangeParticles(send_data, recv_data, mpi_rank, topology);

    mpi::ParticleStorageInterface3D_Uniform_SoA::insertParticles(&particle_array, recv_data, lower_bound, upper_bound, mpi_rank, topology);

    // Validation

    // Total number of particles
    int final_num_particles = 0;
    
    int local_count = particle_array.size();
    MPI_Reduce(&local_count, &final_num_particles, 1, MPI_INT, MPI_SUM, center_rank, topology->getTopologyCommunicator());
    if (mpi_rank == center_rank)
    {
        ASSERT_EQ(num_particles * required_size, final_num_particles);
    }

    // All particles withing their rank's domain
    for (int i = 0; i < particle_array.size(); i++)
    {
        ParticleReadType p = particle_array[i];

        ASSERT_TRUE(lower_bound.x <= p.getPosition().x && p.getPosition().x <= upper_bound.x);
        ASSERT_TRUE(lower_bound.y <= p.getPosition().y && p.getPosition().y <= upper_bound.y);
        ASSERT_TRUE(lower_bound.z <= p.getPosition().z && p.getPosition().z <= upper_bound.z);
    }
}

TEST(MPI_ParticleExchanger, Exhange_Periodical_All_3_3_3)
{
    // MPI Data
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    // Topology
    pfc::Int3 sections(3, 3, 3);

    int required_size = sections.x * sections.y * sections.z;
    if (mpi_size < required_size)
        GTEST_SKIP() << "Not enough processes " << mpi_size << "/" << required_size;

    std::shared_ptr<mpi::Topology> topology = std::make_shared<mpi::Topology>(sections, mpi::Topology::LoopType::LoopXYZ, mpi_size);
    int center_rank = 13;

    // Particle exchanger
    mpi::ParticleExchanger<mpi::ParticleData3D_Uniform_SoA> particleExchanger;

    // Creating particle type & array
    MPI_Datatype particle_type = mpi::ParticleUtils<pfc::Three>::defineParticleType();
    pfc::ParticleArray3d particle_array;    

    // Config
    size_t num_particles = 50;

    pfc::FP generator_out_of_bounds_bias = 1.25;

    pfc::FP3 lower_bound(0, 0, 0);
    pfc::FP3 upper_bound(1, 1, 1);

    pfc::Int3 rank_offset = cubeCornerIdToOffset(mpi_rank);
    pfc::FP3 rank_offset_fp(rank_offset.x, rank_offset.y, rank_offset.z);
    lower_bound += rank_offset_fp;
    upper_bound += rank_offset_fp;

    // Fill with particles
    for (size_t i = 0; i < num_particles; i++)
    {
        ParticleWriteType p;
        p.setVelocity(pfc::FP3(mpi_rank + 1, mpi_rank + 1, mpi_rank + 1));
        pfc::FP3 pos = lower_bound + generator_out_of_bounds_bias * (upper_bound - lower_bound) 
            * pfc::FP3( (float)(rand()) / RAND_MAX, 
                        (float)(rand()) / RAND_MAX, 
                        (float)(rand()) / RAND_MAX);
        p.setPosition(pos);
        
        p.setWeight(mpi_rank);
        particle_array.pushBack(p);
    }

    // Communication
    mpi::ParticleData3D_Uniform_SoA send_data[27];
    mpi::ParticleData3D_Uniform_SoA recv_data[27];

    int expected[27] = {};
    mpi::ParticleStorageInterface3D_Uniform_SoA::collectParticles(&particle_array, send_data, lower_bound, upper_bound, mpi_rank, topology, expected);

    particleExchanger.exchangeParticles(send_data, recv_data, mpi_rank, topology);

    mpi::ParticleStorageInterface3D_Uniform_SoA::insertParticles(&particle_array, recv_data, lower_bound, upper_bound, mpi_rank, topology);

    // Validation

    // Total number of particles
    int final_num_particles = 0;
    
    int local_count = particle_array.size();
    MPI_Reduce(&local_count, &final_num_particles, 1, MPI_INT, MPI_SUM, center_rank, topology->getTopologyCommunicator());
    if (mpi_rank == center_rank)
    {
        ASSERT_EQ(num_particles * required_size, final_num_particles);
    }

    // All particles withing their rank's domain
    for (int i = 0; i < particle_array.size(); i++)
    {
        ParticleReadType p = particle_array[i];

        ASSERT_TRUE(lower_bound.x <= p.getPosition().x && p.getPosition().x <= upper_bound.x);
        ASSERT_TRUE(lower_bound.y <= p.getPosition().y && p.getPosition().y <= upper_bound.y);
        ASSERT_TRUE(lower_bound.z <= p.getPosition().z && p.getPosition().z <= upper_bound.z);
    }
}

// p.setPosition(lower_bound  + generator_out_of_bounds_bias * 
//     (upper_bound - lower_bound) * pfc::FP3( (float)(rand()) / RAND_MAX, 
//                                             (float)(rand()) / RAND_MAX, 
//                                             (float)(rand()) / RAND_MAX));

// for (int r = 0; r < mpi_size; r++)
//     {
//         MPI_Barrier(topology->getTopologyCommunicator());
//         if (mpi_rank == r)
//         {
//             std::cout << "{RANK " << mpi_rank << "} Boundaries" << lower_bound.toString() << " " << upper_bound.toString() << std::endl;
//             for (int i = 0; i < particle_array.size(); i++)
//             {
//                 ParticleReadType p = particle_array[i];

//                 std::cout << "{RANK " << mpi_rank << "} Particle: " << p.getPosition().toString() << " " << p.getWeight() << std::endl;

//                 EXPECT_TRUE(lower_bound.x <= p.getPosition().x && p.getPosition().x <= upper_bound.x);
//                 EXPECT_TRUE(lower_bound.y <= p.getPosition().y && p.getPosition().y <= upper_bound.y);
//                 EXPECT_TRUE(lower_bound.z <= p.getPosition().z && p.getPosition().z <= upper_bound.z);
//             }
//         }
//         MPI_Barrier(topology->getTopologyCommunicator());
//     }