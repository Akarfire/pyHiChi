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
template<class ParticleDataType>
class ParticleExchanger final
{
public:
    ParticleExchanger() 
    {
        // Creating custom mpi types if necessary
        if (ParticleDataType::useCustomHeaderMpiDatatype || ParticleDataType::useCustomPayloadMpiDatatype)
        {
            ParticleDataType::createCustomMpiDatatypes(customMpiDatatypeBuffer);

            // Resolving types
            // Header type
            if (ParticleDataType::useCustomHeaderMpiDatatype)
                headerType = customMpiDatatypeBuffer[ParticleDataType::headerCustomMpiDatatypeIndex];
            else
                headerType = ParticleDataType::headerMpiDataType;

            // Payload type
            if (ParticleDataType::useCustomPayloadMpiDatatype)
                payloadType = customMpiDatatypeBuffer[ParticleDataType::payloadCustomMpiDatatypeIndex];
            else
                payloadType = ParticleDataType::payloadMpiDataType;
        }
    }

    ~ParticleExchanger()
    {
        // Freeing custom mpi data types
        int mpiInit;
        MPI_Initialized(&mpiInit);

        int mpiFinalized;
        MPI_Finalized(&mpiFinalized);

        if (mpiInit && !mpiFinalized)
            for (int i = 0; i < customMpiDatatypeBuffer.size(); i++)
                MPI_Type_free(&customMpiDatatypeBuffer[i]);
    }

private:
    static const int PARTICLE_HEADER_TAG_BASE = 100;
    static const int PARTICLE_PAYLOAD_TAG_BASE = 200;

    // Custom mpi types are created in ParticleData, but stored in particle exchanger
    std::vector<MPI_Datatype> customMpiDatatypeBuffer;

    // Resolved types
    MPI_Datatype headerType = -1;
    MPI_Datatype payloadType = -1;

public:

    // Exchanges out-of-bound particles between MPI ranks
    // send_data is sent to other ranks
    // All data receinved from other ranks is put into recv_data
    // In send_data and recv_data (0, 0, 0) offset is mapped to 13, so it will always be emtpy
    static void exchangeParticles(  const ParticleDataType send_data[27], ParticleDataType recv_data[27],
                                    int mpi_rank, std::shared_ptr<class Topology> topology)
    {
        MPI_Comm communicator = topology->getTopologyCommunicator();

        std::vector<MPI_Request> requests;
        requests.reserve(26 * 2); // 1 send request and 1 recv requst for each direction
        std::vector<MPI_Status> statuses;

        // Exchanging headers
        for (int dir = 0; dir < 27; dir++)
        {
            if (dir == NO_DIRECTION) continue; // Skipping the (0, 0, 0) offset (it is not used)

            pfc::Int3 send_offset = cubeCornerIdToOffset(dir);
            pfc::Int3 recv_offset = -1 * send_offset;

            int send_neighbor = topology->getNeighbor(mpi_rank, send_offset);
            int recv_neighbor = topology->getNeighbor(mpi_rank, recv_offset);

            // Sending
            if (send_neighbor != MPI_INVALID_RANK)
            {
                requests.emplace_back();
                MPI_Isend(  send_data[dir].header, ParticleDataType::headerMpiSize, headerType, 
                            send_neighbor, PARTICLE_HEADER_TAG_BASE + dir, communicator, &requests.back());
            }
            
            // Receiving
            if (recv_neighbor != MPI_INVALID_RANK)
            {
                requests.emplace_back();
                MPI_Irecv(  recv_data[dir].header, ParticleDataType::headerMpiSize, headerType, 
                            recv_neighbor, PARTICLE_HEADER_TAG_BASE + dir, communicator, &requests.back());
            };
        }

        // Waiting for communication to complete
        statuses.resize(requests.size());
        MPI_Waitall(requests.size(), requests.data(), statuses.data());

        // Make sure everyone has completed header exchange
        MPI_Barrier(communicator);

        // Exchanging payload
        requests.clear();

        // Resolving type
        for (int dir = 0; dir < 27; dir++)
        {
            if (dir == NO_DIRECTION) continue; // Skipping the (0, 0, 0) offset (it is not used)

            pfc::Int3 send_offset = cubeCornerIdToOffset(dir);
            pfc::Int3 recv_offset = -1 * send_offset;

            int send_neighbor = topology->getNeighbor(mpi_rank, send_offset);
            int recv_neighbor = topology->getNeighbor(mpi_rank, recv_offset);

            // Sending
            if (send_neighbor != MPI_INVALID_RANK && send_data[dir].getCount() > 0)
            {
                requests.emplace_back();
                MPI_Isend(  send_data[dir].payload, send_data[dir].getPayloadMpiSize(), payloadType, 
                            send_neighbor, PARTICLE_PAYLOAD_TAG_BASE + dir, communicator, &requests.back());
            }
            
            // Receiving
            if (recv_neighbor != MPI_INVALID_RANK && recv_data[dir].getCount() > 0)
            {
                requests.emplace_back();
                MPI_Irecv(  recv_data[dir].payload, recv_data[dir].getPayloadMpiSize(), payloadType, 
                            recv_neighbor, PARTICLE_PAYLOAD_TAG_BASE + dir, communicator, &requests.back());
            };
        }

        // Waiting for communication to complete
        statuses.resize(requests.size());
        MPI_Waitall(requests.size(), requests.data(), statuses.data());
    }
};

}