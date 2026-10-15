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
private:
    // Private constructor to prevent instancing
    ParticleExchanger() {}

private:
    static const int PARTICLE_HEADER_TAG_BASE = 100;
    static const int PARTICLE_PAYLOAD_TAG_BASE = 200;

public:

    // Exchanges out-of-bound particles between MPI ranks
    // send_data is sent to other ranks
    // All data receinved from other ranks is put into recv_data
    // In send_data and recv_data (0, 0, 0) offset is mapped to 13, so it will always be emtpy
    static void exchangeParticles(  const ParticleDataType send_data[27], ParticleDataType recv_data[27],
                                    int mpi_rank, std::shared_ptr<class Topology> topology, MPI_Comm communicator)
    {
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
                MPI_Isend(  send_data[dir].header, ParticleDataType::headerMpiSize, ParticleDataType::headerMpiDataType, 
                            send_neighbor, PARTICLE_HEADER_TAG_BASE + dir, communicator, &requests.back());
            }
            
            // Receiving
            if (recv_neighbor != MPI_INVALID_RANK)
            {
                requests.emplace_back();
                MPI_Irecv(  recv_data[dir].header, ParticleDataType::headerMpiSize, ParticleDataType::headerMpiDataType, 
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
                MPI_Isend(  send_data[dir].payload, ParticleDataType::payloadMpiDataType, send_data[dir].getPayloadMpiSize(), 
                            send_neighbor, PARTICLE_PAYLOAD_TAG_BASE + dir, communicator, &requests.back());
            }
            
            // Receiving
            if (recv_neighbor != MPI_INVALID_RANK && recv_data[dir].getCount() > 0)
            {
                requests.emplace_back();
                MPI_Irecv(  recv_data[dir].payload, ParticleDataType::payloadMpiDataType, recv_data[dir].getPayloadMpiSize(), 
                            recv_neighbor, PARTICLE_PAYLOAD_TAG_BASE + dir, communicator, &requests.back());
            };
        }

        // Waiting for communication to complete
        statuses.resize(requests.size());
        MPI_Waitall(requests.size(), requests.data(), statuses.data());

        // Make sure everyone has completed payload exchange
        MPI_Barrier(communicator);
    }
};

}