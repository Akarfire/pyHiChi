#include "MPI_ParticleExchanger.h"

#include "mpi.h"

#include "MPI_Topology.h"

#include "ParticleTypes.h"

namespace mpi
{

// Returns the cube corner id corresponding to the neighboring rank, that a particle 
// with the specified position shall be sent to
int ParticleExchanger::getResponsibleNeighborID(   FP x, FP y, FP z, 
                                const pfc::FP3& lower_bound, 
                                const pfc::FP3& upper_bound
                                )
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

// Returns an offset corresponding to the neighbouring rank, that a particle 
// with the specified position shall be sent to
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


// Prepares ParticleData_SoA for every direction
void ParticleExchanger::prepareSendData(ParticleData3D out_data[27], std::unique_ptr<Ensemble3d>& ensemble, 
                                        const pfc::FP3& lower_bound, const pfc::FP3& upper_bound)
{
    // Send indices for every direction
    // (0, 0, 0) offset is mapped to 13, so it will always be emtpy
    std::vector<int> send_indices[27][pfc::sizeParticleTypes]; // TO DO : Add a heuristic reservation mechanism
    std::vector<int> deletion_indices[pfc::sizeParticleTypes];

    Ensemble3d& ens = (*ensemble);

    // Loop over the ensemble and select the out-of-bound particles
    for (int type = 0; type < pfc::sizeParticleTypes; type++)
    {
        auto& pos_x = ens[type].getPositions(0);
        auto& pos_y = ens[type].getPositions(1);
        auto& pos_z = ens[type].getPositions(2);

        for (int i = 0; i < ens[type].size(); i++)
        {
            int neighbor_id = getResponsibleNeighborID(pos_x[i], pos_y[i], pos_z[i], lower_bound, upper_bound);
            if (neighbor_id != NO_DIRECTION)
            {
                send_indices[neighbor_id][type].push_back(i);
                deletion_indices[type].push_back(i);
            }
        }
    }

    // Packaging loop
    for (int dir = 0; dir < 27; dir++)
    {
        if (dir == NO_DIRECTION) continue; // Skipping the (0, 0, 0) offset (it is not used)

        // Filling out the header
        out_data[dir].count = 0;
        for (int type = 0; type < pfc::sizeParticleTypes; type++)
        {
            int type_count = static_cast<int>(send_indices[dir][type].size());
            out_data[dir].count += type_count;
            out_data[dir].type_counts[type] = type_count;
        }

        // Reserving space for data
        int count = out_data[dir].count;
        out_data[dir].positions[0].reserve(count);
        out_data[dir].positions[1].reserve(count);
        out_data[dir].positions[2].reserve(count);

        out_data[dir].momentums[0].reserve(count);
        out_data[dir].momentums[1].reserve(count);
        out_data[dir].momentums[2].reserve(count);

        out_data[dir].weights.reserve(count);
        out_data[dir].gammas.reserve(count);  

        // Package positions
        for (int type = 0; type < pfc::sizeParticleTypes; type++)
        {
            std::vector<int>& indices = send_indices[dir][type];

            std::vector<Real>& positions_x = ens[type].getPositions(0);
            std::vector<Real>& positions_y = ens[type].getPositions(1);
            std::vector<Real>& positions_z = ens[type].getPositions(2);

            for (int i = 0; i < indices.size(); i++)
            {
                int index = indices[i];
                out_data[dir].positions[0].push_back(positions_x[index]);
                out_data[dir].positions[1].push_back(positions_y[index]);
                out_data[dir].positions[2].push_back(positions_z[index]);
            }
        }

        // Package momentums
        for (int type = 0; type < pfc::sizeParticleTypes; type++)
        {
            std::vector<Real>& momentums_x = ens[type].getMomentums(0);
            std::vector<Real>& momentums_y = ens[type].getMomentums(1);
            std::vector<Real>& momentums_z = ens[type].getMomentums(2);

            std::vector<int>& indices = send_indices[dir][type];
            for (int i = 0; i < indices.size(); i++)
            {
                int index = indices[i];
                out_data[dir].momentums[0].push_back(momentums_x[index]);
                out_data[dir].momentums[1].push_back(momentums_y[index]);
                out_data[dir].momentums[2].push_back(momentums_z[index]);
            }
        }

        // Package weights
        for (int type = 0; type < pfc::sizeParticleTypes; type++)
        {
            std::vector<int>& indices = send_indices[dir][type];
            for (int i = 0; i < indices.size(); i++)
            {
                int index = indices[i];
                out_data[dir].weights.push_back(ens[type].getWeights()[index]);
            }
        }

        // Package gammas
        for (int type = 0; type < pfc::sizeParticleTypes; type++)
        {
            std::vector<int>& indices = send_indices[dir][type];
            for (int i = 0; i < indices.size(); i++)
            {
                int index = indices[i];
                out_data[dir].gammas.push_back(ens[type].getGammas()[index]);
            }
        }
    }

    // Deletion loop
    for (int type = 0; type < pfc::sizeParticleTypes; type++)
    {
        auto& del = deletion_indices[type];
        for (int i = static_cast<int>(del.size()) - 1; i >= 0; i--)
        {
            ens[type].deleteParticle(del[i]);
        }
    }
}

// Exchanges out-of-bound particles between MPI ranks
void ParticleExchanger::exchangeParticles(  int mpi_rank, std::unique_ptr<Ensemble3d>& ensemble, 
                                            const pfc::FP3& lower_bound, const pfc::FP3& upper_bound,
                                            std::shared_ptr<Topology> topology, MPI_Comm communicator)
{
    ParticleData3D send_data[27]; // (0, 0, 0) offset is mapped to 13, so it will always be emtpy
    prepareSendData(send_data, ensemble, lower_bound, upper_bound);

    std::vector<MPI_Request> requests;
    requests.reserve(26 * 2); // 1 send request and 1 recv requst for each direction
    std::vector<MPI_Status> statuses;

    // Exchanging headers
    constexpr int header_size = pfc::sizeParticleTypes + 1;
    int send_headers[27][header_size];
    int recv_headers[27][header_size] = {};

    for (int dir = 0; dir < 27; dir++)
    {
        if (dir == NO_DIRECTION) continue; // Skipping the (0, 0, 0) offset (it is not used)

        pfc::Int3 send_offset = cubeCornerIdToOffset(dir);
        pfc::Int3 recv_offset = -1 * send_offset;

        int send_neighbor = topology->getNeighbor(mpi_rank, send_offset);
        int recv_neighbor = topology->getNeighbor(mpi_rank, recv_offset);

        // Header send data
        send_headers[dir][0] = send_data[dir].count;
        for (int i = 0; i < pfc::sizeParticleTypes; i++)
            send_headers[dir][i + 1] = send_data[dir].type_counts[i];

        // Sending
        if (send_neighbor != MPI_INVALID_RANK)
        {
            requests.emplace_back();
            MPI_Isend(send_headers[dir], header_size, MPI_INT, send_neighbor, PARTICLE_HEADER_TAG_BASE + dir, communicator, &requests.back());
        }
        
        // Receiving
        if (recv_neighbor != MPI_INVALID_RANK)
        {
            requests.emplace_back();
            MPI_Irecv(recv_headers[dir], header_size, MPI_INT, recv_neighbor, PARTICLE_HEADER_TAG_BASE + dir, communicator, &requests.back());
        };
    }

    // Waiting for communication to complete
    statuses.resize(requests.size());
    MPI_Waitall(requests.size(), requests.data(), statuses.data());

    // Defining recv_data
    ParticleData3D recv_data[27];

    // Unpacking recv_header & reserve space for recv data
    for (int dir = 0; dir < 27; dir++)
    {
        if (dir == NO_DIRECTION) continue; // Skipping the (0, 0, 0) offset (it is not used)

        pfc::Int3 recv_offset = -1 * cubeCornerIdToOffset(dir);
        int recv_neighbor = topology->getNeighbor(mpi_rank, recv_offset);

        if (recv_neighbor == MPI_INVALID_RANK) continue; // Skipping ivalid neighbors

        recv_data[dir].count = recv_headers[dir][0];
        for (int i = 0; i < pfc::sizeParticleTypes; i++)
            recv_data[dir].type_counts[i] = recv_headers[dir][i + 1];

        recv_data[dir].positions[0].resize(recv_data[dir].count);
        recv_data[dir].positions[1].resize(recv_data[dir].count);
        recv_data[dir].positions[2].resize(recv_data[dir].count);
    }

    // ...
}

}