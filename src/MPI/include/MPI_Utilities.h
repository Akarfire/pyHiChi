#pragma once

#include <mpi.h>
#include <vector>

#include "FP.h"
#include "Vectors.h"

namespace mpi
{

// Invalid rank id
static int MPI_INVALID_RANK = -1;

// Invalid topological section value
static pfc::Int3 MPI_INVALID_SECTION = pfc::Int3(MPI_INVALID_RANK, MPI_INVALID_RANK, MPI_INVALID_RANK);

// Same rank (zero offset)
static pfc::Int3 MPI_SAME_RANK = pfc::Int3(0, 0, 0);

// Defines the direction, in which the transmission will be performed (From sender to receiver)
enum class Direction
{
    pX, nX,
    pY, nY,
    pZ, nZ
};

// Conversion map to turn mpi::Direction
// use like this mpi::DirectionToOffset[static_cast<int>(direction)]
static const pfc::Int3 DirectionToOffset[6] = 
{
    pfc::Int3{1, 0, 0}, pfc::Int3{-1, 0, 0},
    pfc::Int3{0, 1, 0}, pfc::Int3{0, -1, 0},
    pfc::Int3{0, 0, 1}, pfc::Int3{0, 0, -1}
};

// Returns the inverse of the specified direction
Direction invertDirection(const Direction& direction);


// Converts (0|+-1, 0|+-1, 0|+-1) type offsets to 0 - 26 indices
static int offsetToCubeCornerID(int x, int y, int z);

// Converts (0|+-1, 0|+-1, 0|+-1) type offsets to 0 - 26 indices
static int offsetToCubeCornerID(const pfc::Int3& offset);

// Converts 0 - 26 indices to (0|+-1, 0|+-1, 0|+-1) type offsets 
static pfc::Int3 cubeCornerIdToOffset(int id);

}