#pragma once

#include <mpi.h>
#include <vector>

#include "FP.h"
#include "Vectors.h"

#define MPI_INVALID_RANK -1

namespace mpi
{

// Defines the direction, in which the transmission will be performed (From sender to receiver)
enum class Direction
{
    positiveX,
    negativeX,
    positiveY,
    negativeY,
    positiveZ,
    negativeZ
};

// Returns the inverse of the specified direction
Direction invertDirection(const Direction& direction);

}