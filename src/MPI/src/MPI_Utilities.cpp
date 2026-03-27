#include "MPI_Utilities.h"

namespace mpi
{

// Returns the inverse of the specified direction
Direction invertDirection(const Direction& direction)
{
    switch (direction)
    {
    case Direction::positiveX: return Direction::negativeX;
    case Direction::negativeX: return Direction::positiveX;
    case Direction::positiveY: return Direction::negativeY;
    case Direction::negativeY: return Direction::positiveY;
    case Direction::positiveZ: return Direction::negativeZ;
    case Direction::negativeZ: return Direction::positiveZ;

    default: return Direction::positiveX;
    }
}

}