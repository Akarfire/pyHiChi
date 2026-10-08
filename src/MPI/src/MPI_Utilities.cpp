#include "MPI_Utilities.h"

namespace mpi
{

// Returns the inverse of the specified direction
Direction invertDirection(const Direction& direction)
{
    switch (direction)
    {
    case Direction::pX: return Direction::nX;
    case Direction::nX: return Direction::pX;
    case Direction::pY: return Direction::nY;
    case Direction::nY: return Direction::pY;
    case Direction::pZ: return Direction::nZ;
    case Direction::nZ: return Direction::pZ;

    default: return Direction::pX;
    }
}

// Converts (0|+-1, 0|+-1, 0|+-1) type offsets to 0 - 26 indices
int offsetToCubeCornerID(const pfc::Int3& offset)
{
    return (offset.x + 1) + (offset.y + 1) * 3 + (offset.z + 1) * 9;
}

// Converts 0 - 26 indices to (0|+-1, 0|+-1, 0|+-1) type offsets 
pfc::Int3 cubeCornerIdToOffset(int id)
{
    pfc::Int3 offset;

    offset.z = id / 9 - 1;
    offset.y = (id % 9) / 3 - 1;
    offset.x = id % 3 - 1;

    return offset;
}

}