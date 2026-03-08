#include "Constants.h"
#include "FieldSolver.h"
#include "Grid.h"
#include "PmlFdtd.h"
#include "FieldBoundaryConditionFdtd.h"
#include "FieldGeneratorFdtd.h"
#include "Vectors.h"

namespace pfc {
using GridType = YeeGrid;

 FP3 eFunc(FP3 coords, FP t) 
 {
    FP3 e;
    for (int i = 0; i < 3; i++)
        e[i] = sin((coords.x + coords.y) * coords.z + i*i);
    return e;
}

FP3 bFunc(FP3 coords, FP t) 
{
    FP3 e;
    for (int i = 0; i < 3; i++)
        e[i] = sin((coords.x + coords.y) * coords.z + i*i);
    return e;
}

void initializeGrid(std::unique_ptr<GridType>& grid) 
{
    Int3 begin = Int3(0, 0, 0);
    Int3 end = grid->numCells;

    for (int i = begin.x; i < end.x; i++)
        for (int j = begin.y; j < end.y; j++)
            for (int k = begin.z; k < end.z; k++) {
                grid->Ex(i, j, k) = eFunc(grid->ExPosition(i, j, k), 0).x;
                grid->Ey(i, j, k) = eFunc(grid->EyPosition(i, j, k), 0).y;
                grid->Ez(i, j, k) = eFunc(grid->EzPosition(i, j, k), 0).z;

                grid->Bx(i, j, k) = bFunc(grid->BxPosition(i, j, k), 0).x;
                grid->By(i, j, k) = bFunc(grid->ByPosition(i, j, k), 0).y;
                grid->Bz(i, j, k) = bFunc(grid->BzPosition(i, j, k), 0).z;
            }
}

int main()
{
    Int3 gridSize = Int3(8, 8, 8);
    FP3 minCoords = FP3(0, 0, 0);
    FP3 maxCoords = FP3(1, 1, 1);
    FP3 gridStep = (maxCoords - minCoords) / (FP3)gridSize;

    std::unique_ptr<GridType> grid;
    grid.reset(new GridType(gridSize, minCoords, gridStep, gridSize));

    initializeGrid(grid);

    std::cout << std::endl << "Hello Grid!" << std::endl;

    return 0;
}

}