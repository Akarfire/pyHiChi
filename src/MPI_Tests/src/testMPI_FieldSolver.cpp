#include "TestingUtility.h"

#include "gtest/gtest.h"

#include "Fdtd.h"
#include "Pstd.h"
#include "Psatd.h"
#include "PsatdTimeStaggered.h"

#include <iostream>
#include <fstream>

#include "HiChi_MPI.h"

template <class TTypeDefinitionsFieldTest>
class MPI_FieldSolverTest : public BaseFixture {
public:
    using FieldSolverType = typename TTypeDefinitionsFieldTest::FieldSolverType;
    using GridType = typename TTypeDefinitionsFieldTest::FieldSolverType::GridType;

    const int dimension = TTypeDefinitionsFieldTest::dimension;
    const CoordinateEnum axis = TTypeDefinitionsFieldTest::axis;

    const int gridSizeLongitudinal = 32;
    const int gridSizeTransverse = 8;

    std::unique_ptr<FieldSolverType> fieldSolver;
    std::unique_ptr<GridType> grid;
    Int3 gridSize;
    FP3 gridStep;
    FP3 mainMinCoords = FP3(0, 0, 0);
    FP3 mainMaxCoords = FP3(0, 0, 0);
    FP timeStep = 0;
    FP3 minCoords, maxCoords;
    int numSteps = 0;

    const FP maxError = 1e-2;

    // MPI
    std::unique_ptr<mpi::Topology> topology;
    std::unique_ptr<mpi::FieldExchanger> fieldExchanger;

    virtual void SetUp() {
        int mpi_rank;
        int mpi_size;
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

        Int3 sections = Int3(4, 1, 1);
        int topologySize = sections.x * sections.y * sections.z;

        if (topologySize > mpi_size)
            GTEST_SKIP() << "Not enough processes for FieldSolverTest";

        if (mpi_rank >= topologySize)
            GTEST_SKIP() << "Rank " << mpi_rank << " is not required";

        topology.reset(new mpi::Topology(sections, mpi::Topology::LoopType::None, mpi_size));

        Int3 mainGridSize = Int3(1, 1, 1);
        for (int d = 0; d < dimension; d++) {
            mainGridSize[d] = gridSizeTransverse;
        }
        mainGridSize[(int)axis] = gridSizeLongitudinal;

        mainMinCoords = FP3(0, 0, 0);
        mainMaxCoords = constants::c * (FP3)mainGridSize;
        this->gridStep = (mainMaxCoords - mainMinCoords) / (FP3)mainGridSize;

        int div = mainGridSize.x / 4;
        std::vector<int> divisions[3] = {
            {div, 2 * div, 3 * div}, 
            {}, 
            {}
        };
        mpi::GridSlicer::getSubGridParameters(  this->minCoords, this->gridSize, 
                                                mainMinCoords, mainGridSize, 
                                                this->gridStep, divisions, 
                                                mpi_rank, *topology);

        // std::cout << "LOCAL SIZE: " << mpi_rank << " : " << this->gridSize.x << " " << this->gridSize.y << " " << this->gridSize.z << std::endl;
        // std::cout << "LOCAL ORIGIN: " << mpi_rank << " : " << this->minCoords.x << " " << this->minCoords.y << " " << this->minCoords.z << std::endl;

        // std::cout << "GLOBAL MIN: " << this->mainMinCoords.x << " " << this->mainMinCoords.y << " " << this->mainMinCoords.z << std::endl;
        // std::cout << "GLOBAL MAX: " << this->mainMaxCoords.x << " " << this->mainMaxCoords.y << " " << this->mainMaxCoords.z << std::endl;
        // std::cout << "GLOBAL SIZE: " << mainGridSize.x << " " << mainGridSize.y << " " << mainGridSize.z << std::endl;

        this->grid.reset(new GridType(this->gridSize, this->minCoords, this->gridStep, this->gridSize));

        fieldExchanger.reset(new mpi::FieldExchanger(grid->numCells, grid->numExternalCells));

        this->timeStep = 0.5 * FieldSolverType::getCourantConditionTimeStep(this->gridStep);
        this->numSteps = (int)((mainMaxCoords - mainMinCoords)[(int)axis] /
            (constants::c * this->timeStep) * 0.2);

        fieldSolver.reset(new FieldSolverType(this->grid.get(), this->timeStep));

        initializeGrid();
    }

    void initializeGrid() {
        Int3 begin = Int3(0, 0, 0);
        Int3 end = grid->numCells;

        for (int i = begin.x; i < end.x; i++)
            for (int j = begin.y; j < end.y; j++)
                for (int k = begin.z; k < end.z; k++) {
                    FP3 coords = grid->ExPosition(i, j, k);
                    grid->Ex(i, j, k) = eFunc(coords.x, coords.y, coords.z, 0).x;
                    coords = grid->EyPosition(i, j, k);
                    grid->Ey(i, j, k) = eFunc(coords.x, coords.y, coords.z, 0).y;
                    coords = grid->EzPosition(i, j, k);
                    grid->Ez(i, j, k) = eFunc(coords.x, coords.y, coords.z, 0).z;
                    coords = grid->BxPosition(i, j, k);
                    grid->Bx(i, j, k) = bFunc(coords.x, coords.y, coords.z, 0).x;
                    coords = grid->ByPosition(i, j, k);
                    grid->By(i, j, k) = bFunc(coords.x, coords.y, coords.z, 0).y;
                    coords = grid->BzPosition(i, j, k);
                    grid->Bz(i, j, k) = bFunc(coords.x, coords.y, coords.z, 0).z;
                }
    }

    FP fieldFunc(FP x, FP y, FP z, FP t) {
        FP3 coord(x, y, z);
        int axis0 = (int)this->axis;
        return sin((FP)2.0 * constants::pi / (mainMaxCoords[axis0] - mainMinCoords[axis0]) *
            (coord[axis0] - constants::c * t - mainMinCoords[axis0]));
    }

    FP3 eFunc(FP x, FP y, FP z, FP t) {
        CoordinateEnum axisE = CoordinateEnum(((int)axis + 1) % 3);
        FP3 e;
        e[(int)axisE] = fieldFunc(x, y, z, t);
        return e;
    }

    FP3 bFunc(FP x, FP y, FP z, FP t) {
        CoordinateEnum axisB = CoordinateEnum(((int)axis + 2) % 3);
        FP3 b;
        b[(int)axisB] = fieldFunc(x, y, z, t);
        return b;
    }
};

// #ifndef __USE_FFT__

typedef ::testing::Types <
    // TypeDefinitionsFieldTest<FDTD, 1, CoordinateEnum::x>,
    // TypeDefinitionsFieldTest<FDTD, 2, CoordinateEnum::x>,
    // TypeDefinitionsFieldTest<FDTD, 2, CoordinateEnum::y>,
    TypeDefinitionsFieldTest<FDTD, 3, CoordinateEnum::x>//,
    // TypeDefinitionsFieldTest<FDTD, 3, CoordinateEnum::y>,
    // TypeDefinitionsFieldTest<FDTD, 3, CoordinateEnum::z>
> types;

// #else

// typedef ::testing::Types <
//     // TypeDefinitionsFieldTest<FDTD, 1, CoordinateEnum::x>,
//     // TypeDefinitionsFieldTest<FDTD, 2, CoordinateEnum::x>,
//     // TypeDefinitionsFieldTest<FDTD, 2, CoordinateEnum::y>,
//     TypeDefinitionsFieldTest<FDTD, 3, CoordinateEnum::x>,
//     TypeDefinitionsFieldTest<FDTD, 3, CoordinateEnum::y>,
//     TypeDefinitionsFieldTest<FDTD, 3, CoordinateEnum::z>,

//     // TypeDefinitionsFieldTest<PSTD, 1, CoordinateEnum::x>,
//     // TypeDefinitionsFieldTest<PSTD, 2, CoordinateEnum::x>,
//     // TypeDefinitionsFieldTest<PSTD, 2, CoordinateEnum::y>,
//     TypeDefinitionsFieldTest<PSTD, 3, CoordinateEnum::x>,
//     TypeDefinitionsFieldTest<PSTD, 3, CoordinateEnum::y>,
//     TypeDefinitionsFieldTest<PSTD, 3, CoordinateEnum::z>,

//     // TypeDefinitionsFieldTest<PSATD, 1, CoordinateEnum::x>,
//     // TypeDefinitionsFieldTest<PSATD, 2, CoordinateEnum::x>,
//     // TypeDefinitionsFieldTest<PSATD, 2, CoordinateEnum::y>,
//     TypeDefinitionsFieldTest<PSATD, 3, CoordinateEnum::x>,
//     TypeDefinitionsFieldTest<PSATD, 3, CoordinateEnum::y>,
//     TypeDefinitionsFieldTest<PSATD, 3, CoordinateEnum::z>,

//     // TypeDefinitionsFieldTest<PSATDTimeStaggered, 1, CoordinateEnum::x>,
//     // TypeDefinitionsFieldTest<PSATDTimeStaggered, 2, CoordinateEnum::x>,
//     // TypeDefinitionsFieldTest<PSATDTimeStaggered, 2, CoordinateEnum::y>,
//     TypeDefinitionsFieldTest<PSATDTimeStaggered, 3, CoordinateEnum::x>,
//     TypeDefinitionsFieldTest<PSATDTimeStaggered, 3, CoordinateEnum::y>,
//     TypeDefinitionsFieldTest<PSATDTimeStaggered, 3, CoordinateEnum::z>
// > types;

// #endif

TYPED_TEST_CASE(MPI_FieldSolverTest, types);


TYPED_TEST(MPI_FieldSolverTest, PeriodicalFieldSolverTest)
{
    int mpi_rank;
    int mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    //this->fieldSolver->setPeriodicalBoundaryConditions();

    for (int step = 0; step < this->numSteps; ++step)
    {
        this->fieldSolver->updateFields();

        // MPI exchange sequence
        this->fieldExchanger->performExchangeSequence(grid->Ex.getData(), *(this->topology), mpi_rank, MPI_COMM_WORLD);
        this->fieldExchanger->performExchangeSequence(grid->Ey.getData(), *(this->topology), mpi_rank, MPI_COMM_WORLD);
        this->fieldExchanger->performExchangeSequence(grid->Ez.getData(), *(this->topology), mpi_rank, MPI_COMM_WORLD);

        this->fieldExchanger->performExchangeSequence(grid->Bx.getData(), *(this->topology), mpi_rank, MPI_COMM_WORLD);
        this->fieldExchanger->performExchangeSequence(grid->By.getData(), *(this->topology), mpi_rank, MPI_COMM_WORLD);
        this->fieldExchanger->performExchangeSequence(grid->Bz.getData(), *(this->topology), mpi_rank, MPI_COMM_WORLD);
    }

    FP finalT = this->fieldSolver->dt * this->numSteps;

    Int3 begin = this->fieldSolver->internalIndexBegin;
    Int3 end = this->fieldSolver->internalIndexEnd;

    std::string file_name = "debug_" + std::to_string(mpi_rank) + ".txt";
    std::ofstream file(file_name);
    if (mpi_rank == 0)
        file << "rank | x | y | z | result" << std::endl;

    for (int i = begin.x; i < end.x; ++i)
        for (int j = begin.y; j < end.y; ++j)
            for (int k = begin.z; k < end.z; ++k)
            {
                //int j = 0;

                FP3 expectedE, actualE;
                FP3 coords = this->grid->ExPosition(i, j, k);
                expectedE.x = this->eFunc(coords.x, coords.y, coords.z, finalT).x;
                coords = this->grid->EyPosition(i, j, k);
                expectedE.y = this->eFunc(coords.x, coords.y, coords.z, finalT).y;
                coords = this->grid->EzPosition(i, j, k);
                expectedE.z = this->eFunc(coords.x, coords.y, coords.z, finalT).z;
                actualE.x = this->grid->Ex(i, j, k);
                actualE.y = this->grid->Ey(i, j, k);
                actualE.z = this->grid->Ez(i, j, k);

                // if (abs(expectedE.norm() - actualE.norm()) > this->maxError)
                // {
                        file << mpi_rank << " | " << coords.x << " | " << coords.y << " | " << coords.z << " | " << actualE.x - expectedE.x << std::endl;
                // }
                    
                // else
                // {
                //     file << mpi_rank << " | " << coords.x << " | " << coords.y << " | " << coords.z << " | " << expectedE.x << std::endl;
                // }
                    
                EXPECT_NEAR(expectedE.norm(), actualE.norm(), this->maxError);
            }

    for (int i = begin.x; i < end.x; ++i)
        for (int j = begin.y; j < end.y; ++j)
            for (int k = begin.z; k < end.z; ++k)
            {
                //int j = 0;

                FP3 expectedB, actualB;
                FP3 coords = this->grid->BxPosition(i, j, k);
                expectedB.x = this->bFunc(coords.x, coords.y, coords.z, finalT).x;
                coords = this->grid->ByPosition(i, j, k);
                expectedB.y = this->bFunc(coords.x, coords.y, coords.z, finalT).y;
                coords = this->grid->BzPosition(i, j, k);
                expectedB.z = this->bFunc(coords.x, coords.y, coords.z, finalT).z;
                actualB.x = this->grid->Bx(i, j, k);
                actualB.y = this->grid->By(i, j, k);
                actualB.z = this->grid->Bz(i, j, k);
                EXPECT_NEAR(expectedB.norm(), actualB.norm(), this->maxError);

                //file << mpi_rank << " | " << coords.x << " | " << coords.y << " | " << coords.z << " | " << actualB.x - expectedB.x << std::endl;
            }
}