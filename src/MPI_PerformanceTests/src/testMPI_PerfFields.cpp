#include <iostream>
#include <chrono>
#include "mpi.h"

#include "Fdtd.h"
#include "HiChi_MPI.h"

#define VERIFY

using namespace pfc;

class PerfTest
{
    // CONFIGURATION

    using FieldSolverType = typename FDTD;
    using GridType = typename YeeGrid;

    const int dimension = 3;
    const CoordinateEnum axis = CoordinateEnum::x;

    const int gridSizeLongitudinal = 128;
    const int gridSizeTransverse = 128;

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

    std::shared_ptr<mpi::Topology> topology;
    std::shared_ptr<mpi::FieldExchanger> fieldExchanger;

    // HELPERS

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
        int axis0 = (int)axis;
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

public:
    // Constructor / Runner
    PerfTest(const Int3& sections)
    {
        auto total_start = std::chrono::high_resolution_clock::now();

        // Setup 
        
        auto setup_start = std::chrono::high_resolution_clock::now();

        int mpi_size; 
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

        int topologySize = sections.x * sections.y * sections.z;
        if (topologySize > mpi_size)
        {
            std::cout << "Not enough processes for the test: " << mpi_size << "/" << topologySize << std::endl;
            return;
        }

        topology = std::make_shared<mpi::Topology>(sections, mpi::Topology::LoopType::LoopXYZ, mpi_size);

        if (!topology->isValidOnThisRank())
            return;

        int mpi_rank;
        MPI_Comm_rank(topology->getTopologyCommunicator(), &mpi_rank);

        Int3 mainGridSize = Int3(1, 1, 1);
        for (int d = 0; d < dimension; d++) 
        {
            mainGridSize[d] = gridSizeTransverse;
        }
        mainGridSize[(int)axis] = gridSizeLongitudinal;

        mainMinCoords = FP3(0, 0, 0);
        mainMaxCoords = constants::c * (FP3)mainGridSize;
        gridStep = (mainMaxCoords - mainMinCoords) / (FP3)mainGridSize;

        std::vector<int> divisions[3] = {
            {},
            {}, 
            {}
        };
        for (int i = 1; i < sections.x; i++)
            divisions[0].push_back((mainGridSize.x / sections.x) * i);

        for (int i = 1; i < sections.y; i++)
            divisions[1].push_back((mainGridSize.y / sections.y) * i);

        for (int i = 1; i < sections.z; i++)
            divisions[2].push_back((mainGridSize.z / sections.z) * i);

        Int3 localIndexOffset;
        mpi::GridSlicer::getSubGridParameters(  minCoords, gridSize, localIndexOffset,
                                                mainMinCoords, mainGridSize, 
                                                gridStep, divisions, 
                                                mpi_rank, topology);

        grid.reset(new GridType(gridSize, minCoords, gridStep, gridSize));

        fieldExchanger = std::make_shared<mpi::FieldExchanger>(grid->numCells, grid->numExternalCells);

        timeStep = 0.5 * FieldSolverType::getCourantConditionTimeStep(gridStep);
        numSteps = (int)((mainMaxCoords - mainMinCoords)[(int)axis] /
            (constants::c * timeStep) * 0.2);

        fieldSolver.reset(new FieldSolverType(grid.get(), timeStep));

        mpi::BoundaryType boundaries[6] = { 
            mpi::BoundaryType::Periodic, // +X
            mpi::BoundaryType::Periodic, // -X
            mpi::BoundaryType::Periodic, // +Y
            mpi::BoundaryType::Periodic, // -Y
            mpi::BoundaryType::Periodic, // +Z
            mpi::BoundaryType::Periodic, // -Z
        };
        using BoundaryManager = mpi::FieldBoundaryManager<FieldSolverType, GridType, ReflectBoundaryConditionMonoDirectionFdtd>;
        BoundaryManager::setupBoundaryConditions(fieldSolver, boundaries, topology, fieldExchanger, mpi_rank);

        auto setup_end = std::chrono::high_resolution_clock::now();
        auto setup_duration = std::chrono::duration<double, std::milli>(setup_end - setup_start);
        if (mpi_rank == 0)
            std::cout << " Setup Time: " << setup_duration.count() << " ms" << std::endl;

        // Grid initialization
        auto grid_start = std::chrono::high_resolution_clock::now();

        initializeGrid();

        auto grid_end = std::chrono::high_resolution_clock::now();
        auto grid_duration = std::chrono::duration<double, std::milli>(grid_end - grid_start);
        if (mpi_rank == 0)
            std::cout << " Grid Initialization Time: " << grid_duration.count() << " ms" << std::endl;

        // Running
        auto run_start = std::chrono::high_resolution_clock::now();

        for (int step = 0; step < this->numSteps; ++step)
        {
            this->fieldSolver->updateFields();
        }

        auto run_end = std::chrono::high_resolution_clock::now();
        auto run_duration = std::chrono::duration<double, std::milli>(run_end - run_start);
        if (mpi_rank == 0)
            std::cout << " Run Time: " << run_duration.count() << " ms" << std::endl;

        // Verification
        #ifdef VERIFY
            FP finalT = this->fieldSolver->dt * this->numSteps;

            Int3 begin = this->fieldSolver->internalIndexBegin;
            Int3 end = this->fieldSolver->internalIndexEnd;

            bool verified = true;

            for (int i = begin.x; i < end.x; ++i)
                for (int j = begin.y; j < end.y; ++j)
                    for (int k = begin.z; k < end.z; ++k)
                    {
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
                            
                        if (abs(expectedE.norm() - actualE.norm()) > this->maxError)
                        {
                            verified = false;
                            break;
                        }
                    }

            for (int i = begin.x; i < end.x; ++i)
                for (int j = begin.y; j < end.y; ++j)
                    for (int k = begin.z; k < end.z; ++k)
                    {
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

                        if (abs(expectedB.norm() - actualB.norm()) > this->maxError)
                        {
                            verified = false;
                            break;
                        }
                    }

            std::cout << "RANK " << mpi_rank << " Verification Result: " << verified << std::endl;
        #endif

        auto total_end = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration<double, std::milli>(total_end - total_start);
        
        MPI_Barrier(topology->getTopologyCommunicator());

        if (mpi_rank == 0)
            std::cout << "Total Time: " << total_duration.count() << " ms" << std::endl;
    }
};




int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
 
    //Int3 sections = Int3(1, 1, 1);

    // Line over X topology
    //sections = Int3(mpi_size, 1, 1);

    // Volume fill topology
    int dims[3] = {0, 0, 0};
    MPI_Dims_create(mpi_size, 3, dims);
    Int3 sections = Int3(dims[0], dims[1], dims[2]);

    // Running the test
    PerfTest perfTest(sections);

    MPI_Finalize();

    return 0;
}