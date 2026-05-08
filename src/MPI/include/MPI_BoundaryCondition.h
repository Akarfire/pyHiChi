#pragma once

#include "FieldBoundaryCondition.h"
#include "MPI_FieldExchanger.h"
#include "MPI_Topology.h"

namespace pfc
{
    template <class GridType>
    class MPI_FieldBoundaryCondition : public FieldBoundaryCondition<GridType>
    {
    protected:
        std::shared_ptr<mpi::FieldExchanger> exchanger;
        std::shared_ptr<mpi::Topology> topology;
        int rank;

    public:

        MPI_FieldBoundaryCondition(GridType* grid,
            Int3 leftBorderIndex, Int3 rightBorderIndex, CoordinateEnum axis, 
            std::shared_ptr<mpi::FieldExchanger> fieldExchanger,
            std::shared_ptr<mpi::Topology> topology, int mpi_rank) :

            FieldBoundaryCondition(grid, leftBorderIndex, rightBorderIndex, axis), 
            exchanger(fieldExchanger), topology(topology), rank(mpi_rank)
        {}

        // constructor for loading
        explicit MPI_FieldBoundaryCondition(GridType* grid,
            Int3 leftBorderIndex, Int3 rightBorderIndex, 
            std::shared_ptr<mpi::FieldExchanger> fieldExchanger,
            std::shared_ptr<mpi::Topology> topology, int mpi_rank) :

            FieldBoundaryCondition(grid, leftBorderIndex, rightBorderIndex), 
            exchanger(fieldExchanger), topology(topology), rank(mpi_rank)
        {}

        void generateB(FP time) override;
        void generateE(FP time) override;

        FieldBoundaryCondition<GridType>* createInstance(
            GridType* grid, Int3 leftBorderIndex, Int3 rightBorderIndex, CoordinateEnum axis) override {
            return new MPI_FieldBoundaryCondition(grid, leftBorderIndex, rightBorderIndex, axis, exchanger, topology, rank);
        }
    };

    template <class GridType>
    inline void MPI_FieldBoundaryCondition<GridType>::generateE(FP time)
    {
        exchanger->performExchangeOverAxis(grid->Ex.getData(), axis, *topology, rank, MPI_COMM_WORLD);
        exchanger->performExchangeOverAxis(grid->Ey.getData(), axis, *topology, rank, MPI_COMM_WORLD);
        exchanger->performExchangeOverAxis(grid->Ez.getData(), axis, *topology, rank, MPI_COMM_WORLD);
    }

    template <class GridType>
    inline void MPI_FieldBoundaryCondition<GridType>::generateB(FP time)
    {
        exchanger->performExchangeOverAxis(grid->Bx.getData(), axis, *topology, rank, MPI_COMM_WORLD);
        exchanger->performExchangeOverAxis(grid->By.getData(), axis, *topology, rank, MPI_COMM_WORLD);
        exchanger->performExchangeOverAxis(grid->Bz.getData(), axis, *topology, rank, MPI_COMM_WORLD);
    }
}