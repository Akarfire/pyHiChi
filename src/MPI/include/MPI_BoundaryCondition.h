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
        SideEnum direction_;
        std::shared_ptr<mpi::FieldExchanger> exchanger;
        std::shared_ptr<mpi::Topology> topology;
        int rank;

        mpi::Direction mpiDirection;

    void checkCommunicator()
    {
        if (!topology->isValidOnThisRank())
            throw std::runtime_error("Attempted to create MPI_FieldBoundaryCondition on rank, that does not participate in topology!");
    }

    public:

        MPI_FieldBoundaryCondition(GridType* grid,
            Int3 leftBorderIndex, Int3 rightBorderIndex, CoordinateEnum axis, SideEnum direction,
            std::shared_ptr<mpi::FieldExchanger> fieldExchanger,
            std::shared_ptr<mpi::Topology> topology, int mpi_rank) :

            FieldBoundaryCondition(grid, leftBorderIndex, rightBorderIndex, axis), direction_(direction),
            exchanger(fieldExchanger), topology(topology), rank(mpi_rank)
        { 
            mpiDirection = static_cast<mpi::Direction>((direction_ == SideEnum::RIGHT) ? static_cast<int>(axis) * 2 : static_cast<int>(axis) * 2 + 1);
            checkCommunicator(); 
        }

        // constructor for loading
        explicit MPI_FieldBoundaryCondition(GridType* grid,
            Int3 leftBorderIndex, Int3 rightBorderIndex, 
            std::shared_ptr<mpi::FieldExchanger> fieldExchanger,
            std::shared_ptr<mpi::Topology> topology, int mpi_rank) :

            FieldBoundaryCondition(grid, leftBorderIndex, rightBorderIndex), 
            exchanger(fieldExchanger), topology(topology), rank(mpi_rank)
        {
            checkCommunicator();
        }

        void generateB(FP time) override;
        void generateE(FP time) override;

        virtual void save(std::ostream& ostr);
        virtual void load(std::istream& istr);

        FieldBoundaryCondition<GridType>* createInstance(
            GridType* grid, Int3 leftBorderIndex, Int3 rightBorderIndex, CoordinateEnum axis) override {
            return new MPI_FieldBoundaryCondition(grid, leftBorderIndex, rightBorderIndex, axis, direction_, exchanger, topology, rank);
        }
    };

    template <class GridType>
    inline void MPI_FieldBoundaryCondition<GridType>::generateE(FP time)
    {
        exchanger->performExchangeInDirection(grid->Ex.getData(), mpiDirection, topology, rank, topology->getTopologyCommunicator());
        exchanger->performExchangeInDirection(grid->Ey.getData(), mpiDirection, topology, rank, topology->getTopologyCommunicator());
        exchanger->performExchangeInDirection(grid->Ez.getData(), mpiDirection, topology, rank, topology->getTopologyCommunicator());
    }

    template <class GridType>
    inline void MPI_FieldBoundaryCondition<GridType>::generateB(FP time)
    {
        exchanger->performExchangeInDirection(grid->Bx.getData(), mpiDirection, topology, rank, topology->getTopologyCommunicator());
        exchanger->performExchangeInDirection(grid->By.getData(), mpiDirection, topology, rank, topology->getTopologyCommunicator());
        exchanger->performExchangeInDirection(grid->Bz.getData(), mpiDirection, topology, rank, topology->getTopologyCommunicator());
    }

    template <class GridType>
    inline void MPI_FieldBoundaryCondition<GridType>::save(std::ostream& ostr)
    {
        FieldBoundaryCondition<YeeGrid>::save(ostr);

        ostr.write((char*)&direction_, sizeof(direction_));
    }

    template <class GridType>
    inline void MPI_FieldBoundaryCondition<GridType>::load(std::istream& istr)
    {
        FieldBoundaryCondition<YeeGrid>::load(istr);

        istr.read((char*)&direction_, sizeof(direction_));
        mpiDirection = static_cast<mpi::Direction>((direction_ == SideEnum::RIGHT) ? static_cast<int>(axis) * 2 : static_cast<int>(axis) * 2 + 1);
    }
}