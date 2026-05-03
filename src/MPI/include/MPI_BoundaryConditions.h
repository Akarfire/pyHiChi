#include "FieldBoundaryCondition.h"
#include "MPI_FieldExchanger.h"
#include "MPI_Topology.h"

namespace pfc
{
    class MPI_FieldBoundaryConditions : public FieldBoundaryCondition<YeeGrid>
    {
    protected:
        std::shared_ptr<mpi::FieldExchanger> exchanger;
        std::shared_ptr<mpi::Topology> topology;
        int rank;

    public:

        MPI_FieldBoundaryConditions(YeeGrid* grid,
            Int3 leftBorderIndex, Int3 rightBorderIndex, CoordinateEnum axis, 
            std::shared_ptr<mpi::FieldExchanger> fieldExchanger,
            std::shared_ptr<mpi::Topology> topology, int mpi_rank) :

            FieldBoundaryCondition(grid, leftBorderIndex, rightBorderIndex, axis), 
            exchanger(fieldExchanger), topology(topology), rank(mpi_rank)
        {}

        // constructor for loading
        explicit MPI_FieldBoundaryConditions(YeeGrid* grid,
            Int3 leftBorderIndex, Int3 rightBorderIndex, 
            std::shared_ptr<mpi::FieldExchanger> fieldExchanger,
            std::shared_ptr<mpi::Topology> topology, int mpi_rank) :

            FieldBoundaryCondition(grid, leftBorderIndex, rightBorderIndex), 
            exchanger(fieldExchanger), topology(topology), rank(mpi_rank)
        {}

        void generateB(FP time) override;
        void generateE(FP time) override;

        FieldBoundaryCondition<YeeGrid>* createInstance(
            YeeGrid* grid, Int3 leftBorderIndex, Int3 rightBorderIndex, CoordinateEnum axis) override {
            return new MPI_FieldBoundaryConditions(grid, leftBorderIndex, rightBorderIndex, axis, exchanger, topology, rank);
        }
    };

    inline void MPI_FieldBoundaryConditions::generateE(FP time)
    {
        exchanger->performExchangeOverAxis(grid->Ex.getData(), axis, *topology, rank, MPI_COMM_WORLD);
        exchanger->performExchangeOverAxis(grid->Ey.getData(), axis, *topology, rank, MPI_COMM_WORLD);
        exchanger->performExchangeOverAxis(grid->Ez.getData(), axis, *topology, rank, MPI_COMM_WORLD);
    }

    inline void MPI_FieldBoundaryConditions::generateB(FP time)
    {
        exchanger->performExchangeOverAxis(grid->Bx.getData(), axis, *topology, rank, MPI_COMM_WORLD);
        exchanger->performExchangeOverAxis(grid->By.getData(), axis, *topology, rank, MPI_COMM_WORLD);
        exchanger->performExchangeOverAxis(grid->Bz.getData(), axis, *topology, rank, MPI_COMM_WORLD);
    }
}