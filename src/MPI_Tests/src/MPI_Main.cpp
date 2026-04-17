#include "mpi.h"
#include "gtest/gtest.h"

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int result = RUN_ALL_TESTS();

    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < size; i++)
    {
        if (rank == i)
            std::cerr << "[Rank " << rank << "] Result: " << ((result == 0) ? "[PASSED]" : "[FAILED]") << ", Num Failed: " << result << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    MPI_Finalize();

    return result;
}