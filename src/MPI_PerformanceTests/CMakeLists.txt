include_directories(
    include
    ../../3rdparty/google/googletest/googletest/include
    ../Core/include
    ../Modules/Particle/include
    ../Modules/Field/include
    ../MPI/include
    ${FFT_INCLUDES})


file(GLOB_RECURSE MPI_SOURCES "../MPI/src/*.cpp")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(mpi_perf_tests
    src/testMPI_PerfFields.cpp

    ${MPI_SOURCES}
)

target_link_libraries(mpi_perf_tests
    MPI::MPI_CXX    
    ${FFT_LIBS}
)


add_executable(reference_perf_test
    src/test_PerfFields_Reference.cpp

    ${MPI_SOURCES}
)

target_link_libraries(reference_perf_test
    MPI::MPI_CXX    
    ${FFT_LIBS}
)