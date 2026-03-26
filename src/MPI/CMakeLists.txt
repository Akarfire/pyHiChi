include_directories(
    include
    ../Core/include
    ../Modules/Particle/include
    ../Modules/Field/include
    ${FFT_INCLUDES})

add_executable(mpi
    src/MPI_Utilities.cpp
    src/Main.cpp)

target_link_libraries(mpi
    MPI::MPI_CXX)

set_target_properties(mpi PROPERTIES FOLDER MPI)