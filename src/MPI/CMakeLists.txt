include_directories(
    include
    ../Core/include
    ../Modules/Particle/include
    ../Modules/Field/include
    ${FFT_INCLUDES})

add_executable(mpi
    ./src/Main.cpp)

set_target_properties(mpi PROPERTIES FOLDER MPI)