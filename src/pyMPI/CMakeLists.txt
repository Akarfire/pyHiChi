option(PYHICHI_WITH_MPI "Build MPI bindings into pyHiChi" ON)

if(NOT PYHICHI_WITH_MPI)
    message(STATUS "pyHiChi: MPI bindings disabled (PYHICHI_WITH_MPI=OFF)")
    return()
endif()

if(POLICY CMP0079)
    cmake_policy(SET CMP0079 NEW)
endif()

if(NOT PYTHON_EXECUTABLE)
    find_package(PythonInterp QUIET)
endif()

if(NOT PYTHON_EXECUTABLE)
    message(FATAL_ERROR
        "PYTHON_EXECUTABLE is not set. Either pass /python <path> to "
        "build_windows_custom.bat, or configure with "
        "-DPYTHON_EXECUTABLE:FILEPATH=/path/to/python.")
endif()

message(STATUS "pyHiChi-MPI: using Python interpreter at ${PYTHON_EXECUTABLE}")


execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c
        "import mpi4py, sys; sys.stdout.write(mpi4py.get_include())"
    OUTPUT_VARIABLE MPI4PY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE MPI4PY_QUERY_RESULT
)

if(NOT MPI4PY_QUERY_RESULT EQUAL 0 OR NOT MPI4PY_INCLUDE_DIR)
    message(FATAL_ERROR
        "Could not locate mpi4py headers using Python at\n"
        "  ${PYTHON_EXECUTABLE}\n"
        "Is mpi4py installed in that interpreter? Verify with:\n"
        "  \"${PYTHON_EXECUTABLE}\" -c \"import mpi4py; print(mpi4py.get_include())\"")
endif()

message(STATUS "pyHiChi-MPI: mpi4py include dir = ${MPI4PY_INCLUDE_DIR}")


find_package(MPI REQUIRED COMPONENTS CXX)


target_sources(pyHiChi PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/pyMPI.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/../MPI/src/MPI_Utilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../MPI/src/MPI_Topology.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../MPI/src/MPI_FieldUtils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../MPI/src/MPI_FieldExchanger.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../MPI/src/MPI_GridSlicer.cpp
)

target_include_directories(pyHiChi PRIVATE
    ../MPI/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${MPI_CXX_INCLUDE_DIRS}
    ${MPI4PY_INCLUDE_DIR}
)

target_link_libraries(pyHiChi PUBLIC MPI::MPI_CXX)

target_compile_definitions(pyHiChi PRIVATE PYHICHI_WITH_MPI=1)