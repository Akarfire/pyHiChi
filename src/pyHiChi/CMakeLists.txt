cmake_minimum_required(VERSION 3.1)

project(pyHiChi)

include_directories(include
    ../Core/include
    ../Modules/Particle/include
    ../Modules/Field/include
    ../Modules/Field/include/Mapping
    ../Modules/Field/include/FieldConfiguration
    ${FFT_INCLUDES})

if (USE_OMP)
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()

set(pyHiChi_source
    include/pyField.h
    include/pyFieldInterface.h
    include/pyFieldMacroses.h
    src/pyHiChi.cpp)

pybind11_add_module(pyHiChi ${pyHiChi_source})

if (APPLE)
    target_link_libraries(pyHiChi PUBLIC ${FFT_LIBS})
elseif (UNIX)
    target_link_libraries(pyHiChi PUBLIC -Wl,--start-group ${FFT_LIBS} -Wl,--end-group)
else()
    target_link_libraries(pyHiChi PUBLIC ${FFT_LIBS})
endif()

# MPI integration
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../pyMPI
                 ${CMAKE_CURRENT_BINARY_DIR}/pyMPI)

# Stub (.pyi) file generation

# Copying fftw if needed
if(WIN32 AND FFTW_DIR)
    file(GLOB FFTW_RUNTIME_DLLS
        "${FFTW_DIR}/bin/*.dll"
        "${FFTW_DIR}/*.dll"
        "${FFTW_DIR}/lib/*.dll")
    if(FFTW_RUNTIME_DLLS)
        add_custom_command(TARGET pyHiChi POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    ${FFTW_RUNTIME_DLLS}
                    "$<TARGET_FILE_DIR:pyHiChi>"
            COMMENT "Copying FFTW runtime DLLs next to pyHiChi.pyd"
        )
    else()
        message(WARNING
            "FFTW_DIR is set (${FFTW_DIR}) but no .dll was found there. "
            "pyHiChi may fail to import; add FFTW's bin/ to PATH as a fallback.")
    endif()
endif()


if(DEFINED PYTHON_EXECUTABLE AND NOT PYTHON_EXECUTABLE STREQUAL "" AND
   NOT PYTHON_EXECUTABLE STREQUAL "python")
    set(PYHICHI_STUBGEN_PYTHON "${PYTHON_EXECUTABLE}")
elseif(WIN32 AND EXISTS "${CMAKE_SOURCE_DIR}/../.venv/Scripts/python.exe")
    set(PYHICHI_STUBGEN_PYTHON "${CMAKE_SOURCE_DIR}/../.venv/Scripts/python.exe")
elseif(EXISTS "${CMAKE_SOURCE_DIR}/../.venv/bin/python")
    set(PYHICHI_STUBGEN_PYTHON "${CMAKE_SOURCE_DIR}/../.venv/bin/python")
else()
    set(PYHICHI_STUBGEN_PYTHON "python")
endif()

message(STATUS "pyHiChi: stubgen will use Python at ${PYHICHI_STUBGEN_PYTHON}")

if(WIN32)
    set(PYHICHI_PYD_DIR "${CMAKE_BINARY_DIR}/src/pyHiChi/Release")
else()
    set(PYHICHI_PYD_DIR "${CMAKE_BINARY_DIR}/src/pyHiChi")
endif()

add_custom_command(TARGET pyHiChi POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E env
            "PYTHONPATH=${PYHICHI_PYD_DIR}"
            "${PYHICHI_STUBGEN_PYTHON}" -m pybind11_stubgen pyHiChi
            -o "${CMAKE_BINARY_DIR}/src/pyHiChi/Release"
            --ignore-all-errors
    WORKING_DIRECTORY "${PYHICHI_PYD_DIR}"
    COMMENT "Generating .pyi stub for pyHiChi"
    VERBATIM
)