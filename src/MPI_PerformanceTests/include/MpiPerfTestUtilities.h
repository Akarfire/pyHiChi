#pragma once

#include <cstring>
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include "mpi.h"

#include "Vectors.h"

struct ProgramOptions
{
    int gridSizeLongitudinal = 32;
    int gridSizeTransverse   = 32;

    enum class SplitMode { Uniform, Line };
    SplitMode splitMode = SplitMode::Uniform;

    static void printUsage(const char* progName)
    {
        std::cout <<
            "Usage: " << progName << " [options]\n"
            "Options:\n"
            "  --grid-longitudinal N   Grid size along x-axis (longitudinal) [default: 32]\n"
            "  --grid-transverse   N   Grid size along y,z axes (transverse)    [default: 32]\n"
            "  --split-mode MODE       Section split mode: 'uniform' or 'x'     [default: uniform]\n"
            "  --help, -h              Show this help message\n";
    }

    static bool parse(int argc, char** argv, ProgramOptions& opts)
    {
        for (int i = 1; i < argc; ++i)
        {
            std::string arg = argv[i];

            if (arg == "--help" || arg == "-h")
            {
                printUsage(argv[0]);
                return false;
            }
            else if (arg == "--grid-longitudinal" && i + 1 < argc)
            {
                opts.gridSizeLongitudinal = std::atoi(argv[++i]);
            }
            else if (arg == "--grid-transverse" && i + 1 < argc)
            {
                opts.gridSizeTransverse = std::atoi(argv[++i]);
            }
            else if (arg == "--split-mode" && i + 1 < argc)
            {
                std::string mode = argv[++i];
                if (mode == "uniform")
                    opts.splitMode = SplitMode::Uniform;
                else if (mode == "line")
                    opts.splitMode = SplitMode::Line;
                else
                {
                    if (mpi_rank_is_zero_safe())
                        std::cerr << "Unknown split mode: " << mode << "\n";
                    return false;
                }
            }
            else
            {
                if (mpi_rank_is_zero_safe())
                    std::cerr << "Unknown argument: " << arg << "\n";
                return false;
            }
        }

        if (opts.gridSizeLongitudinal <= 0 || opts.gridSizeTransverse <= 0)
        {
            if (mpi_rank_is_zero_safe())
                std::cerr << "Grid sizes must be positive\n";
            return false;
        }
        return true;
    }

private:
    static bool mpi_rank_is_zero_safe()
    {
        int rank = 0, initialized = 0;
        MPI_Initialized(&initialized);
        if (initialized) MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        return rank == 0;
    }
};


// Opens "perf_log_L<long>_T<trans>_<mode>_sec<sx>x<sy>x<sz>_rank<r>.log"
inline std::ofstream openPerfLog(const ProgramOptions& opts,
                                 const pfc::Int3& sections,
                                 int mpiRank)
{
    std::ostringstream name;
    name << "perf_log"
         << "_L" << opts.gridSizeLongitudinal
         << "_T" << opts.gridSizeTransverse
         << "_"  << (opts.splitMode == ProgramOptions::SplitMode::Line ? "line" : "uniform")
         << "_sec" << sections.x << "x" << sections.y << "x" << sections.z
         << "_rank" << mpiRank
         << ".log";
    return std::ofstream(name.str());
}