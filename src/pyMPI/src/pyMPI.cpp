#include "pyMPI.h"

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include <pybind11/operators.h>

#include <mpi4py/mpi4py.h>

#include "HiChi_MPI.h"

namespace py = pybind11;

// Helper functions

// In python invalid ranks are returned as `None`
py::object pythonizeRank(int rank)
{
    return (rank == mpi::MPI_INVALID_RANK) ? py::none() : py::cast(rank);
}

// In python invalid sections are `None` type
py::object pythonizeSection(const pfc::Int3& section)
{
    if (section == mpi::MPI_INVALID_SECTION)
        return py::none();
    return py::cast(section);
}

// Checks if mpi is initialized and throws an error if not
void checkMpiInitialization()
{
    int inited = 0;
    MPI_Initialized(&inited);
    if (!inited) {
        throw std::runtime_error(
            "MPI PYTHON BINDING : MPI is not initialized!");
    }
}



// Initializes python bindings for MPI related primitives
void initMpiBinding(py::module& module)
{
    // Int3 binding

    py::class_<pfc::Int3>(module, "IntVector3d")
        .def(py::init<>())
        .def(py::init<int, int, int>(), py::arg("x"), py::arg("y"), py::arg("z"))
        .def("__str__", &pfc::Int3::toString)

        .def(py::self + py::self)
        .def(py::self += py::self)
        .def(py::self - py::self)
        .def(py::self -= py::self)
        .def(-py::self)
        .def(int() * py::self)
        .def(py::self * int())
        .def(py::self *= int())
        .def(py::self * py::self)
        .def(py::self *= py::self)
        .def(py::self / py::self)
        .def(py::self /= py::self)
        .def(py::self / int())
        .def(py::self /= int())

        .def_readwrite("x", &pfc::Int3::x)
        .def_readwrite("y", &pfc::Int3::y)
        .def_readwrite("z", &pfc::Int3::z)
        ;


    // mpi4py setup
    if (import_mpi4py() < 0)
        throw py::error_already_set();

    // Loop Type
    py::enum_<mpi::Topology::LoopType>(module, "MPI_LoopType")
        .value("None",    mpi::Topology::LoopType::None)
        .value("LoopX",   mpi::Topology::LoopType::LoopX)
        .value("LoopY",   mpi::Topology::LoopType::LoopY)
        .value("LoopZ",   mpi::Topology::LoopType::LoopZ)
        .value("LoopXY",  mpi::Topology::LoopType::LoopXY)
        .value("LoopYZ",  mpi::Topology::LoopType::LoopYZ)
        .value("LoopXZ",  mpi::Topology::LoopType::LoopXZ)
        .value("LoopXYZ", mpi::Topology::LoopType::LoopXYZ)
        .export_values();

    // Direction
    py::enum_<mpi::Direction>(module, "MPI_Direction")
        .value("pX", mpi::Direction::pX)
        .value("nX", mpi::Direction::nX)
        .value("pY", mpi::Direction::pY)
        .value("nY", mpi::Direction::nY)
        .value("pZ", mpi::Direction::pZ)
        .value("nZ", mpi::Direction::nZ)
        .export_values();

    module.def("mpi_invert_direction", &mpi::invertDirection,
        pybind11::arg("direction"),
        "Returns the opposite of the given direction.");


    // Topology

    py::class_<mpi::Topology>(module, "Topology")

        // Constructor from Loop Type
        .def(py::init([](pfc::Int3 sections,
                mpi::Topology::LoopType loop_type,
                int node_count) 
        {
            checkMpiInitialization();
            return new mpi::Topology(sections, loop_type, node_count);
        }),
        py::arg("sections"), py::arg("loop_type"), py::arg("node_count"))

        // Constructor from Loop Mask
        .def(py::init([](pfc::Int3 sections,
                        py::tuple loop_mask,
                        int node_count) {
                checkMpiInitialization();
                if (loop_mask.size() != 3)
                    throw py::value_error("loop_mask must be a length-3 tuple");
                bool mask[3] = {
                    loop_mask[0].cast<bool>(),
                    loop_mask[1].cast<bool>(),
                    loop_mask[2].cast<bool>()
                };
                return new mpi::Topology(sections, mask, node_count);
            }),
            py::arg("sections"), py::arg("loop_mask"), py::arg("node_count"))

        .def_property_readonly("sections",
            [](const mpi::Topology& topology) { return topology.getSections(); })

        .def("get_responsible_node",
            [](const mpi::Topology& topology, const pfc::Int3& section) 
            {
                return pythonizeRank(topology.getResponsibleNode(section));
            },
            py::arg("section"))

        .def("get_node_section",
            [](const mpi::Topology& topology, int rank) 
            {
                return pythonizeSection(topology.getNodeSection(rank));
            },
            py::arg("rank"))

        .def("get_neighbor",
            [](const mpi::Topology& topology, int rank, const pfc::Int3& offset) 
            {
                return pythonizeRank(topology.getNeighbor(rank, offset));
            },
            py::arg("node_rank"), py::arg("offset"));
}