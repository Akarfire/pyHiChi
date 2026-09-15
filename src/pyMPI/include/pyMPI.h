#pragma once
#include "pybind11/pybind11.h"

// Initializes python bindings for MPI related primitives
void initMpiBinding(pybind11::module& module);