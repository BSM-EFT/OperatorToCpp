#include <pybind11/pybind11.h>
#include <pybind11/cast.h>
#include <pybind11/native_enum.h>
#include <pybind11/detail/common.h>
#include <pybind11/stl.h>
#include <unordered_map>
#include <string>
#include "MSSM.h"

namespace py = pybind11;

PYBIND11_MODULE(matching, m, py::mod_gil_not_used()) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    py::class_<MSSM>(m, "mssm")
        .def(py::init<>())
        .def(py::init<std::unordered_map<std::string, double>>())
        .def("cH", &MSSM::cH, py::arg("mubarsq"), py::arg("hbar"));
}
