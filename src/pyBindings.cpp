#include <pybind11/pybind11.h>
#include <pybind11/cast.h>
#include <pybind11/detail/common.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/pytypes.h>
#include <stdexcept>
#include <unordered_map>
#include <string>
#include <complex>
#include "MSSM.h"

namespace py = pybind11;

PYBIND11_MODULE(match_to_py, m, py::mod_gil_not_used()) {
    auto task = py::class_<Task>(m, "Task");

    py::class_<MSSM>(m, "MSSM")
        .def(py::init<std::unordered_map<std::string, std::complex<double> >, double, bool>())
        .def("updateParams", &MSSM::updateParams, py::arg("param_dict"))
        .def("getScale", &MSSM::getScale)
        .def("setScale", &MSSM::setScale, py::arg("scale"))
        .def("loopContributions", &MSSM::loopContributions, py::arg("loop"))
        .def("getParams", &MSSM::getParams)

        // batch evaluator
        .def("batch_eval", &MSSM::batch_eval, py::arg("tasks"), py::call_guard<py::gil_scoped_release>())

        // wrapper lambdas
        .def("wrap_0f", [](MSSM& self, std::string name, std::string method_name) {
            static const std::unordered_map<std::string, std::complex<double> (MSSM::*)()> map_0f = {
                {"cG", &MSSM::cG},
                {"cGt", &MSSM::cGt},
                {"cH", &MSSM::cH},
                {"cHB", &MSSM::cHB},
                {"cHBox", &MSSM::cHBox},
                {"cHBt", &MSSM::cHBt},
                {"cHD", &MSSM::cHD},
                {"cHG", &MSSM::cHG},
                {"cHGt", &MSSM::cHGt},
                {"cHW", &MSSM::cHW},
                {"cHWB", &MSSM::cHWB},
                {"cHWt", &MSSM::cHWt},
                {"cHWtB", &MSSM::cHWtB},
                {"cW", &MSSM::cW},
                {"cWt", &MSSM::cWt},
            };
            auto it = map_0f.find(method_name);
            if (it == map_0f.end()) throw std::runtime_error("Method not found: " + method_name);

            auto fn_ptr = it->second;
            return Task{name, [&self, fn_ptr](){ return (self.*fn_ptr)(); }};
        }, py::arg("name"), py::arg("method_name"))

        .def("wrap_2f", [](MSSM& self, std::string name, std::string method_name, int i1, int i2) {
            static const std::unordered_map<std::string, std::complex<double> (MSSM::*)(int, int)> map_2f = {
                {"cdB", &MSSM::cdB},
                {"cdG", &MSSM::cdG},
                {"cdH", &MSSM::cdH},
                {"cdW", &MSSM::cdW},
                {"ceB", &MSSM::ceB},
                {"ceH", &MSSM::ceH},
                {"ceW", &MSSM::ceW},
                {"cHd", &MSSM::cHd},
                {"cHe", &MSSM::cHe},
                {"cHl1", &MSSM::cHl1},
                {"cHl3", &MSSM::cHl3},
                {"cHq1", &MSSM::cHq1},
                {"cHq3", &MSSM::cHq3},
                {"cHu", &MSSM::cHu},
                {"cHud", &MSSM::cHud},
                {"cllHH", &MSSM::cllHH},
                {"cuB", &MSSM::cuB},
                {"cuG", &MSSM::cuG},
                {"cuH", &MSSM::cuH},
                {"cuW", &MSSM::cuW},
            };
            auto it = map_2f.find(method_name);
            if (it == map_2f.end()) throw std::runtime_error("Method not found: " + method_name);

            auto fn_ptr = it->second;
            return Task{name, [&self, fn_ptr, i1, i2](){ return (self.*fn_ptr)(i1, i2); }};
        }, py::arg("name"), py::arg("method_name"), py::arg("i1"), py::arg("i2"))

        .def("wrap_4f", [](MSSM& self, std::string name, std::string method_name, int i1, int i2, int i3, int i4) {
            static const std::unordered_map<std::string, std::complex<double> (MSSM::*)(int, int, int, int)> map_4f = {
                {"cdd", &MSSM::cdd},
                {"cduq", &MSSM::cduq},
                {"cduu", &MSSM::cduu},
                {"ced", &MSSM::ced},
                {"cee", &MSSM::cee},
                {"ceu", &MSSM::ceu},
                {"cld", &MSSM::cld},
                {"cle", &MSSM::cle},
                {"cledq", &MSSM::cledq},
                {"clequ1", &MSSM::clequ1},
                {"clequ3", &MSSM::clequ3},
                {"cll", &MSSM::cll},
                {"clq1", &MSSM::clq1},
                {"clq3", &MSSM::clq3},
                {"clu", &MSSM::clu},
                {"cqd1", &MSSM::cqd1},
                {"cqd8", &MSSM::cqd8},
                {"cqe", &MSSM::cqe},
                {"cqq1", &MSSM::cqq1},
                {"cqq3", &MSSM::cqq3},
                {"cqqq", &MSSM::cqqq},
                {"cqqu", &MSSM::cqqu},
                {"cqu1", &MSSM::cqu1},
                {"cqu8", &MSSM::cqu8},
                {"cquqd1", &MSSM::cquqd1},
                {"cquqd8", &MSSM::cquqd8},
                {"cud1", &MSSM::cud1},
                {"cud8", &MSSM::cud8},
                {"cuu", &MSSM::cuu},
            };
            auto it = map_4f.find(method_name);
            if (it == map_4f.end()) throw std::runtime_error("Method not found: " + method_name);

            auto fn_ptr = it->second;
            return Task{name, [&self, fn_ptr, i1, i2, i3, i4](){ return (self.*fn_ptr)(i1, i2, i3, i4); }};
        }, py::arg("name"), py::arg("method_name"), py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"))

        // Wilson coefficient methods
        .def("cG", &MSSM::cG, py::call_guard<py::gil_scoped_release>())
        .def("cGt", &MSSM::cGt, py::call_guard<py::gil_scoped_release>())
        .def("cH", &MSSM::cH, py::call_guard<py::gil_scoped_release>())
        .def("cHB", &MSSM::cHB, py::call_guard<py::gil_scoped_release>())
        .def("cHBox", &MSSM::cHBox, py::call_guard<py::gil_scoped_release>())
        .def("cHBt", &MSSM::cHBt, py::call_guard<py::gil_scoped_release>())
        .def("cHD", &MSSM::cHD, py::call_guard<py::gil_scoped_release>())
        .def("cHG", &MSSM::cHG, py::call_guard<py::gil_scoped_release>())
        .def("cHGt", &MSSM::cHGt, py::call_guard<py::gil_scoped_release>())
        .def("cHW", &MSSM::cHW, py::call_guard<py::gil_scoped_release>())
        .def("cHWB", &MSSM::cHWB, py::call_guard<py::gil_scoped_release>())
        .def("cHWt", &MSSM::cHWt, py::call_guard<py::gil_scoped_release>())
        .def("cHWtB", &MSSM::cHWtB, py::call_guard<py::gil_scoped_release>())
        .def("cW", &MSSM::cW, py::call_guard<py::gil_scoped_release>())
        .def("cWt", &MSSM::cWt, py::call_guard<py::gil_scoped_release>())
        .def("cdB", &MSSM::cdB, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cdG", &MSSM::cdG, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cdH", &MSSM::cdH, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cdW", &MSSM::cdW, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("ceB", &MSSM::ceB, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("ceH", &MSSM::ceH, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("ceW", &MSSM::ceW, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cHd", &MSSM::cHd, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cHe", &MSSM::cHe, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cHl1", &MSSM::cHl1, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cHl3", &MSSM::cHl3, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cHq1", &MSSM::cHq1, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cHq3", &MSSM::cHq3, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cHu", &MSSM::cHu, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cHud", &MSSM::cHud, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cllHH", &MSSM::cllHH, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cuB", &MSSM::cuB, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cuG", &MSSM::cuG, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cuH", &MSSM::cuH, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cuW", &MSSM::cuW, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cdd", &MSSM::cdd, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cduq", &MSSM::cduq, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cduu", &MSSM::cduu, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("ced", &MSSM::ced, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cee", &MSSM::cee, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("ceu", &MSSM::ceu, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cld", &MSSM::cld, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cle", &MSSM::cle, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cledq", &MSSM::cledq, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("clequ1", &MSSM::clequ1, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("clequ3", &MSSM::clequ3, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cll", &MSSM::cll, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("clq1", &MSSM::clq1, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("clq3", &MSSM::clq3, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("clu", &MSSM::clu, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cqd1", &MSSM::cqd1, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cqd8", &MSSM::cqd8, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cqe", &MSSM::cqe, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cqq1", &MSSM::cqq1, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cqq3", &MSSM::cqq3, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cqqq", &MSSM::cqqq, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cqqu", &MSSM::cqqu, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cqu1", &MSSM::cqu1, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cqu8", &MSSM::cqu8, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cquqd1", &MSSM::cquqd1, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cquqd8", &MSSM::cquqd8, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cud1", &MSSM::cud1, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cud8", &MSSM::cud8, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cuu", &MSSM::cuu, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>());
}
