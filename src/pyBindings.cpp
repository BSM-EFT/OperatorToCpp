#include <pybind11/pybind11.h>
#include <pybind11/cast.h>
#include <pybind11/detail/common.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <stdexcept>
#include <unordered_map>
#include <string>
#include <complex>
#include "modelName.h"
#include "pybind11/pytypes.h"

namespace py = pybind11;

PYBIND11_MODULE(match_to_py, m, py::mod_gil_not_used()) {
    auto task = py::class_<Task>(m, "Task");

    py::class_<Model>(m, py_class)
        .def(py::init<std::unordered_map<std::string, std::complex<double> > >())
        .def("updateParams", &Model::updateParams, py::arg("param_dict"))
        .def("getScale", &Model::getScale)
        .def("setScale", &Model::setScale, py::arg("scale"))
        .def("loopContributions", &Model::loopContributions, py::arg("loop"))
        .def("getParams", &Model::getParams)

        // batch evaluator
        .def("batch_eval", &Model::batch_eval, py::arg("tasks"), py::call_guard<py::gil_scoped_release>())

        // wrapper lambdas
        .def("wrap_0f", [](Model& self, std::string name, std::string method_name) {
            static const std::unordered_map<std::string, std::complex<double> (Model::*)()> map_0f = {
                {"cG", &Model::cG},
                {"cGt", &Model::cGt},
                {"cW", &Model::cW},
                {"cWt", &Model::cWt},
                {"cH", &Model::cH},
                {"cHD", &Model::cHD},
                {"cHBox", &Model::cHBox},
                {"cHG", &Model::cHG},
                {"cHGt", &Model::cHGt},
                {"cHW", &Model::cHW},
                {"cHWt", &Model::cHWt},
                {"cHB", &Model::cHB},
                {"cHBt", &Model::cHBt},
                {"cHWB", &Model::cHWB},
                {"cHWtB", &Model::cHWtB},
            };
            auto it = map_0f.find(method_name);
            if (it == map_0f.end()) throw std::runtime_error("Method not found: " + method_name);

            auto fn_ptr = it->second;
            return Task{name, [&self, fn_ptr](){
                return (self.*fn_ptr)();
            }};
        }, py::arg("name"), py::arg("method_name"))

        .def("wrap_2f", [](Model& self, std::string name, std::string method_name, int i1, int i2) {
            static const std::unordered_map<std::string, std::complex<double> (Model::*)(int, int)> map_2f = {
                {"cllHH", &Model::cllHH},
                {"ceH", &Model::ceH},
                {"cdH", &Model::cdH},
                {"cuH", &Model::cuH},
                {"ceB", &Model::ceB},
                {"ceW", &Model::ceW},
                {"cdB", &Model::cdB},
                {"cdW", &Model::cdW},
                {"cdG", &Model::cdG},
                {"cuB", &Model::cuB},
                {"cuW", &Model::cuW},
                {"cuG", &Model::cuG},
                {"cHe", &Model::cHe},
                {"cHu", &Model::cHu},
                {"cHd", &Model::cHd},
                {"cHud", &Model::cHud},
                {"cHl1", &Model::cHl1},
                {"cHl3", &Model::cHl3},
                {"cHq1", &Model::cHq1},
                {"cHq3", &Model::cHq3}
            };
            auto it = map_2f.find(method_name);
            if (it == map_2f.end()) throw std::runtime_error("Method not found: " + method_name);

            auto fn_ptr = it->second;
            return Task{name, [&self, fn_ptr, i1, i2](){
                return (self.*fn_ptr)(i1, i2);
            }};
        }, py::arg("name"), py::arg("method_name"), py::arg("i1"), py::arg("i2"))

        .def("wrap_4f", [](Model& self, std::string name, std::string method_name, int i1, int i2, int i3, int i4) {
            static const std::unordered_map<std::string, std::complex<double> (Model::*)(int, int, int, int)> map_4f = {
                {"cll", &Model::cll},
                {"cqq1", &Model::cqq1},
                {"cqq3", &Model::cqq3},
                {"clq1", &Model::clq1},
                {"clq3", &Model::clq3},
                {"cee", &Model::cee},
                {"cuu", &Model::cuu},
                {"cdd", &Model::cdd},
                {"ceu", &Model::ceu},
                {"ced", &Model::ced},
                {"cle", &Model::cle},
                {"clu", &Model::clu},
                {"cld", &Model::cld},
                {"cqe", &Model::cqe},
                {"cqu1", &Model::cqu1},
                {"cqu8", &Model::cqu8},
                {"cqd1", &Model::cqd1},
                {"cqd8", &Model::cqd8},
                {"cud1", &Model::cud1},
                {"cud8", &Model::cud8},
                {"cquqd1", &Model::cquqd1},
                {"cquqd8", &Model::cquqd8},
                {"clequ1", &Model::clequ1},
                {"clequ3", &Model::clequ3},
                {"cledq", &Model::cledq},
                {"cqqq", &Model::cqqq},
                {"cduq", &Model::cduq},
                {"cduu", &Model::cduu},
                {"cqqu", &Model::cqqu},
            };
            auto it = map_4f.find(method_name);
            if (it == map_4f.end()) throw std::runtime_error("Method not found: " + method_name);

            auto fn_ptr = it->second;
            return Task{name, [&self, fn_ptr, i1, i2, i3, i4](){
                return (self.*fn_ptr)(i1, i2, i3, i4);
            }};
        }, py::arg("name"), py::arg("method_name"), py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"))

        // Wilson coefficient methods
        .def("cllHH", &Model::cllHH, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cG", &Model::cG, py::call_guard<py::gil_scoped_release>())
        .def("cW", &Model::cW, py::call_guard<py::gil_scoped_release>())
        .def("cGt", &Model::cGt, py::call_guard<py::gil_scoped_release>())
        .def("cWt", &Model::cWt, py::call_guard<py::gil_scoped_release>())
        .def("cH", &Model::cH, py::call_guard<py::gil_scoped_release>())
        .def("cHBox", &Model::cHBox, py::call_guard<py::gil_scoped_release>())
        .def("cHD", &Model::cHD, py::call_guard<py::gil_scoped_release>())
        .def("cHG", &Model::cHG, py::call_guard<py::gil_scoped_release>())
        .def("cHW", &Model::cHW, py::call_guard<py::gil_scoped_release>())
        .def("cHB", &Model::cHB, py::call_guard<py::gil_scoped_release>())
        .def("cHWB", &Model::cHWB, py::call_guard<py::gil_scoped_release>())
        .def("cHGt", &Model::cHGt, py::call_guard<py::gil_scoped_release>())
        .def("cHWt", &Model::cHWt, py::call_guard<py::gil_scoped_release>())
        .def("cHBt", &Model::cHBt, py::call_guard<py::gil_scoped_release>())
        .def("cHWtB", &Model::cHWtB, py::call_guard<py::gil_scoped_release>())
        .def("ceH", &Model::ceH, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cuH", &Model::cuH, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cdH", &Model::cdH, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("ceW", &Model::ceW, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("ceB", &Model::ceB, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cuG", &Model::cuG, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cuW", &Model::cuW, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cuB", &Model::cuB, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cdG", &Model::cdG, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cdW", &Model::cdW, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cdB", &Model::cdB, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cHl1", &Model::cHl1, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cHl3", &Model::cHl3, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cHe", &Model::cHe, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cHq1", &Model::cHq1, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cHq3", &Model::cHq3, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cHu", &Model::cHu, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cHd", &Model::cHd, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cHud", &Model::cHud, py::arg("i1"), py::arg("i2"), py::call_guard<py::gil_scoped_release>())
        .def("cll", &Model::cll, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cqq1", &Model::cqq1, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cqq3", &Model::cqq3, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("clq1", &Model::clq1, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("clq3", &Model::clq3, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cee", &Model::cee, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cuu", &Model::cuu, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cdd", &Model::cdd, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("ceu", &Model::ceu, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("ced", &Model::ced, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cud1", &Model::cud1, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cud8", &Model::cud8, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cle", &Model::cle, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("clu", &Model::clu, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cld", &Model::cld, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cqe", &Model::cqe, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cqu1", &Model::cqu1, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cqu8", &Model::cqu8, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cqd1", &Model::cqd1, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cqd8", &Model::cqd8, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cledq", &Model::cledq, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cquqd1", &Model::cquqd1, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cquqd8", &Model::cquqd8, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("clequ1", &Model::clequ1, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("clequ3", &Model::clequ3, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cduq", &Model::cduq, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cqqu", &Model::cqqu, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cqqq", &Model::cqqq, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>())
        .def("cduu", &Model::cduu, py::arg("i1"), py::arg("i2"), py::arg("i3"), py::arg("i4"), py::call_guard<py::gil_scoped_release>());
}
