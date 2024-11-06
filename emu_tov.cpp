#include "emu_tov.h"
#include <pybind11/embed.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace interpm_py {

p_wrapper::p_wrapper() {
    initialize_python();
    py::module_ emu_tov = py::module_::import("emu_tov");  // Load Python module
    interpm_tf_dnn = emu_tov.attr("interpm_tf_dnn")();     // Instantiate the Python class
}

p_wrapper::~p_wrapper() {
    finalize_python();
}

void p_wrapper::initialize_python() {
    py::initialize_interpreter();  // Start the Python interpreter
}

void p_wrapper::finalize_python() {
    py::finalize_interpreter();    // Stop the Python interpreter
}

void p_wrapper::set_data(const std::vector<double>& in_data, const std::vector<double>& out_data) {
    interpm_tf_dnn.attr("set_data")(in_data, out_data);  // Call Python set_data method
}

std::vector<double> p_wrapper::eval(const std::vector<double>& input) {
    return interpm_tf_dnn.attr("eval")(input).cast<std::vector<double>>();  // Call Python eval and return results
}

void p_wrapper::set_data_str(const std::vector<double>& in_data, const std::vector<double>& out_data, const std::string& options) {
    interpm_tf_dnn.attr("set_data_str")(in_data, out_data, options);  // Call Python set_data_str method
}

}  // namespace interpm_py

// c_wrapper Implementation

c_wrapper::c_wrapper() : pWrapper(new interpm_py::p_wrapper()) {}

c_wrapper::~c_wrapper() {
    delete pWrapper;
}

void c_wrapper::set_data(const std::vector<double>& in_data, const std::vector<double>& out_data) {
    pWrapper->set_data(in_data, out_data);
}

std::vector<double> c_wrapper::eval(const std::vector<double>& input) {
    return pWrapper->eval(input);
}

void c_wrapper::set_data_str(const std::vector<double>& in_data, const std::vector<double>& out_data, const std::string& options) {
    pWrapper->set_data_str(in_data, out_data, options);
}
