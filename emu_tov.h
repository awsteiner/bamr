#ifndef EMU_TOV_H
#define EMU_TOV_H

#include <vector>
#include <string>
#include <pybind11/pybind11.h>  // Add this to include pybind11

namespace interpm_py {

class p_wrapper {
public:
    p_wrapper();
    ~p_wrapper();

    void set_data(const std::vector<double>& in_data, const std::vector<double>& out_data);
    std::vector<double> eval(const std::vector<double>& input);
    void set_data_str(const std::vector<double>& in_data, const std::vector<double>& out_data, const std::string& options);

private:
    void initialize_python();
    void finalize_python();
    pybind11::object interpm_tf_dnn;  // Holds the Python class instance
};

}  // namespace interpm_py

class c_wrapper {
public:
    c_wrapper();
    ~c_wrapper();

    void set_data(const std::vector<double>& in_data, const std::vector<double>& out_data);
    std::vector<double> eval(const std::vector<double>& input);
    void set_data_str(const std::vector<double>& in_data, const std::vector<double>& out_data, const std::string& options);

private:
    interpm_py::p_wrapper* pWrapper;
};

#endif  // EMU_TOV_H
