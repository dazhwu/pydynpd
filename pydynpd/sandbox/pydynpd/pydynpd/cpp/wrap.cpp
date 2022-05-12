#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "dynpd.h"
#include "info.h"
#include "GMM.h"
#include "List_Variables.h"
#include "variable.h"
#include "Command.h"
#include "instruments.h"
#include "Common_Functions.h"
#include "dynamic_model.h"

namespace py = pybind11;

/*

struct _options {
    int steps = 2;
    bool level = true;
    bool beginner = false;
    bool timedumm = false;
    bool collapse = false;
    string mmsc = "bic";
    string transformation = "fd";

*/

PYBIND11_MODULE(gmm_module, m)
{

    using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    py::class_<model_options>(m, "model_options")
        .def(py::init<>())
        .def_readwrite("steps", &model_options::steps)
        .def_readwrite("level", &model_options::level)
        .def_readwrite("beginner", &model_options::beginner)
        .def_readwrite("timedumm", &model_options::timedumm)
        .def_readwrite("collapse", &model_options::collapse)
        .def_readwrite("mmsc", &model_options::mmsc)
        .def_readwrite("transformation", &model_options::transformation);

    py::class_<Regression>(m, "Regression")
        .def(py::init<RowMatrixXd &, RowMatrixXd &, RowMatrixXd &, Ref<RowMatrixXd>, Ref<RowMatrixXd>, Ref<RowMatrixXd>>())
        .def_readwrite("regression_table", &Regression::regression_table)
        .def_readwrite("vcov", &Regression::vcov)
        .def_readwrite("weighting", &Regression::weighting)
        .def_readwrite("Z", &Regression::Z)
        .def_readwrite("Y", &Regression::Cy)
        .def_readwrite("X", &Regression::Cx);

    py::class_<basic_info>(m, "basic_info")
        .def(py::init<int, int, int, int, int, int, int, int, double, vector<regular_variable>, bool>())
        .def_readwrite("dep", &basic_info::dep)
        .def_readwrite("indep", &basic_info::indep)
        .def_readwrite("N", &basic_info::N)
        .def_readwrite("T", &basic_info::T)
        .def_readwrite("num_obs", &basic_info::num_obs)
        .def_readwrite("num_instr", &basic_info::num_instr)
        .def_readwrite("num_indep", &basic_info::num_indep)
        .def_readwrite("diff_width", &basic_info::diff_width)
        .def_readwrite("max_obs", &basic_info::max_obs)
        .def_readwrite("min_obs", &basic_info::min_obs)
        .def_readwrite("avg_obs", &basic_info::avg_obs)
        .def_readwrite("SS", &basic_info::SS);

    py::class_<z_info>(m, "z_info")
        .def(py::init<>())
        .def(py::init<int, int, int, int, int, int, int, int>())
        .def_readwrite("diff_width", &z_info::diff_width)
        .def_readwrite("diff_height", &z_info::diff_height)
        .def_readwrite("level_width", &z_info::level_width)
        .def_readwrite("level_height", &z_info::level_height)
        .def_readwrite("z_width", &z_info::z_width)
        .def_readwrite("z_height", &z_info::z_height)

        .def_readwrite("num_Dgmm_instr", &z_info::num_Dgmm_instr)
        .def_readwrite("num_Lgmm_instr", &z_info::num_Lgmm_instr)
        .def_readwrite("num_instr", &z_info::num_instr);

    py::class_<List_Variables>(m, "List_Variables")
        .def(py::init<>())
        .def_readwrite("names", &List_Variables::names)
        .def_readwrite("lags", &List_Variables::lags)
        .def_readwrite("min_lags", &List_Variables::min_lags)
        .def_readwrite("max_lags", &List_Variables::max_lags)
        .def_readwrite("adjustable_min_lags", &List_Variables::adjustable_min_lags)
        .def_readwrite("adjustable_max_lags", &List_Variables::adjustable_max_lags);

    py::class_<Step_Result>(m, "Step_Result")
        .def(py::init<RowMatrixXd, RowMatrixXd, RowMatrixXd, RowMatrixXd,
                      RowMatrixXd, RowMatrixXd, RowMatrixXd, RowMatrixXd, RowMatrixXd, RowMatrixXd,
                      RowMatrixXd, RowMatrixXd, RowMatrixXd, RowMatrixXd, double>())
        .def_readwrite("residual", &Step_Result::residual)
        .def_readwrite("_residual_t", &Step_Result::_residual_t)
        .def_readwrite("XZ_W", &Step_Result::XZ_W)
        .def_readwrite("W", &Step_Result::W)
        .def_readwrite("W_inv", &Step_Result::W_inv)
        .def_readwrite("W_next", &Step_Result::W_next)
        .def_readwrite("_M_XZ_W", &Step_Result::_M_XZ_W)
        .def_readwrite("zs", &Step_Result::zs)
        .def_readwrite("ZuuZ", &Step_Result::ZuuZ)
        .def_readwrite("vcov", &Step_Result::vcov)
        .def_readwrite("M", &Step_Result::M)
        .def_readwrite("_zs_list", &Step_Result::_zs_list)
        .def_readwrite("beta", &Step_Result::beta)
        .def_readwrite("std_err", &Step_Result::std_err)
        .def_readwrite("SS", &Step_Result::SS);

    py::class_<AR_test_info>(m, "AR_test_info")
        .def(py::init<int, double, double>())
        .def_readwrite("lag", &AR_test_info::lag)
        .def_readwrite("AR", &AR_test_info::AR)
        .def_readwrite("P_value", &AR_test_info::P_value);

    py::class_<Hansen_test_info>(m, "Hansen_test_info")
        .def(py::init<double, int, double, double>())
        .def_readwrite("test_value", &Hansen_test_info::test_value)
        .def_readwrite("df", &Hansen_test_info::df)
        .def_readwrite("critical_value", &Hansen_test_info::critical_value)
        .def_readwrite("P_value", &Hansen_test_info::P_value);

    py::class_<df_info>(m, "df_info")
        .def_readwrite("N", &df_info::N)
        .def_readwrite("T", &df_info::T)
        //.def_readwrite("ids", &df_info::ids)
        .def_readwrite("first_diff_index", &df_info::first_diff_index)
        .def_readwrite("last_diff_index", &df_info::last_diff_index)
        .def_readwrite("first_level_index", &df_info::first_level_index)
        .def_readwrite("last_level_index", &df_info::last_level_index)
        .def_readwrite("max_lag", &df_info::max_lag);

    py::class_<regular_variable>(m, "regular_variable")
        //.def(py::init<string, int>())
        .def_readwrite("name", &regular_variable::name)
        .def_readwrite("lag", &regular_variable::lag);

    py::class_<gmm_var>(m, "gmm_var")
        //.def(py::init<string, int, int, int>())
        .def_readwrite("name", &gmm_var::name)
        .def_readwrite("min_lag", &gmm_var::min_lag)
        .def_readwrite("max_lag", &gmm_var::max_lag)
        .def_readwrite("lag", &gmm_var::lag);

    m.def("AR_test", &AR_test);
    m.def("process_command", &process_command);
    m.def("regular_process", &regular_process);
    m.def("get_z_table", &get_z_table);
    m.def("prepare_data", &prepare_data);
}

/*
g++  -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` -I/usr/local/include -I/usr/include/eigen3 `python3-config --ldflags` -o gmm_module.so GMM.cpp Specification_Test.cpp Step_Result.cpp Common_Functions.cpp wrap.cpp Command.cpp List_Variables.cpp instruments.cpp dynamic_model.cpp variable.cpp



*/