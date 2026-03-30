#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include <rhbm_gem/core/command/CommandApi.hpp>

namespace py = pybind11;

namespace rhbm_gem::bindings {

void BindCommonTypes(py::module_ & module);

namespace {

template <typename Request, typename FieldSpec>
void BindRequestField(py::class_<Request> & py_request, const FieldSpec & field)
{
    py_request.def_readwrite(field.python_name, field.member);
}

template <typename Request>
void BindRequestType(py::module_ & module, const char * type_name)
{
    auto py_request{ py::class_<Request>(module, type_name) };
    py_request.def(py::init<>());
    Request::VisitFields([&](const auto & field)
    {
        BindRequestField(py_request, field);
    });
}

} // namespace

void BindCommandApi(py::module_ & module)
{
    BindRequestType<CommonCommandRequest>(module, "CommonCommandRequest");

#define RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION, PROFILE)                           \
    BindRequestType<COMMAND_ID##Request>(module, #COMMAND_ID "Request");
#include <rhbm_gem/core/command/CommandList.def>
#undef RHBM_GEM_COMMAND

    py::class_<ExecutionReport>(module, "ExecutionReport")
        .def(py::init<>())
        .def_readonly("prepared", &ExecutionReport::prepared)
        .def_readonly("executed", &ExecutionReport::executed)
        .def_readonly("validation_issues", &ExecutionReport::validation_issues);

#define RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION, PROFILE)                           \
    module.def("Run" #COMMAND_ID, &Run##COMMAND_ID);
#include <rhbm_gem/core/command/CommandList.def>
#undef RHBM_GEM_COMMAND
}

} // namespace rhbm_gem::bindings

PYBIND11_MODULE(rhbm_gem_module, module)
{
    rhbm_gem::bindings::BindCommonTypes(module);
    rhbm_gem::bindings::BindCommandApi(module);
}
