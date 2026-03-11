#include "BindingHelpers.hpp"

#include <rhbm_gem/core/command/ResultDumpCommand.hpp>

namespace rhbm_gem::bindings {

template <>
void BindCommand<ResultDumpCommand>(py::module_ & module)
{
    auto result_dump{ BindBuiltInCommand<ResultDumpCommand>(module) };
    result_dump
        .def(py::init<const DataIoServices &>())
        .def("Execute", &ResultDumpCommand::Execute)
        .def("SetPrinterChoice", &ResultDumpCommand::SetPrinterChoice)
        .def("SetModelKeyTagList", &ResultDumpCommand::SetModelKeyTagList)
        .def("SetMapFilePath", &ResultDumpCommand::SetMapFilePath);
    BindCommonCommandSetters(result_dump);
    BindCommandDiagnostics(result_dump);
}

} // namespace rhbm_gem::bindings
