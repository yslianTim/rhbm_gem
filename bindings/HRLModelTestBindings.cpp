#include "BindingHelpers.hpp"

#include <rhbm_gem/core/command/HRLModelTestCommand.hpp>

namespace rhbm_gem::bindings {

void BindHRLModelTest(py::module_ & module)
{
    auto model_test{ BindBuiltInCommand<HRLModelTestCommand>(module) };
    model_test
        .def(py::init<const DataIoServices &>())
        .def("Execute", &HRLModelTestCommand::Execute)
        .def("SetThreadSize", &HRLModelTestCommand::SetThreadSize)
        .def("SetTesterChoice", &HRLModelTestCommand::SetTesterChoice)
        .def("SetFitRangeMinimum", &HRLModelTestCommand::SetFitRangeMinimum)
        .def("SetFitRangeMaximum", &HRLModelTestCommand::SetFitRangeMaximum)
        .def("SetAlphaR", &HRLModelTestCommand::SetAlphaR)
        .def("SetAlphaG", &HRLModelTestCommand::SetAlphaG);
    BindCommandDiagnostics(model_test);
}

} // namespace rhbm_gem::bindings
