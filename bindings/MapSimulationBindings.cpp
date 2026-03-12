#include "BindingHelpers.hpp"

#include <rhbm_gem/core/command/MapSimulationCommand.hpp>

namespace rhbm_gem::bindings {

template <>
void BindCommand<MapSimulationCommand>(py::module_ & module)
{
    auto map_simulation{ BindCommandClass<MapSimulationCommand>(module) };
    map_simulation
        .def(py::init<>())
        .def("Execute", &MapSimulationCommand::Execute)
        .def("SetModelFilePath", &MapSimulationCommand::SetModelFilePath)
        .def("SetPotentialModelChoice", &MapSimulationCommand::SetPotentialModelChoice)
        .def("SetPartialChargeChoice", &MapSimulationCommand::SetPartialChargeChoice)
        .def("SetCutoffDistance", &MapSimulationCommand::SetCutoffDistance)
        .def("SetGridSpacing", &MapSimulationCommand::SetGridSpacing)
        .def("SetBlurringWidthList", &MapSimulationCommand::SetBlurringWidthList);
    BindCommonCommandSetters(map_simulation);
    BindCommandDiagnostics(map_simulation);
}

} // namespace rhbm_gem::bindings
