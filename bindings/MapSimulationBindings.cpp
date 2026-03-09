#include "BindingHelpers.hpp"

#include <rhbm_gem/core/command/MapSimulationCommand.hpp>

namespace rhbm_gem::bindings {

void BindMapSimulation(py::module_ & module)
{
    auto map_simulation{ BindBuiltInCommand<MapSimulationCommand>(module) };
    map_simulation
        .def(py::init<>())
        .def("Execute", &MapSimulationCommand::Execute)
        .def("SetFolderPath", &MapSimulationCommand::SetFolderPath)
        .def("SetModelFilePath", &MapSimulationCommand::SetModelFilePath)
        .def("SetThreadSize", &MapSimulationCommand::SetThreadSize)
        .def("SetPotentialModelChoice", &MapSimulationCommand::SetPotentialModelChoice)
        .def("SetPartialChargeChoice", &MapSimulationCommand::SetPartialChargeChoice)
        .def("SetCutoffDistance", &MapSimulationCommand::SetCutoffDistance)
        .def("SetGridSpacing", &MapSimulationCommand::SetGridSpacing)
        .def("SetBlurringWidthList", &MapSimulationCommand::SetBlurringWidthList);
    BindCommandDiagnostics(map_simulation);
}

} // namespace rhbm_gem::bindings
