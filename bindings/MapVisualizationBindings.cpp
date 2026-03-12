#include "BindingHelpers.hpp"

#include <rhbm_gem/core/command/MapVisualizationCommand.hpp>

namespace rhbm_gem::bindings {

template <>
void BindCommand<MapVisualizationCommand>(py::module_ & module)
{
    auto map_visualization{ BindCommandClass<MapVisualizationCommand>(module) };
    map_visualization
        .def(py::init<const DataIoServices &>())
        .def("Execute", &MapVisualizationCommand::Execute)
        .def("SetModelFilePath", &MapVisualizationCommand::SetModelFilePath)
        .def("SetMapFilePath", &MapVisualizationCommand::SetMapFilePath)
        .def("SetAtomSerialID", &MapVisualizationCommand::SetAtomSerialID)
        .def("SetSamplingSize", &MapVisualizationCommand::SetSamplingSize)
        .def("SetWindowSize", &MapVisualizationCommand::SetWindowSize);
    BindCommonCommandSetters(map_visualization);
    BindCommandDiagnostics(map_visualization);
}

} // namespace rhbm_gem::bindings
