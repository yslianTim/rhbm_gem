#include "BindingHelpers.hpp"

#include <rhbm_gem/core/command/PositionEstimationCommand.hpp>

namespace rhbm_gem::bindings {

template <>
void BindCommand<PositionEstimationCommand>(py::module_ & module)
{
    auto position_estimation{ BindCommandClass<PositionEstimationCommand>(module) };
    position_estimation
        .def(py::init<>())
        .def("Execute", &PositionEstimationCommand::Execute)
        .def("SetMapFilePath", &PositionEstimationCommand::SetMapFilePath)
        .def("SetIterationCount", &PositionEstimationCommand::SetIterationCount)
        .def("SetKNNSize", &PositionEstimationCommand::SetKNNSize)
        .def("SetAlpha", &PositionEstimationCommand::SetAlpha)
        .def("SetThresholdRatio", &PositionEstimationCommand::SetThresholdRatio)
        .def("SetDedupTolerance", &PositionEstimationCommand::SetDedupTolerance);
    BindCommonCommandSetters(position_estimation);
    BindCommandDiagnostics(position_estimation);
}

} // namespace rhbm_gem::bindings
