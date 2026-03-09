#include "BindingHelpers.hpp"

#include <rhbm_gem/core/command/PositionEstimationCommand.hpp>

namespace rhbm_gem::bindings {

void BindPositionEstimation(py::module_ & module)
{
    auto position_estimation{ BindBuiltInCommand<PositionEstimationCommand>(module) };
    position_estimation
        .def(py::init<>())
        .def("Execute", &PositionEstimationCommand::Execute)
        .def("SetFolderPath", &PositionEstimationCommand::SetFolderPath)
        .def("SetThreadSize", &PositionEstimationCommand::SetThreadSize)
        .def("SetMapFilePath", &PositionEstimationCommand::SetMapFilePath)
        .def("SetIterationCount", &PositionEstimationCommand::SetIterationCount)
        .def("SetKNNSize", &PositionEstimationCommand::SetKNNSize)
        .def("SetAlpha", &PositionEstimationCommand::SetAlpha)
        .def("SetThresholdRatio", &PositionEstimationCommand::SetThresholdRatio)
        .def("SetDedupTolerance", &PositionEstimationCommand::SetDedupTolerance);
    BindCommandDiagnostics(position_estimation);
}

} // namespace rhbm_gem::bindings
