#include "BindingHelpers.hpp"

#include <rhbm_gem/core/command/PotentialAnalysisCommand.hpp>

namespace rhbm_gem::bindings {

void BindPotentialAnalysis(py::module_ & module)
{
    auto potential_analysis{ BindBuiltInCommand<PotentialAnalysisCommand>(module) };
    potential_analysis
        .def(py::init<const DataIoServices &>())
        .def("Execute", &PotentialAnalysisCommand::Execute)
        .def("SetAsymmetryFlag", &PotentialAnalysisCommand::SetAsymmetryFlag)
        .def("SetFitRangeMinimum", &PotentialAnalysisCommand::SetFitRangeMinimum)
        .def("SetFitRangeMaximum", &PotentialAnalysisCommand::SetFitRangeMaximum)
        .def("SetAlphaR", &PotentialAnalysisCommand::SetAlphaR)
        .def("SetAlphaG", &PotentialAnalysisCommand::SetAlphaG)
        .def("SetDatabasePath", &PotentialAnalysisCommand::SetDatabasePath)
        .def("SetModelFilePath", &PotentialAnalysisCommand::SetModelFilePath)
        .def("SetMapFilePath", &PotentialAnalysisCommand::SetMapFilePath)
        .def("SetSavedKeyTag", &PotentialAnalysisCommand::SetSavedKeyTag)
        .def("SetTrainingReportDir", &PotentialAnalysisCommand::SetTrainingReportDir)
        .def("SetThreadSize", &PotentialAnalysisCommand::SetThreadSize)
        .def("SetSamplingSize", &PotentialAnalysisCommand::SetSamplingSize)
        .def("SetSamplingRangeMinimum", &PotentialAnalysisCommand::SetSamplingRangeMinimum)
        .def("SetSamplingRangeMaximum", &PotentialAnalysisCommand::SetSamplingRangeMaximum)
        .def("SetSimulationFlag", &PotentialAnalysisCommand::SetSimulationFlag)
        .def("SetSimulatedMapResolution", &PotentialAnalysisCommand::SetSimulatedMapResolution);
    BindCommandDiagnostics(potential_analysis);
}

} // namespace rhbm_gem::bindings
