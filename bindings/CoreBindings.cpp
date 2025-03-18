#include <pybind11/pybind11.h>
#include "PotentialAnalysisCommand.hpp"

namespace py = pybind11;

PYBIND11_MODULE(cpp_module, m)
{
    py::class_<PotentialAnalysisCommand>(m, "PotentialAnalysisCommand")
        .def(py::init<>())
        .def("Execute", &PotentialAnalysisCommand::Execute)
        .def("SetFitRangeMinimum", &PotentialAnalysisCommand::SetFitRangeMinimum)
        .def("SetFitRangeMaximum", &PotentialAnalysisCommand::SetFitRangeMaximum)
        .def("SetAlphaR", &PotentialAnalysisCommand::SetAlphaR)
        .def("SetAlphaG", &PotentialAnalysisCommand::SetAlphaG)
        .def("SetDatabasePath", &PotentialAnalysisCommand::SetDatabasePath)
        .def("SetModelFilePath", &PotentialAnalysisCommand::SetModelFilePath)
        .def("SetMapFilePath", &PotentialAnalysisCommand::SetMapFilePath)
        .def("SetSavedKeyTag", &PotentialAnalysisCommand::SetSavedKeyTag)
        .def("SetThreadSize", &PotentialAnalysisCommand::SetThreadSize)
        .def("SetSamplingSize", &PotentialAnalysisCommand::SetSamplingSize)
        .def("SetSamplingRangeMinimum", &PotentialAnalysisCommand::SetSamplingRangeMinimum)
        .def("SetSamplingRangeMaximum", &PotentialAnalysisCommand::SetSamplingRangeMaximum)
        .def("SetPickChainID", &PotentialAnalysisCommand::SetPickChainID)
        .def("SetPickIndicator", &PotentialAnalysisCommand::SetPickIndicator)
        .def("SetPickResidueType", &PotentialAnalysisCommand::SetPickResidueType)
        .def("SetPickElementType", &PotentialAnalysisCommand::SetPickElementType)
        .def("SetPickRemotenessType", &PotentialAnalysisCommand::SetPickRemotenessType)
        .def("SetPickBranchType", &PotentialAnalysisCommand::SetPickBranchType)
        .def("SetVetoChainID", &PotentialAnalysisCommand::SetVetoChainID)
        .def("SetVetoIndicator", &PotentialAnalysisCommand::SetVetoIndicator)
        .def("SetVetoResidueType", &PotentialAnalysisCommand::SetVetoResidueType)
        .def("SetVetoElementType", &PotentialAnalysisCommand::SetVetoElementType)
        .def("SetVetoRemotenessType", &PotentialAnalysisCommand::SetVetoRemotenessType)
        .def("SetVetoBranchType", &PotentialAnalysisCommand::SetVetoBranchType)
        ;
}