#include <pybind11/pybind11.h>
#include "PotentialAnalysisCommand.hpp"

PYBIND11_MODULE(cpp_module, m)
{
    pybind11::class_<PotentialAnalysisCommand>(m, "PotentialAnalysisCommand")
        .def(pybind11::init<>())
        .def("Execute",                 &PotentialAnalysisCommand::Execute)
        .def("SetFitRangeMinimum",      &PotentialAnalysisCommand::SetFitRangeMinimum)
        .def("SetFitRangeMaximum",      &PotentialAnalysisCommand::SetFitRangeMaximum)
        .def("SetAlphaR",               &PotentialAnalysisCommand::SetAlphaR)
        .def("SetAlphaG",               &PotentialAnalysisCommand::SetAlphaG)
        .def("SetDatabasePath",         &PotentialAnalysisCommand::SetDatabasePath)
        .def("SetModelFilePath",        &PotentialAnalysisCommand::SetModelFilePath)
        .def("SetMapFilePath",          &PotentialAnalysisCommand::SetMapFilePath)
        .def("SetSavedKeyTag",          &PotentialAnalysisCommand::SetSavedKeyTag)
        .def("SetThreadSize",           &PotentialAnalysisCommand::SetThreadSize)
        .def("SetSamplingSize",         &PotentialAnalysisCommand::SetSamplingSize)
        .def("SetSamplingRangeMinimum", &PotentialAnalysisCommand::SetSamplingRangeMinimum)
        .def("SetSamplingRangeMaximum", &PotentialAnalysisCommand::SetSamplingRangeMaximum)
        ;
}