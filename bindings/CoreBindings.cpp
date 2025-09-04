#include <pybind11/pybind11.h>
#include <pybind11/stl/filesystem.h>
#include "PotentialAnalysisCommand.hpp"
#include "PotentialDisplayCommand.hpp"
#include "ResultDumpCommand.hpp"
#include "MapSimulationCommand.hpp"
#include "OptionEnumClass.hpp"

PYBIND11_MODULE(cpp_module, m)
{
    pybind11::enum_<PainterType>(m, "PainterType")
        .value("ATOM",       PainterType::ATOM)
        .value("MODEL",      PainterType::MODEL)
        .value("COMPARISON", PainterType::COMPARISON)
        .value("DEMO",       PainterType::DEMO);
    pybind11::implicitly_convertible<int, PainterType>();

    pybind11::enum_<PrinterType>(m, "PrinterType")
        .value("ATOM_POSITION",  PrinterType::ATOM_POSITION)
        .value("MAP_VALUE",      PrinterType::MAP_VALUE)
        .value("GAUS_ESTIMATES", PrinterType::GAUS_ESTIMATES);
    pybind11::implicitly_convertible<int, PrinterType>();

    pybind11::enum_<PotentialModel>(m, "PotentialModel")
        .value("SINGLE_GAUS",      PotentialModel::SINGLE_GAUS)
        .value("FIVE_GAUS_CHARGE", PotentialModel::FIVE_GAUS_CHARGE)
        .value("SINGLE_GAUS_USER", PotentialModel::SINGLE_GAUS_USER);
    pybind11::implicitly_convertible<int, PotentialModel>();

    pybind11::enum_<PartialCharge>(m, "PartialCharge")
        .value("NEUTRAL", PartialCharge::NEUTRAL)
        .value("PARTIAL", PartialCharge::PARTIAL)
        .value("AMBER",   PartialCharge::AMBER);
    pybind11::implicitly_convertible<int, PartialCharge>();

    pybind11::class_<PotentialAnalysisCommand>(m, "PotentialAnalysisCommand")
        .def(pybind11::init<>())
        .def("Execute",                 &PotentialAnalysisCommand::Execute)
        .def("SetAsymmetryFlag",        &PotentialAnalysisCommand::SetAsymmetryFlag)
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
        .def("SetSimulationFlag",       &PotentialAnalysisCommand::SetSimulationFlag)
        .def("SetSimulatedMapResolution",&PotentialAnalysisCommand::SetSimulatedMapResolution)
        ;

    pybind11::class_<PotentialDisplayCommand>(m, "PotentialDisplayCommand")
        .def(pybind11::init<>())
        .def("Execute",                 &PotentialDisplayCommand::Execute)
        .def("SetPainterChoice",        &PotentialDisplayCommand::SetPainterChoice)
        .def("SetModelKeyTagList",      &PotentialDisplayCommand::SetModelKeyTagList)
        .def("SetRefModelKeyTagListMap",&PotentialDisplayCommand::SetRefModelKeyTagListMap)
        .def("SetDatabasePath",         &PotentialDisplayCommand::SetDatabasePath)
        .def("SetFolderPath",           &PotentialDisplayCommand::SetFolderPath)
        .def("SetPickChainID",          &PotentialDisplayCommand::SetPickChainID)
        .def("SetPickResidueType",      &PotentialDisplayCommand::SetPickResidueType)
        .def("SetPickElementType",      &PotentialDisplayCommand::SetPickElementType)
        .def("SetVetoChainID",          &PotentialDisplayCommand::SetVetoChainID)
        .def("SetVetoResidueType",      &PotentialDisplayCommand::SetVetoResidueType)
        .def("SetVetoElementType",      &PotentialDisplayCommand::SetVetoElementType)
        ;

    pybind11::class_<ResultDumpCommand>(m, "ResultDumpCommand")
        .def(pybind11::init<>())
        .def("Execute",                 &ResultDumpCommand::Execute)
        .def("SetPrinterChoice",        &ResultDumpCommand::SetPrinterChoice)
        .def("SetModelKeyTagList",      &ResultDumpCommand::SetModelKeyTagList)
        .def("SetDatabasePath",         &ResultDumpCommand::SetDatabasePath)
        .def("SetMapFilePath",          &ResultDumpCommand::SetMapFilePath)
        .def("SetFolderPath",           &ResultDumpCommand::SetFolderPath)
        ;

    pybind11::class_<MapSimulationCommand>(m, "MapSimulationCommand")
        .def(pybind11::init<>())
        .def("Execute",                 &MapSimulationCommand::Execute)
        .def("SetFolderPath",           &MapSimulationCommand::SetFolderPath)
        .def("SetModelFilePath",        &MapSimulationCommand::SetModelFilePath)
        .def("SetThreadSize",           &MapSimulationCommand::SetThreadSize)
        .def("SetPotentialModelChoice", &MapSimulationCommand::SetPotentialModelChoice)
        .def("SetPartialChargeChoice",  &MapSimulationCommand::SetPartialChargeChoice)
        .def("SetCutoffDistance",       &MapSimulationCommand::SetCutoffDistance)
        .def("SetGridSpacing",          &MapSimulationCommand::SetGridSpacing)
        .def("SetBlurringWidthList",    &MapSimulationCommand::SetBlurringWidthList)
        ;
}
