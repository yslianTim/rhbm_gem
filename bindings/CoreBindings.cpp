#include <pybind11/pybind11.h>
#include <pybind11/stl/filesystem.h>
#include "PotentialAnalysisCommand.hpp"
#include "PotentialDisplayCommand.hpp"
#include "ResultDumpCommand.hpp"
#include "MapSimulationCommand.hpp"
#include "OptionEnumClass.hpp"

PYBIND11_MODULE(rhbm_gem_module, m)
{
    pybind11::enum_<rhbm_gem::PainterType>(m, "PainterType")
        .value("GAUS",       rhbm_gem::PainterType::GAUS)
        .value("ATOM",       rhbm_gem::PainterType::ATOM)
        .value("MODEL",      rhbm_gem::PainterType::MODEL)
        .value("COMPARISON", rhbm_gem::PainterType::COMPARISON)
        .value("DEMO",       rhbm_gem::PainterType::DEMO);
    pybind11::implicitly_convertible<int, rhbm_gem::PainterType>();

    pybind11::enum_<rhbm_gem::PrinterType>(m, "PrinterType")
        .value("ATOM_POSITION",  rhbm_gem::PrinterType::ATOM_POSITION)
        .value("MAP_VALUE",      rhbm_gem::PrinterType::MAP_VALUE)
        .value("GAUS_ESTIMATES", rhbm_gem::PrinterType::GAUS_ESTIMATES);
    pybind11::implicitly_convertible<int, rhbm_gem::PrinterType>();

    pybind11::enum_<rhbm_gem::PotentialModel>(m, "PotentialModel")
        .value("SINGLE_GAUS",      rhbm_gem::PotentialModel::SINGLE_GAUS)
        .value("FIVE_GAUS_CHARGE", rhbm_gem::PotentialModel::FIVE_GAUS_CHARGE)
        .value("SINGLE_GAUS_USER", rhbm_gem::PotentialModel::SINGLE_GAUS_USER);
    pybind11::implicitly_convertible<int, rhbm_gem::PotentialModel>();

    pybind11::enum_<rhbm_gem::PartialCharge>(m, "PartialCharge")
        .value("NEUTRAL", rhbm_gem::PartialCharge::NEUTRAL)
        .value("PARTIAL", rhbm_gem::PartialCharge::PARTIAL)
        .value("AMBER",   rhbm_gem::PartialCharge::AMBER);
    pybind11::implicitly_convertible<int, rhbm_gem::PartialCharge>();

    pybind11::class_<rhbm_gem::PotentialAnalysisCommand>(m, "PotentialAnalysisCommand")
        .def(pybind11::init<>())
        .def("Execute",                  &rhbm_gem::PotentialAnalysisCommand::Execute)
        .def("SetAsymmetryFlag",         &rhbm_gem::PotentialAnalysisCommand::SetAsymmetryFlag)
        .def("SetFitRangeMinimum",       &rhbm_gem::PotentialAnalysisCommand::SetFitRangeMinimum)
        .def("SetFitRangeMaximum",       &rhbm_gem::PotentialAnalysisCommand::SetFitRangeMaximum)
        .def("SetAlphaR",                &rhbm_gem::PotentialAnalysisCommand::SetAlphaR)
        .def("SetAlphaG",                &rhbm_gem::PotentialAnalysisCommand::SetAlphaG)
        .def("SetDatabasePath",          &rhbm_gem::PotentialAnalysisCommand::SetDatabasePath)
        .def("SetModelFilePath",         &rhbm_gem::PotentialAnalysisCommand::SetModelFilePath)
        .def("SetMapFilePath",           &rhbm_gem::PotentialAnalysisCommand::SetMapFilePath)
        .def("SetSavedKeyTag",           &rhbm_gem::PotentialAnalysisCommand::SetSavedKeyTag)
        .def("SetThreadSize",            &rhbm_gem::PotentialAnalysisCommand::SetThreadSize)
        .def("SetSamplingSize",          &rhbm_gem::PotentialAnalysisCommand::SetSamplingSize)
        .def("SetSamplingRangeMinimum",  &rhbm_gem::PotentialAnalysisCommand::SetSamplingRangeMinimum)
        .def("SetSamplingRangeMaximum",  &rhbm_gem::PotentialAnalysisCommand::SetSamplingRangeMaximum)
        .def("SetSimulationFlag",        &rhbm_gem::PotentialAnalysisCommand::SetSimulationFlag)
        .def("SetSimulatedMapResolution",&rhbm_gem::PotentialAnalysisCommand::SetSimulatedMapResolution)
        ;

    pybind11::class_<rhbm_gem::PotentialDisplayCommand>(m, "PotentialDisplayCommand")
        .def(pybind11::init<>())
        .def("Execute",                 &rhbm_gem::PotentialDisplayCommand::Execute)
        .def("SetPainterChoice",        &rhbm_gem::PotentialDisplayCommand::SetPainterChoice)
        .def("SetModelKeyTagList",      &rhbm_gem::PotentialDisplayCommand::SetModelKeyTagList)
        .def("SetRefModelKeyTagListMap",&rhbm_gem::PotentialDisplayCommand::SetRefModelKeyTagListMap)
        .def("SetDatabasePath",         &rhbm_gem::PotentialDisplayCommand::SetDatabasePath)
        .def("SetFolderPath",           &rhbm_gem::PotentialDisplayCommand::SetFolderPath)
        .def("SetPickChainID",          &rhbm_gem::PotentialDisplayCommand::SetPickChainID)
        .def("SetPickResidueType",      &rhbm_gem::PotentialDisplayCommand::SetPickResidueType)
        .def("SetPickElementType",      &rhbm_gem::PotentialDisplayCommand::SetPickElementType)
        .def("SetVetoChainID",          &rhbm_gem::PotentialDisplayCommand::SetVetoChainID)
        .def("SetVetoResidueType",      &rhbm_gem::PotentialDisplayCommand::SetVetoResidueType)
        .def("SetVetoElementType",      &rhbm_gem::PotentialDisplayCommand::SetVetoElementType)
        ;

    pybind11::class_<rhbm_gem::ResultDumpCommand>(m, "ResultDumpCommand")
        .def(pybind11::init<>())
        .def("Execute",                 &rhbm_gem::ResultDumpCommand::Execute)
        .def("SetPrinterChoice",        &rhbm_gem::ResultDumpCommand::SetPrinterChoice)
        .def("SetModelKeyTagList",      &rhbm_gem::ResultDumpCommand::SetModelKeyTagList)
        .def("SetDatabasePath",         &rhbm_gem::ResultDumpCommand::SetDatabasePath)
        .def("SetMapFilePath",          &rhbm_gem::ResultDumpCommand::SetMapFilePath)
        .def("SetFolderPath",           &rhbm_gem::ResultDumpCommand::SetFolderPath)
        ;

    pybind11::class_<rhbm_gem::MapSimulationCommand>(m, "MapSimulationCommand")
        .def(pybind11::init<>())
        .def("Execute",                 &rhbm_gem::MapSimulationCommand::Execute)
        .def("SetFolderPath",           &rhbm_gem::MapSimulationCommand::SetFolderPath)
        .def("SetModelFilePath",        &rhbm_gem::MapSimulationCommand::SetModelFilePath)
        .def("SetThreadSize",           &rhbm_gem::MapSimulationCommand::SetThreadSize)
        .def("SetPotentialModelChoice", &rhbm_gem::MapSimulationCommand::SetPotentialModelChoice)
        .def("SetPartialChargeChoice",  &rhbm_gem::MapSimulationCommand::SetPartialChargeChoice)
        .def("SetCutoffDistance",       &rhbm_gem::MapSimulationCommand::SetCutoffDistance)
        .def("SetGridSpacing",          &rhbm_gem::MapSimulationCommand::SetGridSpacing)
        .def("SetBlurringWidthList",    &rhbm_gem::MapSimulationCommand::SetBlurringWidthList)
        ;
}
