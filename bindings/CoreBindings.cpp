#include <pybind11/pybind11.h>
#include <pybind11/stl/filesystem.h>
#include <string>

#include "PotentialAnalysisCommand.hpp"
#include "PotentialDisplayCommand.hpp"
#include "ResultDumpCommand.hpp"
#include "MapSimulationCommand.hpp"
#include "OptionEnumTraits.hpp"
#include "OptionEnumClass.hpp"

namespace {

template <typename EnumType>
void BindEnumEntries(pybind11::enum_<EnumType> & py_enum)
{
    for (const auto & entry : rhbm_gem::EnumOptionTraits<EnumType>::kBindingEntries)
    {
        const std::string name{ entry.token };
        py_enum.value(name.c_str(), entry.value);
    }
}

} // namespace

PYBIND11_MODULE(rhbm_gem_module, m)
{
    auto painter_type{ pybind11::enum_<rhbm_gem::PainterType>(m, "PainterType") };
    BindEnumEntries(painter_type);
    pybind11::implicitly_convertible<int, rhbm_gem::PainterType>();

    auto printer_type{ pybind11::enum_<rhbm_gem::PrinterType>(m, "PrinterType") };
    BindEnumEntries(printer_type);
    pybind11::implicitly_convertible<int, rhbm_gem::PrinterType>();

    auto potential_model{ pybind11::enum_<rhbm_gem::PotentialModel>(m, "PotentialModel") };
    BindEnumEntries(potential_model);
    pybind11::implicitly_convertible<int, rhbm_gem::PotentialModel>();

    auto partial_charge{ pybind11::enum_<rhbm_gem::PartialCharge>(m, "PartialCharge") };
    BindEnumEntries(partial_charge);
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
