#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>
#include <stdexcept>
#include <string>

#include "BuiltInCommandCatalog.hpp"
#include "MapSimulationCommand.hpp"
#include "OptionEnumClass.hpp"
#include "OptionEnumTraits.hpp"
#include "PotentialAnalysisCommand.hpp"
#include "PotentialDisplayCommand.hpp"
#include "ResultDumpCommand.hpp"

namespace py = pybind11;

namespace {

template <typename EnumType>
void BindEnumEntries(py::enum_<EnumType> & py_enum)
{
    for (const auto & entry : rhbm_gem::EnumOptionTraits<EnumType>::kBindingEntries)
    {
        const std::string name{ entry.token };
        py_enum.value(name.c_str(), entry.value);
    }
}

template <typename CommandType>
py::class_<CommandType> BindPythonPublicCommand(
    py::module_ & module,
    rhbm_gem::CommandId command_id)
{
    const auto & descriptor{ rhbm_gem::FindCommandDescriptor(command_id) };
    if (!rhbm_gem::IsPythonPublic(descriptor.binding_exposure)
        || descriptor.python_binding_name.empty())
    {
        throw std::runtime_error(
            "Command descriptor is not configured for Python public binding: "
            + std::string(descriptor.name));
    }

    return py::class_<CommandType>(
        module,
        std::string(descriptor.python_binding_name).c_str());
}

template <typename CommandType>
void BindCommandDiagnostics(py::class_<CommandType> & py_command)
{
    py_command
        .def("PrepareForExecution",
            [](CommandType & command)
            {
                return command.PrepareForExecution();
            })
        .def("HasValidationErrors",
            [](const CommandType & command)
            {
                return command.HasValidationErrors();
            })
        .def(
            "GetValidationIssues",
            [](const CommandType & command) -> const std::vector<rhbm_gem::ValidationIssue> &
            {
                return command.GetValidationIssues();
            },
            py::return_value_policy::reference_internal);
}

} // namespace

PYBIND11_MODULE(rhbm_gem_module, m)
{
    py::enum_<LogLevel>(m, "LogLevel")
        .value("Error", LogLevel::Error)
        .value("Warning", LogLevel::Warning)
        .value("Notice", LogLevel::Notice)
        .value("Info", LogLevel::Info)
        .value("Debug", LogLevel::Debug);

    py::enum_<rhbm_gem::ValidationPhase>(m, "ValidationPhase")
        .value("Parse", rhbm_gem::ValidationPhase::Parse)
        .value("Prepare", rhbm_gem::ValidationPhase::Prepare)
        .value("Runtime", rhbm_gem::ValidationPhase::Runtime);

    py::class_<rhbm_gem::ValidationIssue>(m, "ValidationIssue")
        .def_readonly("option_name", &rhbm_gem::ValidationIssue::option_name)
        .def_readonly("phase", &rhbm_gem::ValidationIssue::phase)
        .def_readonly("level", &rhbm_gem::ValidationIssue::level)
        .def_readonly("message", &rhbm_gem::ValidationIssue::message)
        .def_readonly("auto_corrected", &rhbm_gem::ValidationIssue::auto_corrected);

    auto painter_type{ py::enum_<rhbm_gem::PainterType>(m, "PainterType") };
    BindEnumEntries(painter_type);
    py::implicitly_convertible<int, rhbm_gem::PainterType>();

    auto printer_type{ py::enum_<rhbm_gem::PrinterType>(m, "PrinterType") };
    BindEnumEntries(printer_type);
    py::implicitly_convertible<int, rhbm_gem::PrinterType>();

    auto potential_model{ py::enum_<rhbm_gem::PotentialModel>(m, "PotentialModel") };
    BindEnumEntries(potential_model);
    py::implicitly_convertible<int, rhbm_gem::PotentialModel>();

    auto partial_charge{ py::enum_<rhbm_gem::PartialCharge>(m, "PartialCharge") };
    BindEnumEntries(partial_charge);
    py::implicitly_convertible<int, rhbm_gem::PartialCharge>();

    auto potential_analysis{
        BindPythonPublicCommand<rhbm_gem::PotentialAnalysisCommand>(
            m, rhbm_gem::CommandId::PotentialAnalysis)
    };
    potential_analysis
        .def(py::init<>())
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
        .def("SetSimulatedMapResolution",&rhbm_gem::PotentialAnalysisCommand::SetSimulatedMapResolution);
    BindCommandDiagnostics(potential_analysis);

    auto potential_display{
        BindPythonPublicCommand<rhbm_gem::PotentialDisplayCommand>(
            m, rhbm_gem::CommandId::PotentialDisplay)
    };
    potential_display
        .def(py::init<>())
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
        .def("SetVetoElementType",      &rhbm_gem::PotentialDisplayCommand::SetVetoElementType);
    BindCommandDiagnostics(potential_display);

    auto result_dump{
        BindPythonPublicCommand<rhbm_gem::ResultDumpCommand>(
            m, rhbm_gem::CommandId::ResultDump)
    };
    result_dump
        .def(py::init<>())
        .def("Execute",                 &rhbm_gem::ResultDumpCommand::Execute)
        .def("SetPrinterChoice",        &rhbm_gem::ResultDumpCommand::SetPrinterChoice)
        .def("SetModelKeyTagList",      &rhbm_gem::ResultDumpCommand::SetModelKeyTagList)
        .def("SetDatabasePath",         &rhbm_gem::ResultDumpCommand::SetDatabasePath)
        .def("SetMapFilePath",          &rhbm_gem::ResultDumpCommand::SetMapFilePath)
        .def("SetFolderPath",           &rhbm_gem::ResultDumpCommand::SetFolderPath);
    BindCommandDiagnostics(result_dump);

    auto map_simulation{
        BindPythonPublicCommand<rhbm_gem::MapSimulationCommand>(
            m, rhbm_gem::CommandId::MapSimulation)
    };
    map_simulation
        .def(py::init<>())
        .def("Execute",                 &rhbm_gem::MapSimulationCommand::Execute)
        .def("SetFolderPath",           &rhbm_gem::MapSimulationCommand::SetFolderPath)
        .def("SetModelFilePath",        &rhbm_gem::MapSimulationCommand::SetModelFilePath)
        .def("SetThreadSize",           &rhbm_gem::MapSimulationCommand::SetThreadSize)
        .def("SetPotentialModelChoice", &rhbm_gem::MapSimulationCommand::SetPotentialModelChoice)
        .def("SetPartialChargeChoice",  &rhbm_gem::MapSimulationCommand::SetPartialChargeChoice)
        .def("SetCutoffDistance",       &rhbm_gem::MapSimulationCommand::SetCutoffDistance)
        .def("SetGridSpacing",          &rhbm_gem::MapSimulationCommand::SetGridSpacing)
        .def("SetBlurringWidthList",    &rhbm_gem::MapSimulationCommand::SetBlurringWidthList);
    BindCommandDiagnostics(map_simulation);
}
