#pragma once

#include <memory>
#include <string>
#include <string_view>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include <rhbm_gem/core/command/CommandMetadata.hpp>
#include <rhbm_gem/core/command/OptionEnumTraits.hpp>

namespace py = pybind11;

namespace rhbm_gem {
#define RHBM_GEM_COMMAND(COMMAND_ID, COMMAND_TYPE, CLI_NAME, DESCRIPTION, PROFILE, PYTHON_BINDING_NAME) \
    class COMMAND_TYPE;
#include "internal/CommandList.def"
#undef RHBM_GEM_COMMAND
} // namespace rhbm_gem

namespace rhbm_gem::bindings {

template <typename CommandType>
struct CommandBindingName;

#include "internal/CommandBindingNames.generated.inc"

template <typename EnumType>
void BindEnumEntries(py::enum_<EnumType> & py_enum)
{
    for (const auto & entry : EnumOptionTraits<EnumType>::kBindingEntries)
    {
        const std::string name{ entry.token };
        py_enum.value(name.c_str(), entry.value);
    }
}

template <typename CommandType>
py::class_<CommandType, std::shared_ptr<CommandType>> BindCommandClass(py::module_ & module)
{
    const auto binding_name{ std::string(CommandBindingName<CommandType>::value) };
    return py::class_<CommandType, std::shared_ptr<CommandType>>(
        module,
        binding_name.c_str());
}

template <typename CommandType, typename HolderType>
void BindCommonCommandSetters(py::class_<CommandType, HolderType> & py_command)
{
    if constexpr (HasCommonOption(CommandType::kCommonOptions, CommonOption::Threading))
    {
        py_command.def("SetThreadSize", &CommandType::SetThreadSize);
    }
    if constexpr (HasCommonOption(CommandType::kCommonOptions, CommonOption::Verbose))
    {
        py_command.def("SetVerboseLevel", &CommandType::SetVerboseLevel);
    }
    if constexpr (HasCommonOption(CommandType::kCommonOptions, CommonOption::Database))
    {
        py_command.def("SetDatabasePath", &CommandType::SetDatabasePath);
    }
    if constexpr (HasCommonOption(CommandType::kCommonOptions, CommonOption::OutputFolder))
    {
        py_command.def("SetFolderPath", &CommandType::SetFolderPath);
    }
}

template <typename CommandType, typename HolderType>
void BindCommandDiagnostics(py::class_<CommandType, HolderType> & py_command)
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
            [](const CommandType & command)
            {
                return command.GetValidationIssues();
            });
}

template <typename CommandType>
void BindCommand(py::module_ & module);

void BindCommonTypes(py::module_ & module);

} // namespace rhbm_gem::bindings
