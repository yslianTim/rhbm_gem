#pragma once

#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include <rhbm_gem/core/command/OptionEnumTraits.hpp>

#include "internal/BuiltInCommandBindingInternal.hpp"

namespace py = pybind11;

namespace rhbm_gem::bindings {

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
py::class_<CommandType, std::shared_ptr<CommandType>> BindBuiltInCommand(py::module_ & module)
{
    constexpr auto binding_name{ BuiltInCommandBindingName<CommandType>::value };
    return py::class_<CommandType, std::shared_ptr<CommandType>>(
        module,
        std::string(binding_name).c_str());
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

void BindCommonTypes(py::module_ & module);
void BindPotentialAnalysis(py::module_ & module);
void BindPotentialDisplay(py::module_ & module);
void BindResultDump(py::module_ & module);
void BindMapSimulation(py::module_ & module);
void BindMapVisualization(py::module_ & module);
void BindPositionEstimation(py::module_ & module);
void BindHRLModelTest(py::module_ & module);

} // namespace rhbm_gem::bindings
