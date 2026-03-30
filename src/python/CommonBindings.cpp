#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include <rhbm_gem/core/command/CommandContract.hpp>
#include <rhbm_gem/core/command/CommandEnumClass.hpp>

namespace py = pybind11;

namespace rhbm_gem::bindings {
namespace {

template <typename EnumType>
void BindEnumEntries(py::enum_<EnumType> & py_enum)
{
    const auto binding_entries{ GetCommandEnumBindingEntries<EnumType>() };
    for (const auto & entry : binding_entries)
    {
        const std::string name{ entry.token };
        py_enum.value(name.c_str(), entry.value);
    }
}

} // namespace

void BindCommonTypes(py::module_ & module)
{
    py::enum_<LogLevel>(module, "LogLevel")
        .value("Error", LogLevel::Error)
        .value("Warning", LogLevel::Warning)
        .value("Notice", LogLevel::Notice)
        .value("Info", LogLevel::Info)
        .value("Debug", LogLevel::Debug);

    py::enum_<ValidationPhase>(module, "ValidationPhase")
        .value("Parse", ValidationPhase::Parse)
        .value("Prepare", ValidationPhase::Prepare);

    py::class_<ValidationIssue>(module, "ValidationIssue")
        .def_readonly("option_name", &ValidationIssue::option_name)
        .def_readonly("phase", &ValidationIssue::phase)
        .def_readonly("level", &ValidationIssue::level)
        .def_readonly("message", &ValidationIssue::message)
        .def_readonly("auto_corrected", &ValidationIssue::auto_corrected);

    auto painter_type{ py::enum_<PainterType>(module, "PainterType") };
    BindEnumEntries(painter_type);
    py::implicitly_convertible<int, PainterType>();

    auto printer_type{ py::enum_<PrinterType>(module, "PrinterType") };
    BindEnumEntries(printer_type);
    py::implicitly_convertible<int, PrinterType>();

    auto potential_model{ py::enum_<PotentialModel>(module, "PotentialModel") };
    BindEnumEntries(potential_model);
    py::implicitly_convertible<int, PotentialModel>();

    auto partial_charge{ py::enum_<PartialCharge>(module, "PartialCharge") };
    BindEnumEntries(partial_charge);
    py::implicitly_convertible<int, PartialCharge>();

    auto tester_type{ py::enum_<TesterType>(module, "TesterType") };
    BindEnumEntries(tester_type);
    py::implicitly_convertible<int, TesterType>();
}

} // namespace rhbm_gem::bindings
