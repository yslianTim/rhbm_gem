#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include <rhbm_gem/core/command/CommandSystem.hpp>

#include "command/detail/CommandRequestSchema.hpp"

#include <type_traits>

namespace py = pybind11;

namespace rhbm_gem::bindings {

namespace {

template <typename EnumType>
void BindEnumEntries(py::enum_<EnumType> & py_enum)
{
    for (const auto & option : internal::CommandEnumTraits<EnumType>::kOptions)
    {
        const std::string name{ option.binding_token };
        py_enum.value(name.c_str(), option.value);
    }
}

template <typename Request>
void BindRequestType(py::module_ & module, const char * type_name)
{
    auto py_request{ py::class_<Request, CommandRequestBase>(module, type_name) };
    py_request.def(py::init<>());
    command_internal::CommandRequestSchema<Request>::Visit([&](const auto & field)
    {
        py_request.def_readwrite(field.python_name, field.member);
    });
}

} // namespace

void BindCommonTypes(py::module_ & module)
{
    py::class_<ValidationIssue>(module, "ValidationIssue")
        .def_readonly("option_name", &ValidationIssue::option_name)
        .def_readonly("message", &ValidationIssue::message);

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

void BindCommandSystem(py::module_ & module)
{
    auto base_request{ py::class_<CommandRequestBase>(module, "CommandRequestBase") };
    base_request.def(py::init<>());
    command_internal::CommandRequestSchema<CommandRequestBase>::Visit([&](const auto & field)
    {
        base_request.def_readwrite(field.python_name, field.member);
    });

    command::VisitCommands([&](const auto & entry)
    {
        using RequestType = typename std::decay_t<decltype(entry)>::Request;
        BindRequestType<RequestType>(module, entry.request_type_name.data());
    });

    py::class_<CommandResult>(module, "CommandResult")
        .def(py::init<>())
        .def_readonly("succeeded", &CommandResult::succeeded)
        .def_readonly("issues", &CommandResult::issues);

    command::VisitCommands([&](const auto & entry)
    {
        using RequestType = typename std::decay_t<decltype(entry)>::Request;
        module.def(
            "RunCommand",
            static_cast<CommandResult (*)(const RequestType &)>(&RunCommand<RequestType>));
    });
}

} // namespace rhbm_gem::bindings

PYBIND11_MODULE(rhbm_gem_module, module)
{
    rhbm_gem::bindings::BindCommonTypes(module);
    rhbm_gem::bindings::BindCommandSystem(module);
}
