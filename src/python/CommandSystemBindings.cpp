#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include <rhbm_gem/core/command/CommandSystem.hpp>

#include "command/detail/CommandCatalog.hpp"
#include "command/detail/CommandEnumCatalog.hpp"

#include <type_traits>

namespace py = pybind11;

namespace rhbm_gem::bindings {

namespace {

template <typename EnumType>
void BindEnumEntries(py::enum_<EnumType> & py_enum)
{
    for (const auto & option : command_internal::CommandEnumTraits<EnumType>::kOptions)
    {
        const std::string name{ option.binding_token };
        py_enum.value(name.c_str(), option.value);
    }
}

} // namespace

void BindCommandTypes(py::module_ & module)
{
    py::class_<CommandDiagnostic>(module, "CommandDiagnostic")
        .def_readonly("option_name", &CommandDiagnostic::option_name)
        .def_readonly("message", &CommandDiagnostic::message);

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

    auto sphere_sampling_profile_choice{
        py::enum_<SphereSamplingProfileChoice>(module, "SphereSamplingProfileChoice")
    };
    BindEnumEntries(sphere_sampling_profile_choice);
    py::implicitly_convertible<int, SphereSamplingProfileChoice>();
}

void BindCommandSystem(py::module_ & module)
{
    auto base_request{ py::class_<CommandRequestBase>(module, "CommandRequestBase") };
    base_request.def(py::init<>());
    command_internal::RequestFieldCatalog<CommandRequestBase>::Visit([&](const auto & field)
    {
        base_request.def_readwrite(field.field_name, field.member);
    });

    command_internal::VisitCommandCatalog([&](const auto & entry)
    {
        using RequestType = typename std::decay_t<decltype(entry)>::Request;
        auto py_request{
            py::class_<RequestType, CommandRequestBase>(module, entry.request_type_name.data())
        };
        py_request.def(py::init<>());
        command_internal::RequestFieldCatalog<RequestType>::Visit([&](const auto & field)
        {
            py_request.def_readwrite(field.field_name, field.member);
        });
    });

    py::class_<CommandResult>(module, "CommandResult")
        .def(py::init<>())
        .def_readonly("succeeded", &CommandResult::succeeded)
        .def_readonly("issues", &CommandResult::issues);

    command_internal::VisitCommandCatalog([&](const auto & entry)
    {
        using RequestType = typename std::decay_t<decltype(entry)>::Request;
        module.def(
            "RunCommand",
            [execute = entry.execute](const RequestType & request)
            {
                return execute(request);
            });
    });
}

} // namespace rhbm_gem::bindings

PYBIND11_MODULE(rhbm_gem_module, module)
{
    rhbm_gem::bindings::BindCommandTypes(module);
    rhbm_gem::bindings::BindCommandSystem(module);
}
