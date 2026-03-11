#include "BindingHelpers.hpp"

#include <rhbm_gem/core/command/PotentialDisplayCommand.hpp>

namespace rhbm_gem::bindings {

template <>
void BindCommand<PotentialDisplayCommand>(py::module_ & module)
{
    auto potential_display{ BindBuiltInCommand<PotentialDisplayCommand>(module) };
    potential_display
        .def(py::init<const DataIoServices &>())
        .def("Execute", &PotentialDisplayCommand::Execute)
        .def("SetPainterChoice", &PotentialDisplayCommand::SetPainterChoice)
        .def("SetModelKeyTagList", &PotentialDisplayCommand::SetModelKeyTagList)
        .def("SetRefModelKeyTagListMap", &PotentialDisplayCommand::SetRefModelKeyTagListMap)
        .def("SetPickChainID", &PotentialDisplayCommand::SetPickChainID)
        .def("SetPickResidueType", &PotentialDisplayCommand::SetPickResidueType)
        .def("SetPickElementType", &PotentialDisplayCommand::SetPickElementType)
        .def("SetVetoChainID", &PotentialDisplayCommand::SetVetoChainID)
        .def("SetVetoResidueType", &PotentialDisplayCommand::SetVetoResidueType)
        .def("SetVetoElementType", &PotentialDisplayCommand::SetVetoElementType);
    BindCommonCommandSetters(potential_display);
    BindCommandDiagnostics(potential_display);
}

} // namespace rhbm_gem::bindings
