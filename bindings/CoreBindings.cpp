#include "BindingHelpers.hpp"

PYBIND11_MODULE(rhbm_gem_module, module)
{
    rhbm_gem::bindings::BindCommonTypes(module);
    rhbm_gem::bindings::BindPotentialAnalysis(module);
    rhbm_gem::bindings::BindPotentialDisplay(module);
    rhbm_gem::bindings::BindResultDump(module);
    rhbm_gem::bindings::BindMapSimulation(module);
    rhbm_gem::bindings::BindMapVisualization(module);
    rhbm_gem::bindings::BindPositionEstimation(module);
    rhbm_gem::bindings::BindHRLModelTest(module);
}
