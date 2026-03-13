#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace rhbm_gem::bindings {

void BindCommonTypes(py::module_ & module);
void BindCommandApi(py::module_ & module);

} // namespace rhbm_gem::bindings

PYBIND11_MODULE(rhbm_gem_module, module)
{
    rhbm_gem::bindings::BindCommonTypes(module);
    rhbm_gem::bindings::BindCommandApi(module);
}
