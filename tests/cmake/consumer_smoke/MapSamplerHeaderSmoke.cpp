#include <rhbm_gem/core/MapSampler.hpp>

#include <string_view>
#include <type_traits>
#include <vector>

using DefaultAtomSampler = LocalPotentialSampleList (*)(
    const rhbm_gem::MapObject &, const rhbm_gem::AtomObject &);
using ExplicitAtomSampler = LocalPotentialSampleList (*)(
    const rhbm_gem::MapObject &, const rhbm_gem::AtomObject &,
    SphereSamplingMethod);
using AtomDomainWorkflow = void (*)(
    rhbm_gem::MapObject &, rhbm_gem::ModelObject &,
    const std::vector<rhbm_gem::AtomObject *> &, int, std::string_view, bool);
using SelectedAtomsWorkflow = void (*)(
    rhbm_gem::MapObject &, rhbm_gem::ModelObject &, int);

static_assert(std::is_same_v<decltype(static_cast<DefaultAtomSampler>(
    &rhbm_gem::core::SampleAtomMapValues)), DefaultAtomSampler>);
static_assert(std::is_same_v<decltype(static_cast<ExplicitAtomSampler>(
    &rhbm_gem::core::SampleAtomMapValues)), ExplicitAtomSampler>);
static_assert(std::is_same_v<decltype(static_cast<AtomDomainWorkflow>(
    &rhbm_gem::core::RunPotentialSamplingWorkflow)), AtomDomainWorkflow>);
static_assert(std::is_same_v<decltype(static_cast<SelectedAtomsWorkflow>(
    &rhbm_gem::core::RunPotentialSamplingWorkflow)), SelectedAtomsWorkflow>);
