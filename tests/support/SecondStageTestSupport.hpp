#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

#include "core/detail/second_stage/SecondStageState.hpp"
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

namespace second_stage_test {

struct JointPolishFixture
{
    rhbm_gem::core::detail::SecondStageContext context{};
    rhbm_gem::core::detail::FitState state{};
    std::vector<rhbm_gem::core::detail::SampleRef>
        sample_ref_list{};
};

rhbm_gem::core::FitOptions MakeSecondStageOptions();

double Distance(
    const std::array<double, 3> & lhs,
    const std::array<double, 3> & rhs);

std::unique_ptr<rhbm_gem::AtomObject> MakeAtom(
    int serial_id,
    Spot spot,
    Element element,
    const std::array<double, 3> & position);

rhbm_gem::LocalGaussianResult MakeGaussianResult(const rhbm_gem::GaussianModel3D & model);

LocalPotentialSampleList BuildSamples(
    const rhbm_gem::AtomObject & target_atom,
    const std::vector<rhbm_gem::AtomObject *> & atom_list,
    const std::vector<rhbm_gem::GaussianModel3D> & truth_model_list);

JointPolishFixture BuildJointPolishFixture(
    const std::vector<rhbm_gem::GaussianModel3D> & base_model_list,
    const std::vector<rhbm_gem::GaussianModel3D> & target_model_list);

std::unique_ptr<rhbm_gem::ModelObject> BuildDefenseModel(
    const std::vector<std::array<double, 3>> & position_list,
    const std::vector<Spot> & spot_list,
    const std::vector<Element> & element_list,
    const std::vector<rhbm_gem::GaussianModel3D> & truth_model_list,
    const rhbm_gem::GaussianModel3D & initial_model);

std::unique_ptr<rhbm_gem::ModelObject> BuildNearCollinearDefenseModel(
    double intensity_scale = 1.0);

std::unique_ptr<rhbm_gem::ModelObject> BuildJointPolishDefenseModel();

std::unique_ptr<rhbm_gem::ModelObject> BuildIndependentOffsetDefenseModel(
    double intensity_scale = 1.0,
    bool alternate_keys = false);

std::unique_ptr<rhbm_gem::ModelObject> BuildSeparatedRollbackDefenseModel();

double CalculateSelectedAtomResponseMeanSquaredError(
    const rhbm_gem::ModelObject & model,
    std::size_t target_begin,
    std::size_t target_end);

double CalculateSelectedAtomResponseMeanSquaredError(const rhbm_gem::ModelObject & model);

rhbm_gem::GaussianModel3D GetEstimateModel(const rhbm_gem::AtomObject & atom);

void ExpectGaussianModelsNear(
    const rhbm_gem::GaussianModel3D & actual,
    const rhbm_gem::GaussianModel3D & expected,
    double tolerance);

void ExpectSelectedAtomEstimatesAreFinite(const rhbm_gem::ModelObject & model);

} // namespace second_stage_test
