#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "core/detail/GaussianEstimatorStages.hpp"
#include "core/detail/LocalFittingAudit.hpp"
#include "core/detail/LocalFittingAndersonRegime.hpp"
#include "core/detail/LocalFittingHealth.hpp"
#include "core/detail/LocalFittingJointOffsetConditioning.hpp"
#include "core/detail/LocalFittingSeedRepair.hpp"
#include "core/detail/LocalFittingTransformedChange.hpp"
#include "core/detail/PostRefitRollback.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/algorithm/ClusteredAndersonAcceleration.hpp>
#include <rhbm_gem/utils/algorithm/ConvergenceFreezeTracker.hpp>
#include <rhbm_gem/utils/algorithm/DependencyThawHysteresisTracker.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

namespace {
namespace alg = rhbm_gem::algorithm;
namespace audit_detail = rhbm_gem::core::detail;
namespace change_detail = rhbm_gem::core::detail;
namespace conditioning_detail = rhbm_gem::core::detail;
namespace health_detail = rhbm_gem::core::detail;
namespace regime_detail = rhbm_gem::core::detail;
namespace rollback_detail = rhbm_gem::core::detail;
namespace seed_detail = rhbm_gem::core::detail;
namespace rt = rhbm_gem::core;
namespace rg = rhbm_gem;

rt::FitOptions MakeSecondStageOptions(bool quiet_mode = true)
{
    rt::FitOptions options;
    options.distance_min = 0.0;
    options.distance_max = 1.0;
    options.thread_size = 1;
    options.quiet_mode = quiet_mode;
    return options;
}

double Distance(
    const std::array<float, 3> & lhs,
    const std::array<float, 3> & rhs)
{
    const auto dx{ static_cast<double>(lhs.at(0) - rhs.at(0)) };
    const auto dy{ static_cast<double>(lhs.at(1) - rhs.at(1)) };
    const auto dz{ static_cast<double>(lhs.at(2) - rhs.at(2)) };
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

std::unique_ptr<rg::AtomObject> MakeAtom(
    int serial_id,
    Spot spot,
    Element element,
    const std::array<float, 3> & position)
{
    auto atom{ std::make_unique<rg::AtomObject>() };
    atom->SetSerialID(serial_id);
    atom->SetComponentKey(1);
    atom->SetAtomKey(static_cast<AtomKey>(spot));
    atom->SetElement(element);
    atom->SetSpot(spot);
    atom->SetPosition(position);
    return atom;
}

rg::LocalGaussianResult MakeGaussianResult(const rg::GaussianModel3D & model)
{
    rg::LocalGaussianResult result;
    result.ols = rg::GaussianModel3DWithUncertainty{
        model,
        rg::GaussianModel3DUncertainty{}
    };
    result.mdpde = rg::GaussianModel3DWithUncertainty{
        model,
        rg::GaussianModel3DUncertainty{}
    };
    return result;
}

LocalPotentialSampleList BuildSamples(
    const rg::AtomObject & target_atom,
    const std::vector<rg::AtomObject *> & atom_list,
    const std::vector<rg::GaussianModel3D> & truth_model_list)
{
    const std::array<std::array<float, 3>, 6> direction_list{
        std::array<float, 3>{ 1.0F, 0.0F, 0.0F },
        std::array<float, 3>{ -1.0F, 0.0F, 0.0F },
        std::array<float, 3>{ 0.0F, 1.0F, 0.0F },
        std::array<float, 3>{ 0.0F, -1.0F, 0.0F },
        std::array<float, 3>{ 0.0F, 0.0F, 1.0F },
        std::array<float, 3>{ 0.0F, 0.0F, -1.0F }
    };
    const std::array<float, 4> radius_list{ 0.15F, 0.35F, 0.65F, 0.95F };

    LocalPotentialSampleList sample_list;
    const auto target_position{ target_atom.GetPosition() };
    for (const auto radius : radius_list)
    {
        for (const auto & direction : direction_list)
        {
            SamplingPoint point;
            point.distance = radius;
            for (std::size_t axis = 0; axis < point.position.size(); axis++)
            {
                point.position.at(axis) =
                    target_position.at(axis) + radius * direction.at(axis);
            }

            double response{ 0.0 };
            for (std::size_t i = 0; i < atom_list.size(); i++)
            {
                response += truth_model_list.at(i).ResponseAtDistance(
                    Distance(point.position, atom_list.at(i)->GetPosition()));
            }
            sample_list.emplace_back(LocalPotentialSample{
                static_cast<float>(response),
                point
            });
        }
    }
    return sample_list;
}

std::unique_ptr<rg::ModelObject> BuildDefenseModel(
    const std::vector<std::array<float, 3>> & position_list,
    const std::vector<Spot> & spot_list,
    const std::vector<Element> & element_list,
    const std::vector<rg::GaussianModel3D> & truth_model_list,
    const rg::GaussianModel3D & initial_model)
{
    std::vector<std::unique_ptr<rg::AtomObject>> atom_list;
    for (std::size_t i = 0; i < position_list.size(); i++)
    {
        atom_list.emplace_back(MakeAtom(
            static_cast<int>(i + 1),
            spot_list.at(i),
            element_list.at(i),
            position_list.at(i)));
    }

    auto model{ std::make_unique<rg::ModelObject>(std::move(atom_list)) };
    model->SelectAllAtoms();
    auto analysis{ model->EditAnalysis() };
    const auto & selected_atoms{ model->GetSelectedAtoms() };
    for (auto * atom : selected_atoms)
    {
        auto local_editor{ analysis.EnsureAtomLocalPotential(*atom) };
        local_editor.SetAlphaR(0.0);
        local_editor.SetGaussianResult(MakeGaussianResult(initial_model));
        local_editor.SetSamplingEntries(
            BuildSamples(*atom, selected_atoms, truth_model_list));
    }
    return model;
}

std::unique_ptr<rg::ModelObject> BuildNearCollinearDefenseModel(
    double intensity_scale = 1.0)
{
    return BuildDefenseModel(
        {
            std::array<float, 3>{ 0.0F, 0.0F, 0.0F },
            std::array<float, 3>{ 1.0e-4F, 0.0F, 0.0F }
        },
        { Spot::C, Spot::O },
        { Element::CARBON, Element::OXYGEN },
        {
            rg::GaussianModel3D{ 6.0 * intensity_scale, 0.55, 0.20 * intensity_scale },
            rg::GaussianModel3D{ 5.5 * intensity_scale, 0.55, -0.15 * intensity_scale }
        },
        rg::GaussianModel3D{ 5.75 * intensity_scale, 0.55, 0.0 });
}

std::unique_ptr<rg::ModelObject> BuildInvalidCandidateDefenseModel()
{
    auto model{ BuildNearCollinearDefenseModel() };
    auto analysis{ model->EditAnalysis() };
    analysis.EnsureAtomLocalPotential(*model->GetSelectedAtoms().front())
        .SetGaussianResult(MakeGaussianResult(
            rg::GaussianModel3D{ 0.0, 0.0, 0.37 }));
    return model;
}

std::unique_ptr<rg::ModelObject> BuildAllInvalidSeedDefenseModel()
{
    auto model{ BuildNearCollinearDefenseModel() };
    auto analysis{ model->EditAnalysis() };
    const auto & atom_list{ model->GetSelectedAtoms() };
    analysis.EnsureAtomLocalPotential(*atom_list.at(0))
        .SetGaussianResult(MakeGaussianResult(
            rg::GaussianModel3D{ 0.0, 0.0, 0.25 }));
    analysis.EnsureAtomLocalPotential(*atom_list.at(1))
        .SetGaussianResult(MakeGaussianResult(
            rg::GaussianModel3D{ -1.0, -0.5, -0.15 }));
    return model;
}

void MakeAtomSamplesSuspicious(rg::ModelObject & model, std::size_t atom_index)
{
    const auto & selected_atoms{ model.GetSelectedAtoms() };
    auto * target_atom{ selected_atoms.at(atom_index) };
    const auto target_position{ target_atom->GetPosition() };

    auto sampling_entries{
        rg::AtomLocalPotentialView::RequireFor(*target_atom).GetSamplingEntries(false)
    };
    sampling_entries.resize(256);
    sampling_entries.front().response = 0.0F;
    sampling_entries.front().point.distance = 0.0F;
    sampling_entries.front().point.position = target_position;
    for (std::size_t i = 1; i < sampling_entries.size(); i++)
    {
        auto & sample{ sampling_entries.at(i) };
        const auto response_scale{
            0.5F + 0.5F * static_cast<float>(i) /
                static_cast<float>(sampling_entries.size())
        };
        sample.response = std::numeric_limits<float>::max() * response_scale;
        sample.point.position = target_position;
        sample.point.position.at(0) += 100.0F;
        sample.point.distance = 100.0F;
    }

    auto analysis{ model.EditAnalysis() };
    analysis.EnsureAtomLocalPotential(*target_atom).SetSamplingEntries(
        std::move(sampling_entries));
}

std::unique_ptr<rg::ModelObject> BuildSeparatedRollbackDefenseModel()
{
    auto model{
        BuildDefenseModel(
            {
                std::array<float, 3>{ 0.0F, 0.0F, 0.0F },
                std::array<float, 3>{ 1.0e-4F, 0.0F, 0.0F },
                std::array<float, 3>{ 10.0F, 0.0F, 0.0F },
                std::array<float, 3>{ 10.0001F, 0.0F, 0.0F }
            },
            { Spot::C, Spot::O, Spot::N, Spot::CA },
            { Element::CARBON, Element::OXYGEN, Element::NITROGEN, Element::CARBON },
            {
                rg::GaussianModel3D{ 6.0, 0.55, 0.20 },
                rg::GaussianModel3D{ 5.5, 0.55, -0.15 },
                rg::GaussianModel3D{ 6.2, 0.55, 0.18 },
                rg::GaussianModel3D{ 5.7, 0.55, -0.12 }
            },
            rg::GaussianModel3D{ 5.8, 0.55, 0.0 })
    };
    MakeAtomSamplesSuspicious(*model, 0);
    return model;
}

std::unique_ptr<rg::ModelObject> BuildSeparatedSystemBuildFailureDefenseModel()
{
    auto model{ BuildSeparatedRollbackDefenseModel() };
    auto * atom{ model->GetSelectedAtoms().front() };
    auto sampling_entries{
        rg::AtomLocalPotentialView::RequireFor(*atom).GetSamplingEntries(false)
    };
    sampling_entries.front().response = std::numeric_limits<float>::quiet_NaN();
    auto analysis{ model->EditAnalysis() };
    analysis.EnsureAtomLocalPotential(*atom).SetSamplingEntries(
        std::move(sampling_entries));
    return model;
}

std::unique_ptr<rg::ModelObject> BuildTerminalWithPersistentLocalRefitFallbackDefenseModel()
{
    auto model{ BuildSeparatedRollbackDefenseModel() };
    auto * atom{ model->GetSelectedAtoms().at(2) };
    auto sampling_entries{
        rg::AtomLocalPotentialView::RequireFor(*atom).GetSamplingEntries(false)
    };
    sampling_entries.resize(1);
    auto analysis{ model->EditAnalysis() };
    analysis.EnsureAtomLocalPotential(*atom).SetSamplingEntries(
        std::move(sampling_entries));
    return model;
}

std::unique_ptr<rg::ModelObject> BuildUnavailableAuditWithPersistentLocalRefitFallbackDefenseModel()
{
    auto model{ BuildSeparatedSystemBuildFailureDefenseModel() };
    auto * atom{ model->GetSelectedAtoms().at(2) };
    auto sampling_entries{
        rg::AtomLocalPotentialView::RequireFor(*atom).GetSamplingEntries(false)
    };
    sampling_entries.resize(1);
    auto analysis{ model->EditAnalysis() };
    analysis.EnsureAtomLocalPotential(*atom).SetSamplingEntries(
        std::move(sampling_entries));
    return model;
}

std::unique_ptr<rg::ModelObject> BuildSeparatedLocalRefitFallbackDefenseModel()
{
    auto model{
        BuildDefenseModel(
            {
                std::array<float, 3>{ 0.0F, 0.0F, 0.0F },
                std::array<float, 3>{ 1.0e-4F, 0.0F, 0.0F },
                std::array<float, 3>{ 10.0F, 0.0F, 0.0F },
                std::array<float, 3>{ 10.0001F, 0.0F, 0.0F }
            },
            { Spot::C, Spot::O, Spot::N, Spot::CA },
            { Element::CARBON, Element::OXYGEN, Element::NITROGEN, Element::CARBON },
            {
                rg::GaussianModel3D{ 6.0, 0.55, 0.20 },
                rg::GaussianModel3D{ 5.5, 0.55, -0.15 },
                rg::GaussianModel3D{ 6.2, 0.55, 0.18 },
                rg::GaussianModel3D{ 5.7, 0.55, -0.12 }
            },
            rg::GaussianModel3D{ 5.8, 0.55, 0.0 })
    };
    auto * atom{ model->GetSelectedAtoms().front() };
    auto sampling_entries{
        rg::AtomLocalPotentialView::RequireFor(*atom).GetSamplingEntries(false)
    };
    sampling_entries.resize(1);
    auto analysis{ model->EditAnalysis() };
    analysis.EnsureAtomLocalPotential(*atom).SetSamplingEntries(
        std::move(sampling_entries));
    return model;
}

std::unique_ptr<rg::ModelObject> BuildSeparatedEmptyJointOffsetDefenseModel()
{
    auto model{
        BuildDefenseModel(
            {
                std::array<float, 3>{ 0.0F, 0.0F, 0.0F },
                std::array<float, 3>{ 10.0F, 0.0F, 0.0F },
                std::array<float, 3>{ 10.0001F, 0.0F, 0.0F }
            },
            { Spot::C, Spot::N, Spot::CA },
            { Element::CARBON, Element::NITROGEN, Element::CARBON },
            {
                rg::GaussianModel3D{ 6.0, 0.55, 0.20 },
                rg::GaussianModel3D{ 6.2, 0.55, 0.18 },
                rg::GaussianModel3D{ 5.7, 0.55, -0.12 }
            },
            rg::GaussianModel3D{ 5.8, 0.55, 0.0 })
    };
    auto analysis{ model->EditAnalysis() };
    analysis.EnsureAtomLocalPotential(*model->GetSelectedAtoms().front())
        .SetSamplingEntries({});
    return model;
}

std::unique_ptr<rg::ModelObject> BuildPostRefitRollbackChainDefenseModel()
{
    auto model{
        BuildDefenseModel(
            {
                std::array<float, 3>{ 0.0F, 0.0F, 0.0F },
                std::array<float, 3>{ 2.0F, 0.0F, 0.0F },
                std::array<float, 3>{ 4.0F, 0.0F, 0.0F },
                std::array<float, 3>{ 6.0F, 0.0F, 0.0F },
                std::array<float, 3>{ 8.0F, 0.0F, 0.0F }
            },
            { Spot::C, Spot::O, Spot::N, Spot::CA, Spot::C },
            {
                Element::CARBON,
                Element::OXYGEN,
                Element::NITROGEN,
                Element::CARBON,
                Element::CARBON
            },
            {
                rg::GaussianModel3D{ 6.0, 0.55, 0.20 },
                rg::GaussianModel3D{ 5.5, 0.55, -0.15 },
                rg::GaussianModel3D{ 6.2, 0.55, 0.18 },
                rg::GaussianModel3D{ 5.7, 0.55, -0.12 },
                rg::GaussianModel3D{ 9.0, 0.55, 0.10 }
            },
            rg::GaussianModel3D{ 5.8, 0.55, 0.0 })
    };
    auto analysis{ model->EditAnalysis() };
    analysis.EnsureAtomLocalPotential(*model->GetSelectedAtoms().front())
        .SetGaussianResult(MakeGaussianResult(
            rg::GaussianModel3D{ 1.0e300, 0.55, 0.0 }));
    return model;
}

std::unique_ptr<rg::ModelObject> BuildNonFiniteJointOffsetDefenseModel()
{
    auto model{
        BuildDefenseModel(
            { std::array<float, 3>{ 0.0F, 0.0F, 0.0F } },
            { Spot::O },
            { Element::OXYGEN },
            { rg::GaussianModel3D{ 8.0, 0.5, -0.1 } },
            rg::GaussianModel3D{ 7.0, 0.5, 0.0 })
    };
    auto * atom{ model->GetSelectedAtoms().front() };
    auto sampling_entries{
        rg::AtomLocalPotentialView::RequireFor(*atom).GetSamplingEntries(false)
    };
    sampling_entries.front().response = std::numeric_limits<float>::quiet_NaN();
    auto analysis{ model->EditAnalysis() };
    analysis.EnsureAtomLocalPotential(*atom).SetSamplingEntries(std::move(sampling_entries));
    return model;
}

std::unique_ptr<rg::ModelObject> BuildEmptyJointOffsetDefenseModel()
{
    auto model{
        BuildDefenseModel(
            { std::array<float, 3>{ 0.0F, 0.0F, 0.0F } },
            { Spot::O },
            { Element::OXYGEN },
            { rg::GaussianModel3D{ 8.0, 0.5, -0.1 } },
            rg::GaussianModel3D{ 7.0, 0.5, 0.25 })
    };
    auto analysis{ model->EditAnalysis() };
    analysis.EnsureAtomLocalPotential(*model->GetSelectedAtoms().front()).SetSamplingEntries({});
    return model;
}

std::unique_ptr<rg::ModelObject> BuildFiniteNonphysicalProfileDefenseModel()
{
    const rg::GaussianModel3D initial_model{ 6.0, 0.55, 0.0 };
    auto model{
        BuildDefenseModel(
            { std::array<float, 3>{ 0.0F, 0.0F, 0.0F } },
            { Spot::O },
            { Element::OXYGEN },
            { initial_model },
            initial_model)
    };
    auto * atom{ model->GetSelectedAtoms().front() };
    auto sampling_entries{
        rg::AtomLocalPotentialView::RequireFor(*atom).GetSamplingEntries(false)
    };
    for (auto & sample : sampling_entries)
    {
        const auto distance{ static_cast<double>(sample.point.distance) };
        const auto outer_bias{ distance > 0.2 ? 12.0 : 8.0 };
        sample.response = static_cast<float>(
            initial_model.SignalAtDistance(distance) + outer_bias);
    }
    auto analysis{ model->EditAnalysis() };
    analysis.EnsureAtomLocalPotential(*atom).SetSamplingEntries(std::move(sampling_entries));
    return model;
}

double CalculateSelectedAtomResponseMeanSquaredError(
    const rg::ModelObject & model,
    std::size_t target_begin,
    std::size_t target_end)
{
    double squared_error_sum{ 0.0 };
    std::size_t sample_count{ 0 };
    const auto & selected_atoms{ model.GetSelectedAtoms() };
    for (std::size_t target_index = target_begin; target_index < target_end; target_index++)
    {
        const auto * atom{ selected_atoms.at(target_index) };
        const auto local_view{ rg::AtomLocalPotentialView::RequireFor(*atom) };
        for (const auto & sample : local_view.GetSamplingEntries(false))
        {
            double fitted_response{ 0.0 };
            for (const auto * fitted_atom : selected_atoms)
            {
                const auto fitted_view{
                    rg::AtomLocalPotentialView::RequireFor(*fitted_atom)
                };
                fitted_response += fitted_view.GetEstimateMDPDE().ResponseAtDistance(
                    Distance(sample.point.position, fitted_atom->GetPosition()));
            }
            const auto residual{ static_cast<double>(sample.response) - fitted_response };
            squared_error_sum += residual * residual;
            sample_count++;
        }
    }
    return squared_error_sum / static_cast<double>(sample_count);
}

double CalculateSelectedAtomResponseMeanSquaredError(const rg::ModelObject & model)
{
    return CalculateSelectedAtomResponseMeanSquaredError(
        model,
        0,
        model.GetSelectedAtomCount());
}

rg::GaussianModel3D GetEstimateModel(const rg::AtomObject & atom)
{
    return rg::AtomLocalPotentialView::RequireFor(atom).GetEstimateMDPDE();
}

void ExpectGaussianModelsNear(
    const rg::GaussianModel3D & actual,
    const rg::GaussianModel3D & expected,
    double tolerance)
{
    EXPECT_NEAR(actual.GetAmplitude(), expected.GetAmplitude(), tolerance);
    EXPECT_NEAR(actual.GetWidth(), expected.GetWidth(), tolerance);
    EXPECT_NEAR(actual.GetOffset(), expected.GetOffset(), tolerance);
}

void ExpectSelectedAtomEstimatesAreFinite(const rg::ModelObject & model)
{
    for (const auto * atom : model.GetSelectedAtoms())
    {
        const auto estimate{
            rg::AtomLocalPotentialView::RequireFor(*atom).GetEstimateMDPDE()
        };
        EXPECT_TRUE(std::isfinite(estimate.GetAmplitude()));
        EXPECT_TRUE(std::isfinite(estimate.GetWidth()));
        EXPECT_TRUE(std::isfinite(estimate.GetOffset()));
    }
}

std::vector<Eigen::VectorXd> MakeAndersonState(double left, double right)
{
    return {
        Eigen::VectorXd::Constant(1, left),
        Eigen::VectorXd::Constant(1, right)
    };
}

regime_detail::LocalFittingAndersonRegimeSignatureMap MakeAndersonRegimeSignatures(
    double global_ridge_ratio,
    double left_multiplier,
    double right_multiplier,
    health_detail::JointOffsetSolveStatus left_status =
        health_detail::JointOffsetSolveStatus::Converged,
    health_detail::JointOffsetSolveStatus right_status =
        health_detail::JointOffsetSolveStatus::Converged)
{
    return {
        { { 0 }, { left_status, global_ridge_ratio, { left_multiplier } } },
        { { 1 }, { right_status, global_ridge_ratio, { right_multiplier } } }
    };
}

std::size_t FindAcceptedCandidateCluster(
    const std::string & output,
    bool find_anderson)
{
    const std::string marker{ "accepted candidate clusters aa/fixed-point = " };
    std::size_t search_position{ 0 };
    while ((search_position = output.find(marker, search_position)) != std::string::npos)
    {
        const auto anderson_count_position{ search_position + marker.size() };
        const auto separator_position{ output.find('/', anderson_count_position) };
        if (separator_position == std::string::npos) return std::string::npos;
        const auto count_position{
            find_anderson ? anderson_count_position : separator_position + 1
        };
        if (count_position < output.size() && output.at(count_position) != '0')
        {
            return search_position;
        }
        search_position = separator_position + 1;
    }
    return std::string::npos;
}

} // namespace

TEST(EstimatorSecondStageDefenseTest, SeedRepairUsesConfiguredFallbackPriority)
{
    const auto make_candidate = [](double amplitude)
    {
        return rg::GaussianModel3DWithUncertainty{
            rg::GaussianModel3D{ amplitude, 0.5, -0.2 },
            rg::GaussianModel3DUncertainty{ 0.1, 0.02, 0.03 }
        };
    };
    const auto invalid_candidate{
        rg::GaussianModel3DWithUncertainty{
            rg::GaussianModel3D{ 0.0, 0.0, 0.0 },
            rg::GaussianModel3DUncertainty{}
        }
    };
    seed_detail::SecondStageSeedRepairCandidates candidates;
    candidates.group_posterior = make_candidate(1.0);
    candidates.group_prior = make_candidate(2.0);
    candidates.local_ols = make_candidate(3.0);
    candidates.group_median = make_candidate(4.0);
    candidates.global_median = make_candidate(5.0);

    const auto expect_source = [&](seed_detail::SecondStageSeedRepairSource source)
    {
        const auto selection{ seed_detail::SelectSecondStageSeedRepair(candidates) };
        ASSERT_TRUE(selection.has_value());
        EXPECT_EQ(selection->source, source);
    };
    expect_source(seed_detail::SecondStageSeedRepairSource::GroupPosterior);
    candidates.group_posterior = invalid_candidate;
    expect_source(seed_detail::SecondStageSeedRepairSource::GroupPrior);
    candidates.group_prior = invalid_candidate;
    expect_source(seed_detail::SecondStageSeedRepairSource::LocalOls);
    candidates.local_ols = invalid_candidate;
    expect_source(seed_detail::SecondStageSeedRepairSource::GroupMedian);
    candidates.group_median = invalid_candidate;
    expect_source(seed_detail::SecondStageSeedRepairSource::GlobalMedian);
    candidates.global_median = invalid_candidate;
    EXPECT_FALSE(seed_detail::SelectSecondStageSeedRepair(candidates).has_value());
}

TEST(EstimatorSecondStageDefenseTest, SeedRepairPreservesFiniteOffsetAndSourceUncertainty)
{
    const seed_detail::SecondStageSeedRepairSelection selection{
        seed_detail::SecondStageSeedRepairSource::GroupPosterior,
        rg::GaussianModel3DWithUncertainty{
            rg::GaussianModel3D{ 6.0, 0.55, -0.2 },
            rg::GaussianModel3DUncertainty{ 0.1, 0.02, 0.03 }
        }
    };

    const auto repaired_with_finite_offset{
        seed_detail::BuildRepairedSecondStageSeed(
            rg::GaussianModel3D{ 0.0, 0.0, 0.37 },
            selection)
    };
    EXPECT_DOUBLE_EQ(repaired_with_finite_offset.GetModel().GetAmplitude(), 6.0);
    EXPECT_DOUBLE_EQ(repaired_with_finite_offset.GetModel().GetWidth(), 0.55);
    EXPECT_DOUBLE_EQ(repaired_with_finite_offset.GetModel().GetOffset(), 0.37);
    EXPECT_DOUBLE_EQ(
        repaired_with_finite_offset.GetStandardDeviationModel().GetAmplitude(),
        0.1);
    EXPECT_DOUBLE_EQ(
        repaired_with_finite_offset.GetStandardDeviationModel().GetWidth(),
        0.02);
    EXPECT_DOUBLE_EQ(
        repaired_with_finite_offset.GetStandardDeviationModel().GetOffset(),
        0.03);

    const auto repaired_with_invalid_offset{
        seed_detail::BuildRepairedSecondStageSeed(
            rg::GaussianModel3D{
                0.0,
                0.0,
                std::numeric_limits<double>::quiet_NaN()
            },
            selection)
    };
    EXPECT_DOUBLE_EQ(repaired_with_invalid_offset.GetModel().GetOffset(), -0.2);
}

TEST(EstimatorSecondStageDefenseTest, AuditObjectiveKeepsEarlierBestOnTie)
{
    EXPECT_TRUE(audit_detail::IsBetterLocalFittingAuditObjective(
        0.8, 1.0, 1.0e-8));
    EXPECT_FALSE(audit_detail::IsBetterLocalFittingAuditObjective(
        1.2, 1.0, 1.0e-8));
    EXPECT_FALSE(audit_detail::IsBetterLocalFittingAuditObjective(
        1.0 - 0.5e-8, 1.0, 1.0e-8));
    EXPECT_FALSE(audit_detail::IsBetterLocalFittingAuditObjective(
        std::numeric_limits<double>::infinity(), 1.0, 1.0e-8));
    EXPECT_TRUE(audit_detail::IsBetterLocalFittingAuditObjective(
        1.0, std::numeric_limits<double>::infinity(), 1.0e-8));
    EXPECT_THROW(
        audit_detail::IsBetterLocalFittingAuditObjective(0.8, 1.0, -1.0),
        std::invalid_argument);
}

TEST(EstimatorSecondStageDefenseTest, ObjectivePenaltyUsesAtomMeanIndependentOfCardinality)
{
    const auto single_atom{
        audit_detail::BuildLocalFittingMeanObjectiveBreakdown(
            0.4,
            25.0,
            9.0,
            16.0,
            1,
            0.01,
            0.02,
            0.03)
    };
    const auto repeated_atoms{
        audit_detail::BuildLocalFittingMeanObjectiveBreakdown(
            0.4,
            2500.0,
            900.0,
            1600.0,
            100,
            0.01,
            0.02,
            0.03)
    };
    const auto single_previous{
        audit_detail::BuildLocalFittingMeanObjectiveBreakdown(
            1.4,
            0.0,
            0.0,
            0.0,
            1,
            0.01,
            0.02,
            0.03)
    };
    const auto repeated_previous{
        audit_detail::BuildLocalFittingMeanObjectiveBreakdown(
            1.4,
            0.0,
            0.0,
            0.0,
            100,
            0.01,
            0.02,
            0.03)
    };

    ASSERT_TRUE(single_atom.has_value());
    ASSERT_TRUE(repeated_atoms.has_value());
    ASSERT_TRUE(single_previous.has_value());
    ASSERT_TRUE(repeated_previous.has_value());
    EXPECT_DOUBLE_EQ(
        single_atom->residual_objective,
        repeated_atoms->residual_objective);
    EXPECT_DOUBLE_EQ(
        single_atom->width_prior_penalty,
        repeated_atoms->width_prior_penalty);
    EXPECT_DOUBLE_EQ(
        single_atom->offset_plausibility_penalty,
        repeated_atoms->offset_plausibility_penalty);
    EXPECT_DOUBLE_EQ(
        single_atom->movement_penalty,
        repeated_atoms->movement_penalty);
    EXPECT_DOUBLE_EQ(
        single_atom->total_objective,
        repeated_atoms->total_objective);
    EXPECT_DOUBLE_EQ(
        single_atom->total_objective,
        single_atom->residual_objective +
            single_atom->width_prior_penalty +
            single_atom->offset_plausibility_penalty +
            single_atom->movement_penalty);
    EXPECT_GT(single_atom->movement_penalty, 0.0);
    EXPECT_DOUBLE_EQ(single_previous->movement_penalty, 0.0);
    EXPECT_EQ(
        audit_detail::IsBetterLocalFittingAuditObjective(
            single_atom->total_objective,
            single_previous->total_objective,
            1.0e-8),
        audit_detail::IsBetterLocalFittingAuditObjective(
            repeated_atoms->total_objective,
            repeated_previous->total_objective,
            1.0e-8));
    EXPECT_FALSE(
        audit_detail::BuildLocalFittingMeanObjectiveBreakdown(
            0.4,
            25.0,
            9.0,
            16.0,
            0,
            0.01,
            0.02,
            0.03).has_value());
    EXPECT_FALSE(
        audit_detail::BuildLocalFittingMeanObjectiveBreakdown(
            std::numeric_limits<double>::infinity(),
            25.0,
            9.0,
            16.0,
            1,
            0.01,
            0.02,
            0.03).has_value());
}

TEST(EstimatorSecondStageDefenseTest, AndersonRegimeTracksGlobalAndClusterEffectiveRidge)
{
    const std::vector<regime_detail::LocalFittingAndersonRegimeClusterKey> key_list{
        { 0 }, { 1 }
    };
    const auto baseline{ MakeAndersonRegimeSignatures(1.0e-3, 1.0, 1.0) };
    regime_detail::LocalFittingAndersonRegimeTracker tracker;
    tracker.Commit(key_list, baseline);

    EXPECT_TRUE(tracker.FindIncompatible(baseline).empty());
    EXPECT_EQ(key_list, tracker.FindIncompatible(
        MakeAndersonRegimeSignatures(8.0e-4, 1.0, 1.0)));
    EXPECT_EQ(
        (std::vector<regime_detail::LocalFittingAndersonRegimeClusterKey>{ { 0 } }),
        tracker.FindIncompatible(
            MakeAndersonRegimeSignatures(1.0e-3, 10.0, 1.0)));
}

TEST(EstimatorSecondStageDefenseTest, AndersonRegimeDetectsRetryAndCollinearityTransitions)
{
    const std::vector<regime_detail::LocalFittingAndersonRegimeClusterKey> key_list{
        { 0 }, { 1 }
    };
    const auto baseline{ MakeAndersonRegimeSignatures(1.0e-3, 1.0, 1.0) };
    const auto guarded{ MakeAndersonRegimeSignatures(1.0e-3, 10.0, 1.0) };
    regime_detail::LocalFittingAndersonRegimeTracker tracker;
    tracker.Commit(key_list, baseline);

    EXPECT_EQ(
        (std::vector<regime_detail::LocalFittingAndersonRegimeClusterKey>{ { 0 } }),
        tracker.FindIncompatible(guarded));
    tracker.Commit({ { 0 } }, guarded);
    EXPECT_EQ(
        (std::vector<regime_detail::LocalFittingAndersonRegimeClusterKey>{ { 0 } }),
        tracker.FindIncompatible(baseline));

    tracker.Reconcile({ { 0, 1 } });
    const regime_detail::LocalFittingAndersonRegimeSignatureMap merged{
        { { 0, 1 }, {
            health_detail::JointOffsetSolveStatus::Converged,
            1.0e-3,
            { 1.0, 1.0 }
        } }
    };
    EXPECT_TRUE(tracker.FindIncompatible(merged).empty());
}

TEST(EstimatorSecondStageDefenseTest, AndersonRegimeRejectsInvalidSignatures)
{
    regime_detail::LocalFittingAndersonRegimeTracker tracker;
    using Status = health_detail::JointOffsetSolveStatus;
    const regime_detail::LocalFittingAndersonRegimeSignatureMap invalid_ridge{
        { { 0 }, { Status::Converged, 0.0, { 1.0 } } }
    };
    const regime_detail::LocalFittingAndersonRegimeSignatureMap invalid_multiplier{
        { { 0 }, {
            Status::Converged,
            1.0e-3,
            { std::numeric_limits<double>::infinity() }
        } }
    };
    const regime_detail::LocalFittingAndersonRegimeSignatureMap wrong_size_signature{
        { { 0 }, { Status::Converged, 1.0e-3, { 1.0, 1.0 } } }
    };
    const regime_detail::LocalFittingAndersonRegimeSignatureMap invalid_status{
        { { 0 }, { static_cast<Status>(-1), 1.0e-3, { 1.0 } } }
    };
    EXPECT_THROW(tracker.FindIncompatible(invalid_ridge), std::invalid_argument);
    EXPECT_THROW(tracker.FindIncompatible(invalid_multiplier), std::invalid_argument);
    EXPECT_THROW(tracker.FindIncompatible(wrong_size_signature), std::invalid_argument);
    EXPECT_THROW(tracker.FindIncompatible(invalid_status), std::logic_error);
}

TEST(EstimatorSecondStageDefenseTest, AndersonRegimeTracksJointOffsetStatusPerCluster)
{
    using Status = health_detail::JointOffsetSolveStatus;
    const std::vector<regime_detail::LocalFittingAndersonRegimeClusterKey> key_list{
        { 0 }, { 1 }
    };
    const auto baseline{ MakeAndersonRegimeSignatures(1.0e-3, 1.0, 1.0) };
    const auto changed{
        MakeAndersonRegimeSignatures(
            1.0e-3,
            1.0,
            1.0,
            Status::IrlsMaximumIterationsReached,
            Status::Converged)
    };
    regime_detail::LocalFittingAndersonRegimeTracker tracker;
    alg::ClusteredAndersonAccelerationHistorySet history{
        alg::AndersonAccelerationOptions{ 5, 100.0, 10.0, 1.0e-12 }
    };
    history.Reconcile(key_list);
    tracker.Commit(key_list, baseline);
    history.Commit(
        key_list,
        MakeAndersonState(0.0, 0.0),
        MakeAndersonState(1.0, 1.0));

    const auto incompatible_key_list{ tracker.FindIncompatible(changed) };
    EXPECT_EQ(
        (std::vector<regime_detail::LocalFittingAndersonRegimeClusterKey>{ { 0 } }),
        incompatible_key_list);
    history.ClearAndSuppress(incompatible_key_list);
    tracker.Invalidate(incompatible_key_list);
    const auto candidate{
        history.BuildCandidate(
            key_list,
            MakeAndersonState(1.0, 1.0),
            MakeAndersonState(1.5, 1.5))
    };

    ASSERT_TRUE(candidate.has_value());
    EXPECT_EQ((std::vector<alg::ClusterKey>{ { 1 } }), candidate->used_cluster_key_list);
}

TEST(EstimatorSecondStageDefenseTest, AndersonRegimeInvalidatesOnlyAffectedHistory)
{
    const std::vector<alg::ClusterKey> key_list{ { 0 }, { 1 } };
    alg::ClusteredAndersonAccelerationHistorySet history{
        alg::AndersonAccelerationOptions{ 5, 100.0, 10.0, 1.0e-12 }
    };
    regime_detail::LocalFittingAndersonRegimeTracker tracker;
    const auto baseline{ MakeAndersonRegimeSignatures(1.0e-3, 1.0, 1.0) };
    const auto changed{ MakeAndersonRegimeSignatures(1.0e-3, 10.0, 1.0) };
    history.Reconcile(key_list);
    tracker.Commit(key_list, baseline);
    history.Commit(key_list, MakeAndersonState(0.0, 0.0), MakeAndersonState(1.0, 1.0));

    const auto incompatible_key_list{ tracker.FindIncompatible(changed) };
    history.ClearAndSuppress(incompatible_key_list);
    tracker.Invalidate(incompatible_key_list);
    const auto candidate{
        history.BuildCandidate(
            key_list,
            MakeAndersonState(1.0, 1.0),
            MakeAndersonState(1.5, 1.5))
    };

    ASSERT_TRUE(candidate.has_value());
    EXPECT_EQ((std::vector<alg::ClusterKey>{ { 1 } }), candidate->used_cluster_key_list);

    history.ReleaseSuppression({ { 0 } });
    const auto before_recommit{
        history.BuildCandidate(
            key_list,
            MakeAndersonState(1.0, 1.0),
            MakeAndersonState(1.5, 1.5))
    };
    ASSERT_TRUE(before_recommit.has_value());
    EXPECT_EQ((std::vector<alg::ClusterKey>{ { 1 } }), before_recommit->used_cluster_key_list);

    history.Commit({ { 0 } }, MakeAndersonState(1.0, 1.0), MakeAndersonState(1.5, 1.5));
    tracker.Commit({ { 0 } }, changed);
    const auto after_recommit{
        history.BuildCandidate(
            key_list,
            MakeAndersonState(1.5, 1.5),
            MakeAndersonState(1.75, 1.75))
    };
    ASSERT_TRUE(after_recommit.has_value());
    EXPECT_NE(
        std::find(
            after_recommit->used_cluster_key_list.begin(),
            after_recommit->used_cluster_key_list.end(),
            alg::ClusterKey{ 0 }),
        after_recommit->used_cluster_key_list.end());
}

TEST(EstimatorSecondStageDefenseTest, JointOffsetConditioningDetectsJointDependence)
{
    Eigen::SparseMatrix<double> design_matrix{ 3, 3 };
    const std::vector<Eigen::Triplet<double>> entries{
        { 0, 0, 1.0 }, { 0, 2, 1.0 },
        { 1, 1, 1.0 }, { 1, 2, 1.0 },
        { 2, 0, 1.0 }, { 2, 1, 1.0 }, { 2, 2, 2.0 }
    };
    design_matrix.setFromTriplets(entries.begin(), entries.end());

    const Eigen::MatrixXd dense{ design_matrix };
    for (Eigen::Index left = 0; left < dense.cols(); left++)
    {
        for (Eigen::Index right = left + 1; right < dense.cols(); right++)
        {
            const auto overlap{
                std::abs(dense.col(left).dot(dense.col(right))) /
                (dense.col(left).norm() * dense.col(right).norm())
            };
            EXPECT_LT(overlap, 0.98);
        }
    }

    const auto diagnostics{
        conditioning_detail::EvaluateLocalFittingJointOffsetConditioning(
            design_matrix,
            1.0e-8)
    };
    EXPECT_TRUE(diagnostics.guard_required);
    EXPECT_LE(diagnostics.pivot_ratio, 1.0e-8);
}

TEST(EstimatorSecondStageDefenseTest, JointOffsetConditioningKeepsIndependentColumns)
{
    Eigen::SparseMatrix<double> design_matrix{ 3, 3 };
    design_matrix.setIdentity();

    const auto diagnostics{
        conditioning_detail::EvaluateLocalFittingJointOffsetConditioning(
            design_matrix,
            1.0e-8)
    };
    EXPECT_FALSE(diagnostics.guard_required);
    EXPECT_NEAR(diagnostics.pivot_ratio, 1.0, 1.0e-12);
}

TEST(EstimatorSecondStageDefenseTest, LocalRefitHealthSeparatesProgressFromStationarity)
{
    EXPECT_TRUE(health_detail::IsLocalGaussianRefitStatusProgressEligible(
        rg::RHBMEstimationStatus::SUCCESS));
    EXPECT_TRUE(health_detail::IsLocalGaussianRefitStatusStationarityEligible(
        rg::RHBMEstimationStatus::SUCCESS));
    EXPECT_TRUE(health_detail::IsLocalGaussianRefitStatusProgressEligible(
        rg::RHBMEstimationStatus::MAX_ITERATIONS_REACHED));
    EXPECT_FALSE(health_detail::IsLocalGaussianRefitStatusStationarityEligible(
        rg::RHBMEstimationStatus::MAX_ITERATIONS_REACHED));
    for (const auto status : {
        rg::RHBMEstimationStatus::NUMERICAL_FALLBACK,
        rg::RHBMEstimationStatus::INSUFFICIENT_DATA,
        rg::RHBMEstimationStatus::SINGLE_MEMBER })
    {
        EXPECT_FALSE(health_detail::IsLocalGaussianRefitStatusProgressEligible(status));
        EXPECT_FALSE(health_detail::IsLocalGaussianRefitStatusStationarityEligible(status));
    }
    EXPECT_THROW(
        health_detail::IsLocalGaussianRefitStatusProgressEligible(
            static_cast<rg::RHBMEstimationStatus>(-1)),
        std::logic_error);
    EXPECT_THROW(
        health_detail::IsLocalGaussianRefitStatusStationarityEligible(
            static_cast<rg::RHBMEstimationStatus>(-1)),
        std::logic_error);
}

TEST(EstimatorSecondStageDefenseTest, JointOffsetHealthSeparatesProgressFromStationarity)
{
    using Status = health_detail::JointOffsetSolveStatus;

    EXPECT_TRUE(health_detail::IsJointOffsetSolveProgressEligible(Status::Converged));
    EXPECT_TRUE(health_detail::IsJointOffsetSolveStationarityEligible(Status::Converged));
    EXPECT_FALSE(health_detail::IsJointOffsetSolveHardFailure(Status::Converged));

    for (const auto status : {
        Status::IrlsObjectiveDeteriorated,
        Status::IrlsMaximumIterationsReached })
    {
        EXPECT_TRUE(health_detail::IsJointOffsetSolveProgressEligible(status));
        EXPECT_FALSE(health_detail::IsJointOffsetSolveStationarityEligible(status));
        EXPECT_FALSE(health_detail::IsJointOffsetSolveHardFailure(status));
    }

    for (const auto status : {
        Status::SystemBuildFailed,
        Status::EmptySystem,
        Status::InitialSolveFailed,
        Status::IrlsSolveFailed })
    {
        EXPECT_FALSE(health_detail::IsJointOffsetSolveProgressEligible(status));
        EXPECT_FALSE(health_detail::IsJointOffsetSolveStationarityEligible(status));
        EXPECT_TRUE(health_detail::IsJointOffsetSolveHardFailure(status));
    }

    const auto invalid_status{ static_cast<Status>(-1) };
    EXPECT_THROW(
        health_detail::IsJointOffsetSolveProgressEligible(invalid_status),
        std::logic_error);
    EXPECT_FALSE(
        health_detail::IsJointOffsetSolveStationarityEligible(invalid_status));
    EXPECT_THROW(
        health_detail::IsJointOffsetSolveHardFailure(invalid_status),
        std::logic_error);
}

TEST(EstimatorSecondStageDefenseTest, TransformedChangeIsIntensityScaleInvariant)
{
    const rg::GaussianModel3D previous{ 8.0, 0.50, -0.10 };
    const rg::GaussianModel3D current{ 8.8, 0.55, -0.12 };
    const auto base_change{
        change_detail::CalculateLocalFittingTransformedChange(current, previous)
    };

    for (const auto scale : { 1.0e-2, 1.0e2 })
    {
        const auto scaled_change{
            change_detail::CalculateLocalFittingTransformedChange(
                rg::GaussianModel3D{
                    current.GetAmplitude() * scale,
                    current.GetWidth(),
                    current.GetOffset() * scale
                },
                rg::GaussianModel3D{
                    previous.GetAmplitude() * scale,
                    previous.GetWidth(),
                    previous.GetOffset() * scale
                })
        };
        ASSERT_EQ(base_change.value_list.size(), scaled_change.value_list.size());
        for (std::size_t i = 0; i < base_change.value_list.size(); i++)
        {
            EXPECT_NEAR(
                base_change.value_list.at(i),
                scaled_change.value_list.at(i),
                1.0e-12);
        }

        alg::DependencyThawHysteresisTracker thaw_tracker{ 1, 2.0, 8.0, 0.9, 5 };
        EXPECT_EQ(
            thaw_tracker.ShouldThaw(
                0,
                alg::GetMaximumParameterChange(base_change),
                1.0e-3),
            thaw_tracker.ShouldThaw(
                0,
                alg::GetMaximumParameterChange(scaled_change),
                1.0e-3));
    }
}

TEST(EstimatorSecondStageDefenseTest, TransformedCoordinatesRoundTrip)
{
    const rg::GaussianModel3D model{ 8.5, 0.65, -0.12 };
    const auto encoded{
        change_detail::EncodeLocalFittingTransformedCoordinates(model)
    };
    ASSERT_TRUE(encoded.has_value());

    const auto decoded{
        change_detail::DecodeLocalFittingTransformedCoordinates(*encoded)
    };
    ASSERT_TRUE(decoded.has_value());
    ExpectGaussianModelsNear(model, *decoded, 1.0e-12);
}

TEST(EstimatorSecondStageDefenseTest, TransformedDampingIsIntensityScaleInvariant)
{
    const rg::GaussianModel3D previous{ 8.0, 0.50, -0.10 };
    const rg::GaussianModel3D current{ 9.0, 0.60, -0.15 };
    constexpr double damping{ 0.25 };

    const auto damp = [&](const rg::GaussianModel3D & lhs,
                          const rg::GaussianModel3D & rhs)
    {
        const auto lhs_coordinates{
            change_detail::EncodeLocalFittingTransformedCoordinates(lhs)
        };
        const auto rhs_coordinates{
            change_detail::EncodeLocalFittingTransformedCoordinates(rhs)
        };
        EXPECT_TRUE(lhs_coordinates.has_value());
        EXPECT_TRUE(rhs_coordinates.has_value());
        return change_detail::DecodeLocalFittingTransformedCoordinates(
            *lhs_coordinates + damping * (*rhs_coordinates - *lhs_coordinates));
    };

    const auto base{ damp(previous, current) };
    ASSERT_TRUE(base.has_value());
    for (const auto scale : { 1.0e-2, 1.0e2 })
    {
        const auto scaled{
            damp(
                rg::GaussianModel3D{
                    previous.GetAmplitude() * scale,
                    previous.GetWidth(),
                    previous.GetOffset() * scale
                },
                rg::GaussianModel3D{
                    current.GetAmplitude() * scale,
                    current.GetWidth(),
                    current.GetOffset() * scale
                })
        };
        ASSERT_TRUE(scaled.has_value());
        EXPECT_NEAR(base->GetAmplitude() * scale, scaled->GetAmplitude(), 1.0e-10);
        EXPECT_NEAR(base->GetWidth(), scaled->GetWidth(), 1.0e-12);
        EXPECT_NEAR(base->GetOffset() * scale, scaled->GetOffset(), 1.0e-12);
    }
}

TEST(EstimatorSecondStageDefenseTest, TransformedExtrapolationKeepsPositiveShape)
{
    const auto left{
        change_detail::EncodeLocalFittingTransformedCoordinates(
            rg::GaussianModel3D{ 8.0, 0.50, -0.10 })
    };
    const auto right{
        change_detail::EncodeLocalFittingTransformedCoordinates(
            rg::GaussianModel3D{ 9.0, 0.60, -0.15 })
    };
    ASSERT_TRUE(left.has_value());
    ASSERT_TRUE(right.has_value());

    const auto extrapolated{
        change_detail::DecodeLocalFittingTransformedCoordinates(
            2.0 * *right - *left)
    };
    ASSERT_TRUE(extrapolated.has_value());
    EXPECT_GT(extrapolated->GetAmplitude(), 0.0);
    EXPECT_GT(extrapolated->GetWidth(), 0.0);
    EXPECT_TRUE(std::isfinite(extrapolated->GetOffset()));
}

TEST(EstimatorSecondStageDefenseTest, TransformedChangeSeparatesPeakHeightAndWidth)
{
    const auto change{
        change_detail::CalculateLocalFittingTransformedChange(
            rg::GaussianModel3D{ 8.0, 1.0, 0.0 },
            rg::GaussianModel3D{ 1.0, 0.5, 0.0 })
    };

    EXPECT_NEAR(
        0.0,
        change.value_list.at(change_detail::kLogPeakHeightChangeIndex),
        1.0e-12);
    EXPECT_NEAR(
        std::log(2.0),
        change.value_list.at(change_detail::kLogWidthChangeIndex),
        1.0e-12);
    EXPECT_DOUBLE_EQ(
        0.0,
        change.value_list.at(change_detail::kOffsetToPeakRatioChangeIndex));
}

TEST(EstimatorSecondStageDefenseTest, InvalidTransformedCoordinatesProduceInfiniteChange)
{
    const rg::GaussianModel3D valid_model{ 1.0, 0.5, 0.0 };
    const std::array<rg::GaussianModel3D, 4> invalid_model_list{
        rg::GaussianModel3D{ 0.0, 0.5, 0.0 },
        rg::GaussianModel3D{ 1.0, 0.0, 0.0 },
        rg::GaussianModel3D{
            std::numeric_limits<double>::quiet_NaN(),
            0.5,
            0.0
        },
        rg::GaussianModel3D{
            std::numeric_limits<double>::denorm_min(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max()
        }
    };

    for (const auto & invalid_model : invalid_model_list)
    {
        EXPECT_FALSE(seed_detail::IsValidSecondStageGaussianModel(invalid_model));
        const auto change{
            change_detail::CalculateLocalFittingTransformedChange(
                invalid_model,
                valid_model)
        };
        ASSERT_EQ(change_detail::kTransformedChangeSize, change.value_list.size());
        for (const auto value : change.value_list)
        {
            EXPECT_TRUE(std::isinf(value));
        }
    }
}

TEST(EstimatorSecondStageDefenseTest, TransformedConvergenceRejectsHiddenMaximumTail)
{
    std::vector<alg::ParameterChange> change_list(
        1000,
        alg::ParameterChange{ std::vector<double>(3, 0.0) });
    change_list.back().value_list.at(change_detail::kLogPeakHeightChangeIndex) =
        2.0e-3;
    std::vector<std::size_t> index_list(change_list.size());
    for (std::size_t i = 0; i < index_list.size(); i++)
    {
        index_list.at(i) = i;
    }

    const auto percentile_stats{
        alg::SummarizeParameterChangeStats(change_list, index_list, 0.99)
    };
    const auto maximum_list{
        change_detail::SummarizeLocalFittingMaximumTransformedChanges(
            change_list,
            index_list)
    };
    EXPECT_LT(
        percentile_stats.percentile_list.at(
            change_detail::kLogPeakHeightChangeIndex),
        1.0e-4);
    EXPECT_FALSE(change_detail::IsLocalFittingTransformedChangeConverged(
        percentile_stats,
        maximum_list,
        1.0e-4,
        1.0e-3));
}

TEST(EstimatorSecondStageDefenseTest, FreezeUsesTransformedRatherThanAbsoluteMovement)
{
    const auto scale_consistent_small_change{
        change_detail::CalculateLocalFittingTransformedChange(
            rg::GaussianModel3D{ 1.0e8 * std::exp(5.0e-5), 0.5, 0.0 },
            rg::GaussianModel3D{ 1.0e8, 0.5, 0.0 })
    };
    const auto relatively_large_tiny_absolute_change{
        change_detail::CalculateLocalFittingTransformedChange(
            rg::GaussianModel3D{ 2.0e-8, 0.5, 0.0 },
            rg::GaussianModel3D{ 1.0e-8, 0.5, 0.0 })
    };
    alg::ConvergenceFreezeTracker stable_tracker{ 1, 1.0e-6, 0.1, 3 };
    alg::ConvergenceFreezeTracker moving_tracker{ 1, 1.0e-6, 0.1, 3 };

    for (int i = 0; i < 3; i++)
    {
        stable_tracker.Update({ scale_consistent_small_change }, { 0 });
        moving_tracker.Update({ relatively_large_tiny_absolute_change }, { 0 });
    }

    EXPECT_TRUE(stable_tracker.IsFrozen(0));
    EXPECT_FALSE(moving_tracker.IsFrozen(0));
}

TEST(EstimatorSecondStageDefenseTest, FreezeEvidenceRequiresSmallRawFixedPointResidual)
{
    const rg::GaussianModel3D previous{ 1.0e8, 0.5, 0.0 };
    const rg::GaussianModel3D accepted{
        previous.GetAmplitude() * std::exp(5.0e-5),
        previous.GetWidth(),
        previous.GetOffset()
    };
    const rg::GaussianModel3D raw_fixed_point{
        previous.GetAmplitude() * std::exp(1.0e-2),
        previous.GetWidth(),
        previous.GetOffset()
    };
    const auto unstable_evidence{
        change_detail::CalculateLocalFittingFreezeEvidenceChange(
            accepted,
            raw_fixed_point,
            previous)
    };
    const auto stable_evidence{
        change_detail::CalculateLocalFittingFreezeEvidenceChange(
            accepted,
            accepted,
            previous)
    };
    alg::ConvergenceFreezeTracker unstable_tracker{ 1, 1.0e-6, 0.1, 3 };
    alg::ConvergenceFreezeTracker stable_tracker{ 1, 1.0e-6, 0.1, 3 };

    for (int i = 0; i < 3; i++)
    {
        unstable_tracker.Update({ unstable_evidence }, { 0 });
        stable_tracker.Update({ stable_evidence }, { 0 });
    }

    EXPECT_FALSE(unstable_tracker.IsFrozen(0));
    EXPECT_TRUE(stable_tracker.IsFrozen(0));
    EXPECT_NEAR(
        1.0e-2,
        unstable_evidence.value_list.at(
            change_detail::kLogPeakHeightChangeIndex),
        1.0e-12);
}

TEST(EstimatorSecondStageDefenseTest, FreezeEvidenceIsScaleInvariantAndRejectsInvalidRawState)
{
    const rg::GaussianModel3D previous{ 8.0, 0.50, -0.10 };
    const rg::GaussianModel3D accepted{ 8.1, 0.51, -0.11 };
    const rg::GaussianModel3D raw_fixed_point{ 8.8, 0.55, -0.12 };
    const auto base_evidence{
        change_detail::CalculateLocalFittingFreezeEvidenceChange(
            accepted,
            raw_fixed_point,
            previous)
    };
    const auto scale{ 1.0e2 };
    const auto scaled_evidence{
        change_detail::CalculateLocalFittingFreezeEvidenceChange(
            rg::GaussianModel3D{
                accepted.GetAmplitude() * scale,
                accepted.GetWidth(),
                accepted.GetOffset() * scale
            },
            rg::GaussianModel3D{
                raw_fixed_point.GetAmplitude() * scale,
                raw_fixed_point.GetWidth(),
                raw_fixed_point.GetOffset() * scale
            },
            rg::GaussianModel3D{
                previous.GetAmplitude() * scale,
                previous.GetWidth(),
                previous.GetOffset() * scale
            })
    };
    ASSERT_EQ(base_evidence.value_list.size(), scaled_evidence.value_list.size());
    for (std::size_t i = 0; i < base_evidence.value_list.size(); i++)
    {
        EXPECT_NEAR(
            base_evidence.value_list.at(i),
            scaled_evidence.value_list.at(i),
            1.0e-12);
    }

    const auto invalid_evidence{
        change_detail::CalculateLocalFittingFreezeEvidenceChange(
            accepted,
            rg::GaussianModel3D{ 0.0, 0.5, 0.0 },
            previous)
    };
    for (const auto value : invalid_evidence.value_list)
    {
        EXPECT_TRUE(std::isinf(value));
    }
}

TEST(EstimatorSecondStageDefenseTest, PostRefitRollbackExpandsCompleteContributorCluster)
{
    const std::vector<std::size_t> active_index_list{ 10, 11, 12, 13, 20 };
    const std::vector<rollback_detail::PostRefitRollbackClusterKey> cluster_key_list{
        { 10, 11, 12, 13 },
        { 20 }
    };
    std::vector<char> suspicious_mask{ 0, 1, 0, 0, 0 };

    const auto affected_position_list{
        rollback_detail::ExpandPostRefitRollbackClusters(
            active_index_list,
            cluster_key_list,
            { 3, 0, 3 },
            suspicious_mask)
    };

    EXPECT_EQ(
        (std::vector<std::size_t>{ 0, 1, 2, 3 }),
        affected_position_list);
    EXPECT_EQ((std::vector<char>{ 1, 1, 1, 1, 0 }), suspicious_mask);

    std::vector<char> permuted_suspicious_mask(5, 0);
    EXPECT_EQ(
        (std::vector<std::size_t>{ 0, 1, 2, 3 }),
        rollback_detail::ExpandPostRefitRollbackClusters(
            { 13, 12, 11, 10, 20 },
            cluster_key_list,
            { 0 },
            permuted_suspicious_mask));
    EXPECT_EQ(
        (std::vector<char>{ 1, 1, 1, 1, 0 }),
        permuted_suspicious_mask);
}

TEST(EstimatorSecondStageDefenseTest, PostRefitRollbackRejectsInconsistentTopology)
{
    std::vector<char> suspicious_mask{ 0, 0 };
    EXPECT_THROW(
        rollback_detail::ExpandPostRefitRollbackClusters(
            { 10, 11 },
            { { 10 } },
            { 0 },
            suspicious_mask),
        std::invalid_argument);
    EXPECT_THROW(
        rollback_detail::ExpandPostRefitRollbackClusters(
            { 10, 11 },
            { { 10, 11 } },
            { 2 },
            suspicious_mask),
        std::invalid_argument);
    EXPECT_THROW(
        rollback_detail::ExpandPostRefitRollbackClusters(
            { 10, 11 },
            { { 10 }, { 10, 11 } },
            { 0 },
            suspicious_mask),
        std::invalid_argument);
}

TEST(EstimatorSecondStageDefenseTest, PostRefitRollbackRestoresCompleteLongChain)
{
    auto model{ BuildPostRefitRollbackChainDefenseModel() };
    const auto & selected_atoms{ model->GetSelectedAtoms() };
    std::vector<rg::GaussianModel3D> previous_model_list;
    previous_model_list.reserve(selected_atoms.size());
    for (const auto * atom : selected_atoms)
    {
        previous_model_list.emplace_back(GetEstimateModel(*atom));
    }

    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions(false));
    const std::string output{ testing::internal::GetCapturedStdout() };
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(previous_log_level);

    for (std::size_t i = 0; i < selected_atoms.size(); i++)
    {
        ExpectGaussianModelsNear(
            GetEstimateModel(*selected_atoms.at(i)),
            previous_model_list.at(i),
            1.0e-12);
    }
    EXPECT_NE(output.find("terminal-suspicious atoms = 5"), std::string::npos);
    EXPECT_NE(
        error_output.find(
            "terminal suspicious rollback fallback clusters/atoms = 1/5"),
        std::string::npos);
    EXPECT_EQ(FindAcceptedCandidateCluster(output, true), std::string::npos);
    EXPECT_EQ(output.find("Iter. 6/200"), std::string::npos);
    EXPECT_EQ(error_output.find("Reached maximum iteration size"), std::string::npos);
}

TEST(EstimatorSecondStageDefenseTest, RunSecondStageLocalFittingFallsBackWhenJointOffsetSamplesAreNonFinite)
{
    auto model{ BuildNonFiniteJointOffsetDefenseModel() };
    auto * atom{ model->GetSelectedAtoms().front() };
    const auto previous_model{ GetEstimateModel(*atom) };

    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions(false));
    const std::string output{ testing::internal::GetCapturedStdout() };
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(previous_log_level);

    ExpectGaussianModelsNear(GetEstimateModel(*atom), previous_model, 1.0e-12);
    EXPECT_NE(output.find("joint-offset = system-build-failed"), std::string::npos);
    EXPECT_NE(
        error_output.find("terminal suspicious rollback fallback clusters/atoms = 1/1"),
        std::string::npos);
    EXPECT_EQ(error_output.find("Reached maximum iteration size"), std::string::npos);
    EXPECT_EQ(FindAcceptedCandidateCluster(output, true), std::string::npos);
    EXPECT_EQ(output.find("Converged after"), std::string::npos);
    EXPECT_EQ(
        error_output.find("Second-stage local fitting fallback summary:"),
        std::string::npos);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, RunSecondStageLocalFittingReportsEmptyJointOffsetSystem)
{
    auto model{ BuildEmptyJointOffsetDefenseModel() };
    auto * atom{ model->GetSelectedAtoms().front() };
    const auto previous_model{ GetEstimateModel(*atom) };

    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions(false));
    const std::string output{ testing::internal::GetCapturedStdout() };
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(previous_log_level);

    EXPECT_NE(output.find("joint-offset = empty-system"), std::string::npos);
    EXPECT_NE(
        output.find("terminal-joint-offset-failure atoms = 1"),
        std::string::npos);
    EXPECT_NE(output.find("Iter. 5/200"), std::string::npos);
    EXPECT_EQ(output.find("Iter. 6/200"), std::string::npos);
    EXPECT_EQ(
        error_output.find("terminal suspicious rollback fallback"),
        std::string::npos);
    EXPECT_NE(
        error_output.find(
            "terminal joint-offset failure fallback clusters/atoms = 1/1"),
        std::string::npos);
    EXPECT_NE(
        error_output.find("statuses = empty-system:1"),
        std::string::npos);
    EXPECT_NE(
        error_output.find("after 5 accepted iterations"),
        std::string::npos);
    EXPECT_EQ(error_output.find("Reached maximum iteration size"), std::string::npos);
    EXPECT_EQ(
        error_output.find("normalized percentile amplitude change"),
        std::string::npos);
    ExpectGaussianModelsNear(GetEstimateModel(*atom), previous_model, 1.0e-12);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, ClusterLocalHealthAllowsRemoteAndersonDuringFallback)
{
    auto model{ BuildSeparatedSystemBuildFailureDefenseModel() };
    const auto initial_remote_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model, 2, 4)
    };

    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    rt::RunSecondStageLocalFitting(
        *model,
        MakeSecondStageOptions(false));
    const std::string output{ testing::internal::GetCapturedStdout() };
    static_cast<void>(testing::internal::GetCapturedStderr());
    Logger::SetLogLevel(previous_log_level);

    EXPECT_NE(
        output.find("joint-offset = system-build-failed"),
        std::string::npos);
    EXPECT_NE(
        output.find("health-unhealthy clusters/atoms = 1/2"),
        std::string::npos);
    EXPECT_NE(
        output.find("joint-offset statuses clusters/atoms = "),
        std::string::npos);
    const auto first_anderson_position{
        FindAcceptedCandidateCluster(output, true)
    };
    const auto terminal_position{
        output.find("terminal-suspicious atoms = 2")
    };
    ASSERT_NE(first_anderson_position, std::string::npos);
    ASSERT_NE(terminal_position, std::string::npos);
    EXPECT_LT(first_anderson_position, terminal_position);
    EXPECT_LT(
        CalculateSelectedAtomResponseMeanSquaredError(*model, 2, 4),
        initial_remote_error);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, LocalRefitFallbackOnlyInvalidatesItsClusterHealth)
{
    auto model{ BuildSeparatedLocalRefitFallbackDefenseModel() };
    const auto initial_remote_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model, 2, 4)
    };

    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions(false));
    const std::string output{ testing::internal::GetCapturedStdout() };
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(previous_log_level);

    EXPECT_NE(
        output.find("health-unhealthy clusters/atoms = 1/2"),
        std::string::npos);
    EXPECT_NE(
        output.find("local-refit-fallback clusters/atoms = 1/1"),
        std::string::npos);
    EXPECT_NE(
        output.find("objective acc./rej. clusters = 2/0, atoms = 4/0"),
        std::string::npos);
    EXPECT_NE(
        output.find("accepted candidate clusters aa/fixed-point = 1/1"),
        std::string::npos);
    EXPECT_NE(FindAcceptedCandidateCluster(output, true), std::string::npos);
    EXPECT_EQ(
        error_output.find("terminal joint-offset failure fallback"),
        std::string::npos);
    EXPECT_NE(
        error_output.find("Reached maximum iteration size"),
        std::string::npos);
    EXPECT_NE(
        error_output.find("applying best validated audit state"),
        std::string::npos);
    EXPECT_NE(
        error_output.find("audit best source = accepted iteration "),
        std::string::npos);
    EXPECT_NE(
        error_output.find("fixed audit objective residual/width/offset/total ="),
        std::string::npos);
    EXPECT_NE(
        error_output.find("offsets finite ="),
        std::string::npos);
    EXPECT_LT(
        CalculateSelectedAtomResponseMeanSquaredError(*model, 2, 4),
        initial_remote_error);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, MaximumIterationAuditPreservesTerminalFallbackAtoms)
{
    auto model{ BuildTerminalWithPersistentLocalRefitFallbackDefenseModel() };
    const auto & selected_atoms{ model->GetSelectedAtoms() };
    const std::array<rg::GaussianModel3D, 2> previous_terminal_model_list{
        GetEstimateModel(*selected_atoms.at(0)),
        GetEstimateModel(*selected_atoms.at(1))
    };

    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions(false));
    const std::string output{ testing::internal::GetCapturedStdout() };
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(previous_log_level);

    for (std::size_t i = 0; i < previous_terminal_model_list.size(); i++)
    {
        ExpectGaussianModelsNear(
            GetEstimateModel(*selected_atoms.at(i)),
            previous_terminal_model_list.at(i),
            1.0e-12);
    }
    EXPECT_NE(output.find("terminal-suspicious atoms = 2"), std::string::npos);
    EXPECT_NE(
        error_output.find("Reached maximum iteration size"),
        std::string::npos);
    EXPECT_NE(
        error_output.find("applying best validated audit state"),
        std::string::npos);
    EXPECT_NE(
        error_output.find("audit best source ="),
        std::string::npos);
    EXPECT_NE(
        error_output.find("fixed audit objective residual/width/offset/total ="),
        std::string::npos);
    EXPECT_NE(
        error_output.find("offsets finite ="),
        std::string::npos);
    EXPECT_NE(
        error_output.find(
            "terminal suspicious rollback fallback clusters/atoms = 1/2"),
        std::string::npos);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, MaximumIterationKeepsLegacyFallbackWhenAuditIsUnavailable)
{
    auto model{ BuildUnavailableAuditWithPersistentLocalRefitFallbackDefenseModel() };

    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions(false));
    static_cast<void>(testing::internal::GetCapturedStdout());
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(previous_log_level);

    EXPECT_NE(
        error_output.find("Reached maximum iteration size"),
        std::string::npos);
    EXPECT_NE(
        error_output.find(
            "applying current accepted candidate because no finite fixed audit state is available"),
        std::string::npos);
    EXPECT_EQ(
        error_output.find("applying best validated audit state"),
        std::string::npos);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, PersistentEmptySystemDoesNotBlockRemoteCluster)
{
    auto model{ BuildSeparatedEmptyJointOffsetDefenseModel() };
    const auto & selected_atoms{ model->GetSelectedAtoms() };
    const auto previous_empty_model{ GetEstimateModel(*selected_atoms.front()) };
    const auto initial_remote_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model, 1, 3)
    };

    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions(false));
    const std::string output{ testing::internal::GetCapturedStdout() };
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(previous_log_level);

    ExpectGaussianModelsNear(
        GetEstimateModel(*selected_atoms.front()),
        previous_empty_model,
        1.0e-12);
    EXPECT_NE(FindAcceptedCandidateCluster(output, true), std::string::npos);
    EXPECT_NE(
        output.find("terminal-joint-offset-failure atoms = 1"),
        std::string::npos);
    EXPECT_NE(
        error_output.find(
            "terminal joint-offset failure fallback clusters/atoms = 1/1"),
        std::string::npos);
    EXPECT_NE(
        error_output.find("statuses = empty-system:1"),
        std::string::npos);
    EXPECT_EQ(
        error_output.find("terminal suspicious rollback fallback"),
        std::string::npos);
    EXPECT_EQ(error_output.find("Reached maximum iteration size"), std::string::npos);
    EXPECT_LT(
        CalculateSelectedAtomResponseMeanSquaredError(*model, 1, 3),
        initial_remote_error);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, RunSecondStageLocalFittingRollsBackFiniteNonphysicalProfile)
{
    auto model{ BuildFiniteNonphysicalProfileDefenseModel() };
    auto * atom{ model->GetSelectedAtoms().front() };
    const auto previous_model{ GetEstimateModel(*atom) };

    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

    ExpectGaussianModelsNear(GetEstimateModel(*atom), previous_model, 1.0e-12);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, RunSecondStageLocalFittingAppliesCollinearRidgeGuard)
{
    auto model{ BuildNearCollinearDefenseModel() };
    const auto initial_error{ CalculateSelectedAtomResponseMeanSquaredError(*model) };

    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

    const auto fitted_error{ CalculateSelectedAtomResponseMeanSquaredError(*model) };
    const auto tolerance{ 1.0e-3 * std::max(initial_error, 1.0) };
    EXPECT_LE(fitted_error, initial_error + tolerance);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, RunSecondStageLocalFittingAcceptsRemoteClusterWhenLocalClusterRollsBack)
{
    auto model{ BuildSeparatedRollbackDefenseModel() };
    const auto & selected_atoms{ model->GetSelectedAtoms() };
    const std::array<rg::GaussianModel3D, 2> previous_left_model_list{
        GetEstimateModel(*selected_atoms.at(0)),
        GetEstimateModel(*selected_atoms.at(1))
    };
    const auto initial_right_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model, 2, 4)
    };

    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions(false));
    const std::string output{ testing::internal::GetCapturedStdout() };
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(previous_log_level);

    for (std::size_t i = 0; i < previous_left_model_list.size(); i++)
    {
        ExpectGaussianModelsNear(
            GetEstimateModel(*selected_atoms.at(i)),
            previous_left_model_list.at(i),
            1.0e-12);
    }
    EXPECT_LT(
        CalculateSelectedAtomResponseMeanSquaredError(*model, 2, 4),
        initial_right_error);
    EXPECT_NE(FindAcceptedCandidateCluster(output, true), std::string::npos);
    EXPECT_NE(output.find("terminal-suspicious atoms = 2"), std::string::npos);
    EXPECT_NE(
        error_output.find("terminal suspicious rollback fallback clusters/atoms = 1/2"),
        std::string::npos);
    EXPECT_EQ(error_output.find("Reached maximum iteration size"), std::string::npos);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, RunSecondStageLocalFittingHandlesMultipleSuspiciousSeedsWithRemoteCluster)
{
    auto model{ BuildSeparatedRollbackDefenseModel() };
    MakeAtomSamplesSuspicious(*model, 1);
    const auto & selected_atoms{ model->GetSelectedAtoms() };
    const std::array<rg::GaussianModel3D, 2> previous_left_model_list{
        GetEstimateModel(*selected_atoms.at(0)),
        GetEstimateModel(*selected_atoms.at(1))
    };
    const auto initial_right_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model, 2, 4)
    };

    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions(false));
    const std::string output{ testing::internal::GetCapturedStdout() };
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(previous_log_level);

    for (std::size_t i = 0; i < previous_left_model_list.size(); i++)
    {
        ExpectGaussianModelsNear(
            GetEstimateModel(*selected_atoms.at(i)),
            previous_left_model_list.at(i),
            1.0e-12);
    }
    EXPECT_LT(
        CalculateSelectedAtomResponseMeanSquaredError(*model, 2, 4),
        initial_right_error);
    EXPECT_NE(output.find("terminal-suspicious atoms = 2"), std::string::npos);
    EXPECT_NE(
        error_output.find("terminal suspicious rollback fallback clusters/atoms = 1/2"),
        std::string::npos);
    EXPECT_EQ(error_output.find("Reached maximum iteration size"), std::string::npos);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, RunSecondStageLocalFittingUsesFixedPointWhenNoAndersonCandidateExists)
{
    auto model{ BuildNearCollinearDefenseModel() };

    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions(false));
    const std::string output{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(previous_log_level);

    EXPECT_NE(
        FindAcceptedCandidateCluster(output, false),
        std::string::npos);
    EXPECT_NE(
        output.find("offset dQ_C p99 raw/accept ="),
        std::string::npos);
    EXPECT_NE(
        output.find("objective acc./rej. clusters = 1/0, atoms = 2/0"),
        std::string::npos);
    EXPECT_EQ(output.find("objective ="), std::string::npos);
    EXPECT_EQ(output.find("d_amplitude ="), std::string::npos);
    EXPECT_EQ(output.find("active/frozen/thawed atoms"), std::string::npos);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, RunSecondStageLocalFittingIsIntensityScaleInvariant)
{
    constexpr double scale{ 100.0 };
    auto base_model{ BuildNearCollinearDefenseModel() };
    auto scaled_model{ BuildNearCollinearDefenseModel(scale) };

    rt::RunSecondStageLocalFitting(*base_model, MakeSecondStageOptions());
    rt::RunSecondStageLocalFitting(*scaled_model, MakeSecondStageOptions());

    const auto & base_atoms{ base_model->GetSelectedAtoms() };
    const auto & scaled_atoms{ scaled_model->GetSelectedAtoms() };
    ASSERT_EQ(base_atoms.size(), scaled_atoms.size());
    for (std::size_t i = 0; i < base_atoms.size(); i++)
    {
        const auto base{ GetEstimateModel(*base_atoms.at(i)) };
        const auto scaled{ GetEstimateModel(*scaled_atoms.at(i)) };
        EXPECT_NEAR(
            base.GetAmplitude() * scale,
            scaled.GetAmplitude(),
            std::max(1.0e-8, std::abs(scaled.GetAmplitude()) * 1.0e-5));
        EXPECT_NEAR(base.GetWidth(), scaled.GetWidth(), 1.0e-6);
        EXPECT_NEAR(
            base.GetOffset() * scale,
            scaled.GetOffset(),
            std::max(1.0e-8, std::abs(scaled.GetOffset()) * 5.0e-5));
    }
}

TEST(EstimatorSecondStageDefenseTest, RunSecondStageLocalFittingKeepsRidgeRetryOnProgressLineWhenPresent)
{
    auto model{ BuildNearCollinearDefenseModel() };

    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions(false));
    const std::string output{ testing::internal::GetCapturedStdout() };
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(previous_log_level);

    const auto retry_position{
        output.find("Objective backtracking rejected all attempts; retrying after")
    };
    if (retry_position != std::string::npos)
    {
        EXPECT_NE(
            output.find("\rObjective backtracking rejected all attempts; retrying after"),
            std::string::npos);
        EXPECT_EQ(
            output.find("\nObjective backtracking rejected all attempts; retrying after"),
            std::string::npos);
        EXPECT_TRUE(
            output.find("increased cluster-local objective ridge", retry_position) != std::string::npos ||
            output.find("increased global ridge ratio", retry_position) != std::string::npos);
        EXPECT_NE(
            output.find("offset dQ_C p99 raw =", retry_position),
            std::string::npos);
        EXPECT_NE(
            output.find(
                "objective acc./rej. clusters = 0/1, atoms = 0/2",
                retry_position),
            std::string::npos);
    }
    EXPECT_EQ(
        error_output.find("best fixed-point candidate"),
        std::string::npos);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, RejectedClusterObjectiveDiagnosticsAreDebugOnly)
{
    auto info_model{ BuildSeparatedLocalRefitFallbackDefenseModel() };
    auto debug_model{ BuildSeparatedLocalRefitFallbackDefenseModel() };

    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    rt::RunSecondStageLocalFitting(*info_model, MakeSecondStageOptions(false));
    const std::string info_output{ testing::internal::GetCapturedStdout() };
    static_cast<void>(testing::internal::GetCapturedStderr());

    Logger::SetLogLevel(LogLevel::Debug);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    rt::RunSecondStageLocalFitting(*debug_model, MakeSecondStageOptions(false));
    const std::string debug_output{ testing::internal::GetCapturedStdout() };
    static_cast<void>(testing::internal::GetCapturedStderr());
    Logger::SetLogLevel(previous_log_level);

    EXPECT_EQ(
        info_output.find("Rejected local fitting cluster objective diagnostics:"),
        std::string::npos);
    EXPECT_NE(
        debug_output.find("Rejected local fitting cluster objective diagnostics:"),
        std::string::npos);
    EXPECT_NE(
        debug_output.find(
            "Rejected local fitting cluster objective diagnostics: atoms = 2"),
        std::string::npos);
    EXPECT_NE(
        debug_output.find("breakdown order = residual/width/offset/movement/total"),
        std::string::npos);
    EXPECT_NE(
        debug_output.find("kind = fixed-point, damping = 1.00e+00 (raw)"),
        std::string::npos);
    for (const auto damping_text : {
        "5.00e-01",
        "2.50e-01",
        "1.25e-01",
        "6.25e-02" })
    {
        EXPECT_NE(
            debug_output.find(
                std::string{ "kind = fixed-point, damping = " } + damping_text),
            std::string::npos);
    }
    EXPECT_NE(debug_output.find("candidate = "), std::string::npos);
    EXPECT_NE(debug_output.find("previous = "), std::string::npos);
    EXPECT_NE(debug_output.find("best = "), std::string::npos);
    EXPECT_NE(debug_output.find("rejected-by = "), std::string::npos);
    const auto & info_atoms{ info_model->GetSelectedAtoms() };
    const auto & debug_atoms{ debug_model->GetSelectedAtoms() };
    ASSERT_EQ(info_atoms.size(), debug_atoms.size());
    for (std::size_t i = 0; i < info_atoms.size(); i++)
    {
        ExpectGaussianModelsNear(
            GetEstimateModel(*info_atoms.at(i)),
            GetEstimateModel(*debug_atoms.at(i)),
            1.0e-12);
    }
    ExpectSelectedAtomEstimatesAreFinite(*info_model);
    ExpectSelectedAtomEstimatesAreFinite(*debug_model);
}

TEST(EstimatorSecondStageDefenseTest, InvalidSeedIsRepairedBeforeConnectedClusterFitting)
{
    auto model{ BuildInvalidCandidateDefenseModel() };

    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Debug);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions(false));
    const std::string output{ testing::internal::GetCapturedStdout() };
    static_cast<void>(testing::internal::GetCapturedStderr());
    Logger::SetLogLevel(previous_log_level);

    EXPECT_NE(
        output.find("Repaired invalid second-stage seed atoms = 1"),
        std::string::npos);
    EXPECT_NE(output.find("source = global-median"), std::string::npos);
    EXPECT_NE(
        output.find("original A/B/C = 0.00e+00/0.00e+00/3.70e-01"),
        std::string::npos);
    EXPECT_NE(
        output.find("repaired A/B/C = 5.75e+00/5.50e-01/3.70e-01"),
        std::string::npos);
    EXPECT_EQ(output.find("reason = non-positive-width"), std::string::npos);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, MissingValidSeedSkipsSecondStageWithoutChangingResults)
{
    auto model{ BuildAllInvalidSeedDefenseModel() };
    std::vector<rg::GaussianModel3D> previous_model_list;
    for (const auto * atom : model->GetSelectedAtoms())
    {
        previous_model_list.emplace_back(GetEstimateModel(*atom));
    }

    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions(false));
    const std::string output{ testing::internal::GetCapturedStdout() };
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(previous_log_level);

    EXPECT_NE(
        (output + error_output).find(
            "Skip 2nd-stage local atom fitting because no valid Gaussian seed"),
        std::string::npos);
    for (std::size_t i = 0; i < model->GetSelectedAtoms().size(); i++)
    {
        ExpectGaussianModelsNear(
            GetEstimateModel(*model->GetSelectedAtoms().at(i)),
            previous_model_list.at(i),
            0.0);
    }
}

TEST(EstimatorSecondStageDefenseTest, RunSecondStageLocalFittingReportsMaximumGlobalRidgeStopWhenReached)
{
    auto model{ BuildNearCollinearDefenseModel() };

    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions(false));
    const std::string output{ testing::internal::GetCapturedStdout() };
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(previous_log_level);

    const auto maximum_ridge_position{
        output.find("next attempt uses increased global ridge ratio = 1.00000")
    };
    const auto stop_warning_position{
        error_output.find("Stopped local fitting because objective backtracking rejected all")
    };
    if (maximum_ridge_position != std::string::npos && stop_warning_position != std::string::npos)
    {
        EXPECT_NE(
            error_output.find("maximum joint-offset ridge ratio", stop_warning_position),
            std::string::npos);
        EXPECT_NE(
            error_output.find("applying best validated audit state", stop_warning_position),
            std::string::npos);
    }
    ExpectSelectedAtomEstimatesAreFinite(*model);
}
