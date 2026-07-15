#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "core/detail/GaussianEstimatorStages.hpp"
#include "core/detail/LocalFittingAudit.hpp"
#include "core/detail/LocalFittingCouplingGraph.hpp"
#include "core/detail/LocalFittingHealth.hpp"
#include "core/detail/LocalFittingJointOffsetConditioning.hpp"
#include "core/detail/LocalFittingJointPolish.hpp"
#include "core/detail/LocalFittingSeedRepair.hpp"
#include "core/detail/LocalFittingTrustRegion.hpp"
#include "core/detail/LocalFittingTransformedChange.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

namespace {
namespace alg = rhbm_gem::algorithm;
namespace audit_detail = rhbm_gem::core::detail;
namespace change_detail = rhbm_gem::core::detail;
namespace conditioning_detail = rhbm_gem::core::detail;
namespace coupling_detail = rhbm_gem::core::detail;
namespace health_detail = rhbm_gem::core::detail;
namespace polish_detail = rhbm_gem::core::detail;
namespace seed_detail = rhbm_gem::core::detail;
namespace trust_detail = rhbm_gem::core::detail;
namespace rt = rhbm_gem::core;
namespace rg = rhbm_gem;

rt::FitOptions MakeSecondStageOptions()
{
    rt::FitOptions options;
    options.distance_min = 0.0;
    options.distance_max = 1.0;
    options.thread_size = 1;
    options.quiet_mode = true;
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

bool HasCouplingNeighbor(
    const coupling_detail::LocalFittingCouplingTopology & topology,
    std::size_t atom_index,
    std::size_t neighbor_index)
{
    const auto & neighbor_index_list{ topology.adjacency_list.at(atom_index) };
    return std::find(
        neighbor_index_list.begin(),
        neighbor_index_list.end(),
        neighbor_index) != neighbor_index_list.end();
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

std::unique_ptr<rg::ModelObject> BuildJointPolishDefenseModel()
{
    return BuildDefenseModel(
        {
            std::array<float, 3>{ 0.0F, 0.0F, 0.0F },
            std::array<float, 3>{ 0.8F, 0.0F, 0.0F }
        },
        { Spot::C, Spot::O },
        { Element::CARBON, Element::OXYGEN },
        {
            rg::GaussianModel3D{ 6.0, 0.45, 0.20 },
            rg::GaussianModel3D{ 4.5, 0.70, -0.12 }
        },
        rg::GaussianModel3D{ 5.25, 0.60, 0.02 });
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

TEST(EstimatorSecondStageDefenseTest, AuditObjectiveProgressGuardChecksPreviousAndBest)
{
    EXPECT_TRUE(audit_detail::IsLocalFittingAuditObjectiveAcceptableForProgress(
        1.0005, 1.0, std::optional<double>{ 1.0 }, 1.0e-3));
    EXPECT_FALSE(audit_detail::IsLocalFittingAuditObjectiveAcceptableForProgress(
        1.002, 1.0, std::optional<double>{ 1.0 }, 1.0e-3));
    EXPECT_FALSE(audit_detail::IsLocalFittingAuditObjectiveAcceptableForProgress(
        1.0, 1.0, std::optional<double>{ 0.99 }, 1.0e-3));
    EXPECT_FALSE(audit_detail::IsLocalFittingAuditObjectiveAcceptableForProgress(
        std::numeric_limits<double>::infinity(),
        1.0,
        std::nullopt,
        1.0e-3));
    EXPECT_FALSE(audit_detail::IsLocalFittingAuditObjectiveAcceptableForProgress(
        std::nullopt, 1.0, std::nullopt, 1.0e-3));
    EXPECT_FALSE(audit_detail::IsLocalFittingAuditObjectiveAcceptableForProgress(
        1.0, std::nullopt, std::nullopt, 1.0e-3));
}

TEST(EstimatorSecondStageDefenseTest, AuditPatienceStopsAfterThreeStableAcceptedIterations)
{
    std::size_t patience_count{ 0 };
    patience_count = audit_detail::AdvanceLocalFittingAuditPatience(
        patience_count, false, false);
    EXPECT_EQ(patience_count, 1U);
    patience_count = audit_detail::AdvanceLocalFittingAuditPatience(
        patience_count, false, false);
    EXPECT_EQ(patience_count, 2U);
    patience_count = audit_detail::AdvanceLocalFittingAuditPatience(
        patience_count, false, false);
    EXPECT_EQ(patience_count, 3U);

    EXPECT_EQ(
        audit_detail::AdvanceLocalFittingAuditPatience(
            patience_count, true, false),
        0U);
    EXPECT_EQ(
        audit_detail::AdvanceLocalFittingAuditPatience(
            patience_count, false, true),
        0U);
}

TEST(EstimatorSecondStageDefenseTest, ScientificObjectiveUsesAtomMeanIndependentOfCardinality)
{
    const auto single_atom{
        audit_detail::BuildLocalFittingMeanObjectiveBreakdown(
            0.4,
            25.0,
            9.0,
            1,
            0.01,
            0.02)
    };
    const auto repeated_atoms{
        audit_detail::BuildLocalFittingMeanObjectiveBreakdown(
            0.4,
            2500.0,
            900.0,
            100,
            0.01,
            0.02)
    };
    const auto single_previous{
        audit_detail::BuildLocalFittingMeanObjectiveBreakdown(
            1.4,
            0.0,
            0.0,
            1,
            0.01,
            0.02)
    };
    const auto repeated_previous{
        audit_detail::BuildLocalFittingMeanObjectiveBreakdown(
            1.4,
            0.0,
            0.0,
            100,
            0.01,
            0.02)
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
        single_atom->total_objective,
        repeated_atoms->total_objective);
    EXPECT_DOUBLE_EQ(
        single_atom->total_objective,
        single_atom->residual_objective +
            single_atom->width_prior_penalty +
            single_atom->offset_plausibility_penalty);
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
            0,
            0.01,
            0.02).has_value());
    EXPECT_FALSE(
        audit_detail::BuildLocalFittingMeanObjectiveBreakdown(
            std::numeric_limits<double>::infinity(),
            25.0,
            9.0,
            1,
            0.01,
            0.02).has_value());
}

TEST(EstimatorSecondStageDefenseTest, TrustRegionDampingCapsLargeTransformedStep)
{
    const std::array<double, 3> scale{ 0.50, 0.35, 1.0 };
    std::vector<Eigen::VectorXd> previous{
        Eigen::VectorXd::Zero(3)
    };
    auto candidate{ previous };
    candidate.at(0) << 1.0, 0.35, 0.5;

    const auto capped{
        trust_detail::LimitLocalFittingTrustRegionDamping(
            previous, candidate, scale, 1.0, 1.0)
    };
    EXPECT_DOUBLE_EQ(capped.effective_damping, 0.5);
    EXPECT_DOUBLE_EQ(capped.step_norm, 1.0);

    const auto inside{
        trust_detail::LimitLocalFittingTrustRegionDamping(
            previous, candidate, scale, 0.25, 1.0)
    };
    EXPECT_DOUBLE_EQ(inside.effective_damping, 0.25);
    EXPECT_DOUBLE_EQ(inside.step_norm, 0.5);
}

TEST(EstimatorSecondStageDefenseTest, TrustRegionDampingIsIntensityScaleInvariant)
{
    const std::array<double, 3> scale{ 0.50, 0.35, 1.0 };
    const auto encode = [](const rg::GaussianModel3D & model)
    {
        const auto estimation{
            change_detail::EncodeLocalFittingTransformedCoordinates(model)
        };
        EXPECT_TRUE(estimation.has_value());
        return estimation.value_or(Eigen::VectorXd{});
    };
    const std::vector<Eigen::VectorXd> base_previous{
        encode(rg::GaussianModel3D{ 2.0, 0.8, 0.2 })
    };
    const std::vector<Eigen::VectorXd> base_candidate{
        encode(rg::GaussianModel3D{ 4.0, 1.0, 0.4 })
    };
    constexpr double intensity_scale{ 1.0e5 };
    const std::vector<Eigen::VectorXd> scaled_previous{
        encode(rg::GaussianModel3D{
            2.0 * intensity_scale,
            0.8,
            0.2 * intensity_scale })
    };
    const std::vector<Eigen::VectorXd> scaled_candidate{
        encode(rg::GaussianModel3D{
            4.0 * intensity_scale,
            1.0,
            0.4 * intensity_scale })
    };

    const auto base{
        trust_detail::LimitLocalFittingTrustRegionDamping(
            base_previous, base_candidate, scale, 1.0, 0.5)
    };
    const auto scaled{
        trust_detail::LimitLocalFittingTrustRegionDamping(
            scaled_previous, scaled_candidate, scale, 1.0, 0.5)
    };
    EXPECT_NEAR(base.effective_damping, scaled.effective_damping, 1.0e-12);
    EXPECT_NEAR(base.step_norm, scaled.step_norm, 1.0e-12);
}

TEST(EstimatorSecondStageDefenseTest, TrustRegionPolishHonorsOuterStepBoundary)
{
    const std::array<double, 3> scale{ 0.50, 0.35, 1.0 };
    std::vector<Eigen::VectorXd> outer_previous{
        Eigen::VectorXd::Zero(3)
    };
    auto boundary_state{ outer_previous };
    boundary_state.at(0)(0) = 0.5;
    auto outward_target{ boundary_state };
    outward_target.at(0)(0) = 1.0;

    const auto outward{
        trust_detail::LimitLocalFittingTrustRegionSubstepDamping(
            outer_previous,
            boundary_state,
            outward_target,
            scale,
            1.0,
            1.0)
    };
    EXPECT_DOUBLE_EQ(outward.effective_damping, 0.0);
    EXPECT_DOUBLE_EQ(outward.step_norm, 1.0);

    const auto inward{
        trust_detail::LimitLocalFittingTrustRegionSubstepDamping(
            outer_previous,
            boundary_state,
            outer_previous,
            scale,
            1.0,
            1.0)
    };
    EXPECT_DOUBLE_EQ(inward.effective_damping, 1.0);
    EXPECT_DOUBLE_EQ(inward.step_norm, 0.0);
}

TEST(EstimatorSecondStageDefenseTest, TrustRegionStateReconcilesShrinksGrowsAndSaturates)
{
    trust_detail::LocalFittingTrustRegionStateSet state;
    const trust_detail::LocalFittingTrustRegionClusterKey key{ 0 };
    state.Reconcile({ key });
    EXPECT_DOUBLE_EQ(state.GetRadius(key), 1.0);

    for (const auto expected : { 0.5, 0.25, 0.125, 0.0625 })
    {
        const auto update{ state.Shrink({ key }) };
        EXPECT_EQ(
            update.changed_key_list,
            std::vector<trust_detail::LocalFittingTrustRegionClusterKey>{ key });
        EXPECT_TRUE(update.saturated_key_list.empty());
        EXPECT_DOUBLE_EQ(state.GetRadius(key), expected);
    }
    const auto saturated{ state.Shrink({ key }) };
    EXPECT_TRUE(saturated.changed_key_list.empty());
    EXPECT_EQ(
        saturated.saturated_key_list,
        std::vector<trust_detail::LocalFittingTrustRegionClusterKey>{ key });

    state.Grow({ key });
    EXPECT_DOUBLE_EQ(state.GetRadius(key), 0.125);

    const trust_detail::LocalFittingTrustRegionClusterKey replacement_key{ 1 };
    state.Reconcile({ replacement_key });
    EXPECT_THROW(state.GetRadius(key), std::invalid_argument);
    EXPECT_DOUBLE_EQ(state.GetRadius(replacement_key), 1.0);
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

TEST(EstimatorSecondStageDefenseTest, LocalRefitHealthTracksStationarity)
{
    EXPECT_TRUE(health_detail::IsLocalGaussianRefitStatusStationarityEligible(
        rg::RHBMEstimationStatus::SUCCESS));
    EXPECT_FALSE(health_detail::IsLocalGaussianRefitStatusStationarityEligible(
        rg::RHBMEstimationStatus::MAX_ITERATIONS_REACHED));
    for (const auto status : {
        rg::RHBMEstimationStatus::NUMERICAL_FALLBACK,
        rg::RHBMEstimationStatus::INSUFFICIENT_DATA,
        rg::RHBMEstimationStatus::SINGLE_MEMBER })
    {
        EXPECT_FALSE(health_detail::IsLocalGaussianRefitStatusStationarityEligible(status));
    }
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
    }
}

TEST(EstimatorSecondStageDefenseTest, JointPolishJacobianMatchesFiniteDifference)
{
    constexpr double step{ 1.0e-6 };
    const std::array<rg::GaussianModel3D, 2> model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 0.20 },
        rg::GaussianModel3D{ 5.5, 0.70, -0.15 }
    };
    const std::array<double, 3> distance_list{ 0.0, 5.0e-6, 0.35 };
    for (const auto & model : model_list)
    {
        const auto transformed{
            change_detail::EncodeLocalFittingTransformedCoordinates(model)
        };
        ASSERT_TRUE(transformed.has_value());
        for (const auto distance : distance_list)
        {
            const auto evaluation{
                polish_detail::EvaluateLocalFittingTransformedResponse(
                    model,
                    distance)
            };
            ASSERT_TRUE(evaluation.has_value());
            EXPECT_TRUE(std::isfinite(evaluation->response));
            EXPECT_TRUE(evaluation->jacobian.allFinite());

            for (Eigen::Index parameter_index = 0;
                parameter_index < transformed->size();
                parameter_index++)
            {
                auto lower{ *transformed };
                auto upper{ *transformed };
                lower(parameter_index) -= step;
                upper(parameter_index) += step;
                const auto lower_model{
                    change_detail::DecodeLocalFittingTransformedCoordinates(lower)
                };
                const auto upper_model{
                    change_detail::DecodeLocalFittingTransformedCoordinates(upper)
                };
                ASSERT_TRUE(lower_model.has_value());
                ASSERT_TRUE(upper_model.has_value());
                const auto finite_difference{
                    (upper_model->ResponseAtDistance(distance) -
                        lower_model->ResponseAtDistance(distance)) /
                    (2.0 * step)
                };
                const auto analytic{
                    evaluation->jacobian(parameter_index)
                };
                const auto tolerance{
                    1.0e-6 * std::max(std::abs(finite_difference), 1.0)
                };
                EXPECT_NEAR(analytic, finite_difference, tolerance);
            }
        }
    }
}

TEST(EstimatorSecondStageDefenseTest, CouplingGraphNormalizesFullJacobianEnergy)
{
    Eigen::Vector3d jacobian;
    jacobian << 1.0, 2.0, 3.0;
    coupling_detail::LocalFittingCouplingGraphBuilder builder{ 2 };
    builder.AddSample({ 0, 0 }, { { 0, jacobian }, { 1, jacobian } });
    builder.AddSample({ 0, 1 }, { { 0, 2.0 * jacobian }, { 1, 2.0 * jacobian } });
    const auto topology{ builder.BuildWeighted(0.05) };
    ASSERT_TRUE(topology.has_value());
    EXPECT_TRUE(HasCouplingNeighbor(*topology, 0, 1));
    EXPECT_NEAR(topology->summary.weight_median, 1.0, 1.0e-12);
    EXPECT_NEAR(topology->summary.weight_percentile_95, 1.0, 1.0e-12);
    EXPECT_NEAR(topology->summary.weight_maximum, 1.0, 1.0e-12);

    coupling_detail::LocalFittingCouplingGraphBuilder scaled_builder{ 2 };
    scaled_builder.AddSample(
        { 0, 0 },
        { { 0, 5.0 * jacobian }, { 1, jacobian } });
    scaled_builder.AddSample(
        { 0, 1 },
        { { 0, 10.0 * jacobian }, { 1, 2.0 * jacobian } });
    const auto scaled_topology{ scaled_builder.BuildWeighted(0.05) };
    ASSERT_TRUE(scaled_topology.has_value());
    EXPECT_TRUE(HasCouplingNeighbor(*scaled_topology, 0, 1));
    EXPECT_NEAR(
        scaled_topology->summary.weight_median,
        topology->summary.weight_median,
        1.0e-12);

    coupling_detail::LocalFittingCouplingGraphBuilder tiny_builder{ 2 };
    tiny_builder.AddSample(
        { 0, 0 },
        { { 0, 1.0e-100 * jacobian }, { 1, 1.0e-100 * jacobian } });
    const auto tiny_topology{ tiny_builder.BuildWeighted(0.05) };
    ASSERT_TRUE(tiny_topology.has_value());
    EXPECT_TRUE(HasCouplingNeighbor(*tiny_topology, 0, 1));
    EXPECT_NEAR(
        tiny_topology->summary.weight_median,
        topology->summary.weight_median,
        1.0e-12);
}

TEST(EstimatorSecondStageDefenseTest, CouplingGraphCutsWeakAndCancelledEdges)
{
    const Eigen::Vector3d unit{ 1.0, 0.0, 0.0 };
    coupling_detail::LocalFittingCouplingGraphBuilder weak_builder{ 2 };
    weak_builder.AddSample({ 0, 0 }, { { 0, unit }, { 1, unit } });
    weak_builder.AddSample({ 0, 1 }, { { 0, 10.0 * unit } });
    weak_builder.AddSample({ 1, 0 }, { { 1, 10.0 * unit } });
    const auto weak_topology{ weak_builder.BuildWeighted(0.05) };
    ASSERT_TRUE(weak_topology.has_value());
    EXPECT_FALSE(HasCouplingNeighbor(*weak_topology, 0, 1));
    EXPECT_EQ(weak_topology->summary.candidate_edge_count, 1U);
    EXPECT_EQ(weak_topology->summary.cut_edge_count, 1U);

    coupling_detail::LocalFittingCouplingGraphBuilder cancelled_builder{ 2 };
    cancelled_builder.AddSample({ 0, 0 }, { { 0, unit }, { 1, unit } });
    cancelled_builder.AddSample({ 0, 1 }, { { 0, unit }, { 1, -unit } });
    const auto cancelled_topology{ cancelled_builder.BuildWeighted(0.05) };
    ASSERT_TRUE(cancelled_topology.has_value());
    EXPECT_FALSE(HasCouplingNeighbor(*cancelled_topology, 0, 1));
}

TEST(EstimatorSecondStageDefenseTest, CouplingGraphReportsThresholdSensitivity)
{
    const Eigen::Vector3d unit{ 1.0, 0.0, 0.0 };
    coupling_detail::LocalFittingCouplingGraphBuilder builder{ 7 };
    const std::array<double, 3> edge_weight_list{ 0.06, 0.12, 0.25 };
    for (std::size_t edge_index = 0; edge_index < edge_weight_list.size(); edge_index++)
    {
        const auto left_index{ 2 * edge_index };
        const auto right_index{ left_index + 1 };
        const auto self_scale{
            std::sqrt(1.0 / edge_weight_list.at(edge_index) - 1.0)
        };
        builder.AddSample(
            { edge_index, 0 },
            { { left_index, unit }, { right_index, unit } });
        builder.AddSample(
            { edge_index, 1 },
            { { left_index, self_scale * unit } });
        builder.AddSample(
            { edge_index, 2 },
            { { right_index, self_scale * unit } });
    }

    const std::vector<double> threshold_list{ 0.05, 0.075, 0.10, 0.15, 0.20, 0.30 };
    const auto topology{ builder.BuildWeighted(0.05, threshold_list) };
    ASSERT_TRUE(topology.has_value());
    ASSERT_EQ(topology->summary.threshold_sensitivity_list.size(), threshold_list.size());

    const std::array<std::size_t, 6> retained_edge_count_list{ 3, 2, 2, 1, 1, 0 };
    for (std::size_t i = 0; i < threshold_list.size(); i++)
    {
        const auto & sensitivity{ topology->summary.threshold_sensitivity_list.at(i) };
        const auto retained_edge_count{ retained_edge_count_list.at(i) };
        EXPECT_DOUBLE_EQ(sensitivity.minimum_weight, threshold_list.at(i));
        EXPECT_EQ(sensitivity.retained_edge_count, retained_edge_count);
        EXPECT_EQ(sensitivity.cut_edge_count, 3U - retained_edge_count);
        EXPECT_EQ(sensitivity.component_count, 7U - retained_edge_count);
        EXPECT_EQ(sensitivity.maximum_component_size, retained_edge_count == 0 ? 1U : 2U);
        EXPECT_NEAR(
            sensitivity.maximum_component_ratio,
            static_cast<double>(sensitivity.maximum_component_size) / 7.0,
            1.0e-12);
    }

    const auto & formal_threshold{ topology->summary.threshold_sensitivity_list.front() };
    EXPECT_EQ(formal_threshold.retained_edge_count, topology->summary.retained_edge_count);
    EXPECT_EQ(
        topology->summary.candidate_edge_count - formal_threshold.retained_edge_count,
        topology->summary.cut_edge_count);
    const auto formal_partition{
        coupling_detail::BuildLocalFittingCouplingPartition(
            *topology,
            { 0, 1, 2, 3, 4, 5, 6 })
    };
    EXPECT_EQ(formal_threshold.component_count, formal_partition.sample_id_list_by_key.size());
    EXPECT_EQ(formal_threshold.maximum_component_size, 2U);
}

TEST(EstimatorSecondStageDefenseTest, CouplingPartitionCutsWeakBridgeAndDuplicatesBoundarySample)
{
    coupling_detail::LocalFittingCouplingTopology topology;
    topology.adjacency_list.resize(3);
    topology.adjacency_list.at(0).push_back(1);
    topology.adjacency_list.at(1).push_back(0);
    topology.sample_dependency_list = {
        { { 0, 0 }, { 0, 1 } },
        { { 1, 0 }, { 1, 2 } }
    };

    const auto partition{
        coupling_detail::BuildLocalFittingCouplingPartition(topology, { 0, 1, 2 })
    };
    ASSERT_EQ(partition.sample_id_list_by_key.size(), 2U);
    EXPECT_EQ(partition.sample_id_list_by_key.count({ 0, 1 }), 1U);
    EXPECT_EQ(partition.sample_id_list_by_key.count({ 2 }), 1U);
    EXPECT_EQ(partition.boundary_sample_count, 1U);
    EXPECT_EQ(partition.sample_id_list_by_key.at({ 0, 1 }).size(), 2U);
    EXPECT_EQ(partition.sample_id_list_by_key.at({ 2 }).size(), 1U);

    const auto inactive_partition{
        coupling_detail::BuildLocalFittingCouplingPartition(topology, { 2, 0 })
    };
    EXPECT_EQ(inactive_partition.sample_id_list_by_key.count({ 0 }), 1U);
    EXPECT_EQ(inactive_partition.sample_id_list_by_key.count({ 2 }), 1U);
    EXPECT_EQ(inactive_partition.boundary_sample_count, 0U);
}

TEST(EstimatorSecondStageDefenseTest, CouplingPartitionKeepsStrongChainAndBinaryFallback)
{
    coupling_detail::LocalFittingCouplingTopology strong_topology;
    strong_topology.adjacency_list = {
        { 1 },
        { 0, 2 },
        { 1 }
    };
    const auto strong_partition{
        coupling_detail::BuildLocalFittingCouplingPartition(
            strong_topology,
            { 0, 1, 2 })
    };
    EXPECT_EQ(strong_partition.sample_id_list_by_key.count({ 0, 1, 2 }), 1U);

    coupling_detail::LocalFittingCouplingGraphBuilder builder{ 2 };
    const auto invalid{
        Eigen::Vector3d::Constant(std::numeric_limits<double>::quiet_NaN())
    };
    builder.AddSample(
        { 0, 0 },
        { { 0, Eigen::Vector3d::Ones() }, { 1, invalid } });
    EXPECT_FALSE(builder.BuildWeighted(0.05, { 0.05, 0.10 }).has_value());
    const auto binary_topology{ builder.BuildBinary() };
    EXPECT_FALSE(binary_topology.summary.uses_weighted_graph);
    EXPECT_TRUE(binary_topology.summary.threshold_sensitivity_list.empty());
    const auto binary_partition{
        coupling_detail::BuildLocalFittingCouplingPartition(
            binary_topology,
            { 0, 1 })
    };
    EXPECT_EQ(binary_partition.sample_id_list_by_key.count({ 0, 1 }), 1U);

    coupling_detail::LocalFittingCouplingGraphBuilder overflow_builder{ 2 };
    const auto huge{ Eigen::Vector3d::Constant(1.0e200) };
    overflow_builder.AddSample({ 0, 0 }, { { 0, huge }, { 1, huge } });
    EXPECT_FALSE(overflow_builder.BuildWeighted(0.05).has_value());
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

    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

    for (std::size_t i = 0; i < selected_atoms.size(); i++)
    {
        ExpectGaussianModelsNear(
            GetEstimateModel(*selected_atoms.at(i)),
            previous_model_list.at(i),
            1.0e-12);
    }
}

TEST(EstimatorSecondStageDefenseTest, RunSecondStageLocalFittingFallsBackWhenJointOffsetSamplesAreNonFinite)
{
    auto model{ BuildNonFiniteJointOffsetDefenseModel() };
    auto * atom{ model->GetSelectedAtoms().front() };
    const auto previous_model{ GetEstimateModel(*atom) };

    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

    ExpectGaussianModelsNear(GetEstimateModel(*atom), previous_model, 1.0e-12);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, SystemBuildFailureDoesNotBlockRemoteCluster)
{
    auto model{ BuildSeparatedSystemBuildFailureDefenseModel() };
    const auto initial_remote_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model, 2, 4)
    };

    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

    EXPECT_LT(
        CalculateSelectedAtomResponseMeanSquaredError(*model, 2, 4),
        initial_remote_error);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, LocalRefitFallbackDoesNotBlockRemoteCluster)
{
    auto model{ BuildSeparatedLocalRefitFallbackDefenseModel() };
    const auto initial_remote_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model, 2, 4)
    };

    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

    EXPECT_LT(
        CalculateSelectedAtomResponseMeanSquaredError(*model, 2, 4),
        initial_remote_error);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, TerminalFallbackPreservesAffectedCluster)
{
    auto model{ BuildTerminalWithPersistentLocalRefitFallbackDefenseModel() };
    const auto & selected_atoms{ model->GetSelectedAtoms() };
    const std::array<rg::GaussianModel3D, 2> previous_terminal_model_list{
        GetEstimateModel(*selected_atoms.at(0)),
        GetEstimateModel(*selected_atoms.at(1))
    };

    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

    for (std::size_t i = 0; i < previous_terminal_model_list.size(); i++)
    {
        ExpectGaussianModelsNear(
            GetEstimateModel(*selected_atoms.at(i)),
            previous_terminal_model_list.at(i),
            1.0e-12);
    }
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

    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

    ExpectGaussianModelsNear(
        GetEstimateModel(*selected_atoms.front()),
        previous_empty_model,
        1.0e-12);
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

TEST(EstimatorSecondStageDefenseTest, RunSecondStageLocalFittingJointlyPolishesClusterParameters)
{
    auto model{ BuildJointPolishDefenseModel() };
    const auto & atom_list{ model->GetSelectedAtoms() };
    std::vector<rg::GaussianModel3D> initial_model_list;
    for (const auto * atom : atom_list)
    {
        initial_model_list.emplace_back(GetEstimateModel(*atom));
    }
    const auto initial_error{ CalculateSelectedAtomResponseMeanSquaredError(*model) };

    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

    const auto fitted_error{ CalculateSelectedAtomResponseMeanSquaredError(*model) };
    EXPECT_LT(fitted_error, initial_error);
    bool amplitude_changed{ false };
    bool width_changed{ false };
    bool offset_changed{ false };
    for (std::size_t i = 0; i < atom_list.size(); i++)
    {
        const auto fitted{ GetEstimateModel(*atom_list.at(i)) };
        const auto & initial{ initial_model_list.at(i) };
        amplitude_changed = amplitude_changed ||
            std::abs(fitted.GetAmplitude() - initial.GetAmplitude()) > 1.0e-6;
        width_changed = width_changed ||
            std::abs(fitted.GetWidth() - initial.GetWidth()) > 1.0e-6;
        offset_changed = offset_changed ||
            std::abs(fitted.GetOffset() - initial.GetOffset()) > 1.0e-6;
    }
    EXPECT_TRUE(amplitude_changed);
    EXPECT_TRUE(width_changed);
    EXPECT_TRUE(offset_changed);
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

    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

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

TEST(EstimatorSecondStageDefenseTest, MissingValidSeedSkipsSecondStageWithoutChangingResults)
{
    auto model{ BuildAllInvalidSeedDefenseModel() };
    std::vector<rg::GaussianModel3D> previous_model_list;
    for (const auto * atom : model->GetSelectedAtoms())
    {
        previous_model_list.emplace_back(GetEstimateModel(*atom));
    }

    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

    for (std::size_t i = 0; i < model->GetSelectedAtoms().size(); i++)
    {
        ExpectGaussianModelsNear(
            GetEstimateModel(*model->GetSelectedAtoms().at(i)),
            previous_model_list.at(i),
            0.0);
    }
}

TEST(EstimatorSecondStageDefenseTest, NonQuietSecondStageLogsEveryOuterAttempt)
{
    auto model{ BuildSeparatedRollbackDefenseModel() };
    auto options{ MakeSecondStageOptions() };
    options.quiet_mode = false;
    const auto previous_log_level{ Logger::GetLogLevel() };

    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    rt::RunSecondStageLocalFitting(*model, options);
    const std::string out{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(previous_log_level);

    const auto count_occurrences = [&](std::string_view text)
    {
        std::size_t count{ 0 };
        for (std::size_t position = 0;
            (position = out.find(text, position)) != std::string::npos;
            position += text.size())
        {
            count++;
        }
        return count;
    };
    EXPECT_EQ(count_occurrences("Try/Acc"), 1U);
    EXPECT_NE(out.find("Atom A/T"), std::string::npos);
    EXPECT_NE(out.find("Cluster A/R"), std::string::npos);
    EXPECT_NE(out.find("Suspicious"), std::string::npos);
    EXPECT_NE(out.find("dMax A/R"), std::string::npos);

    const auto header_start{ out.find("Try/Acc") };
    ASSERT_NE(header_start, std::string::npos);
    const auto header_end{ out.find('\n', header_start) };
    ASSERT_NE(header_end, std::string::npos);
    const std::string_view header{
        out.data() + header_start,
        header_end - header_start
    };
    std::vector<std::string_view> progress_row_list;
    for (std::size_t row_start = out.find('\r');
        row_start != std::string::npos;
        row_start = out.find('\r', row_start))
    {
        row_start++;
        const auto row_end{ out.find_first_of("\r\n", row_start) };
        ASSERT_NE(row_end, std::string::npos);
        progress_row_list.emplace_back(
            out.data() + row_start,
            row_end - row_start);
        row_start = row_end;
    }
    const auto separator_position_list = [](std::string_view row)
    {
        std::vector<std::size_t> position_list;
        for (std::size_t position = row.find('|');
            position != std::string::npos;
            position = row.find('|', position + 1))
        {
            position_list.emplace_back(position);
        }
        return position_list;
    };
    const auto header_separator_position_list{
        separator_position_list(header)
    };
    ASSERT_EQ(header_separator_position_list.size(), 4U);
    ASSERT_EQ(progress_row_list.size(), 6U);
    for (const auto row : progress_row_list)
    {
        EXPECT_EQ(row.size(), header.size());
        EXPECT_EQ(
            separator_position_list(row),
            header_separator_position_list);
    }
    EXPECT_NE(
        progress_row_list.front().find("3.58e-02/4.14e-02"),
        std::string::npos);

    const std::string summary_prefix{
        "Second-stage local fitting summary: accepted_iterations="
    };
    const auto summary_position{ out.find(summary_prefix) };
    ASSERT_NE(summary_position, std::string::npos);
    ASSERT_GT(summary_position, 0U);
    EXPECT_EQ(out.at(summary_position - 1), '\n');
    const auto accepted_start{ summary_position + summary_prefix.size() };
    const auto accepted_end{ out.find(',', accepted_start) };
    ASSERT_NE(accepted_end, std::string::npos);
    const auto accepted_iteration_count{
        static_cast<std::size_t>(std::stoull(
            out.substr(accepted_start, accepted_end - accepted_start)))
    };
    EXPECT_EQ(accepted_iteration_count, 5U);
    EXPECT_NE(out.find("-/"), std::string::npos);
    EXPECT_NE(
        out.find(
            "best_iteration=5, stop_reason=all-rejected-minimum-radius"),
        std::string::npos);
}

TEST(EstimatorSecondStageDefenseTest, QuietSecondStageSuppressesIterationTable)
{
    auto model{ BuildSeparatedRollbackDefenseModel() };
    const auto previous_log_level{ Logger::GetLogLevel() };

    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());
    const std::string out{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(previous_log_level);

    EXPECT_EQ(out.find("Try/Acc"), std::string::npos);
    EXPECT_EQ(std::count(out.begin(), out.end(), '\r'), 0);
}
