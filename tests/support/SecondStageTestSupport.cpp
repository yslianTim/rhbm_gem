#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "support/SecondStageTestSupport.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>

namespace second_stage_test {

namespace detail = rhbm_gem::core::detail;
namespace rg = rhbm_gem;
namespace rt = rhbm_gem::core;
using rhbm_gem::FittingStage;

namespace {

void MakeAtomSamplesSuspicious(rg::ModelObject & model, std::size_t atom_index)
{
    const auto & selected_atoms{ model.GetSelectedAtoms() };
    auto * target_atom{ selected_atoms.at(atom_index) };
    const auto target_position{ target_atom->GetPosition() };

    auto raw_sampling_entries{
        rg::AtomLocalPotentialView::For(*target_atom)
            .GetRawSamplingEntries(false)
    };
    raw_sampling_entries.resize(256);
    raw_sampling_entries.front().response = 0.0;
    raw_sampling_entries.front().point.distance = 0.0;
    raw_sampling_entries.front().point.position = target_position;
    for (std::size_t i = 1; i < raw_sampling_entries.size(); i++)
    {
        auto & sample{ raw_sampling_entries.at(i) };
        const auto response_scale{
            0.5 + 0.5 * static_cast<double>(i) /
                static_cast<double>(raw_sampling_entries.size())
        };
        sample.response = std::numeric_limits<double>::max() * response_scale;
        sample.point.position = target_position;
        sample.point.position.at(0) += 100.0;
        sample.point.distance = 100.0;
    }

    auto analysis{ model.EditAnalysis() };
    analysis.SetAtomLocalRawSamplingEntries(
        *target_atom, std::move(raw_sampling_entries));
}

} // namespace

rt::FitOptions MakeSecondStageOptions()
{
    rt::FitOptions options;
    options.thread_size = 1;
    options.quiet_mode = true;
    return options;
}

double Distance(
    const std::array<double, 3> & lhs,
    const std::array<double, 3> & rhs)
{
    const auto dx{ lhs.at(0) - rhs.at(0) };
    const auto dy{ lhs.at(1) - rhs.at(1) };
    const auto dz{ lhs.at(2) - rhs.at(2) };
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

std::unique_ptr<rg::AtomObject> MakeAtom(
    int serial_id,
    Spot spot,
    Element element,
    const std::array<double, 3> & position)
{
    auto atom{ std::make_unique<rg::AtomObject>() };
    atom->SetSerialID(serial_id);
    atom->SetChainID("A");
    atom->SetSequenceID(serial_id);
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
    const std::array<std::array<double, 3>, 6> direction_list{
        std::array<double, 3>{ 1.0, 0.0, 0.0 },
        std::array<double, 3>{ -1.0, 0.0, 0.0 },
        std::array<double, 3>{ 0.0, 1.0, 0.0 },
        std::array<double, 3>{ 0.0, -1.0, 0.0 },
        std::array<double, 3>{ 0.0, 0.0, 1.0 },
        std::array<double, 3>{ 0.0, 0.0, -1.0 }
    };
    const std::array<double, 4> radius_list{ 0.15, 0.35, 0.65, 0.95 };

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
            sample_list.emplace_back(LocalPotentialSample{ response, point });
        }
    }
    return sample_list;
}

JointPolishFixture BuildJointPolishFixture(
    const std::vector<rg::GaussianModel3D> & base_model_list,
    const std::vector<rg::GaussianModel3D> & target_model_list)
{
    if (base_model_list.size() != target_model_list.size() ||
        base_model_list.empty())
    {
        throw std::invalid_argument(
            "Joint polish fixture input sizes are inconsistent.");
    }

    constexpr std::array<double, 5> distance_list{
        0.0,
        0.15,
        0.30,
        0.45,
        0.60
    };
    JointPolishFixture fixture;
    fixture.context.atom_list.resize(base_model_list.size());
    fixture.state.reserve(base_model_list.size());
    for (std::size_t atom_index = 0;
        atom_index < base_model_list.size();
        atom_index++)
    {
        auto & atom_context{ fixture.context.atom_list.at(atom_index) };
        atom_context.neighbor_atom_sample_offset_list.assign(
            distance_list.size() + 1,
            0);
        for (const auto distance : distance_list)
        {
            atom_context.raw_sampling_entries.emplace_back(
                LocalPotentialSample{
                    target_model_list.at(atom_index).ResponseAtDistance(distance),
                    SamplingPoint{ distance }
                });
            fixture.sample_ref_list.emplace_back(
                detail::SampleRef{
                    atom_index,
                    atom_context.raw_sampling_entries.size() - 1
                });
        }
        fixture.state.emplace_back(MakeGaussianResult(
            base_model_list.at(atom_index)));

    }
    return fixture;
}

std::unique_ptr<rg::ModelObject> BuildDefenseModel(
    const std::vector<std::array<double, 3>> & position_list,
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
        analysis.SetAtomLocalAlphaR(FittingStage::Second, *atom, 0.0);
        analysis.SetAtomLocalGaussianResult(
            FittingStage::Second,
            *atom,
            MakeGaussianResult(initial_model));
        analysis.SetAtomLocalRawSamplingEntries(
            *atom, BuildSamples(*atom, selected_atoms, truth_model_list));
    }
    return model;
}

std::unique_ptr<rg::ModelObject> BuildNearCollinearDefenseModel(
    double intensity_scale)
{
    return BuildDefenseModel(
        {
            std::array<double, 3>{ 0.0, 0.0, 0.0 },
            std::array<double, 3>{ 1.0e-4, 0.0, 0.0 }
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
    auto model{ BuildDefenseModel(
        {
            std::array<double, 3>{ 0.0, 0.0, 0.0 },
            std::array<double, 3>{ 0.8, 0.0, 0.0 }
        },
        { Spot::C, Spot::O },
        { Element::CARBON, Element::OXYGEN },
        {
            rg::GaussianModel3D{ 6.0, 0.45, 0.20 },
            rg::GaussianModel3D{ 4.5, 0.70, -0.12 }
        },
        rg::GaussianModel3D{ 5.25, 0.60, 0.02 }) };
    auto analysis{ model->EditAnalysis() };
    const auto & atom_list{ model->GetSelectedAtoms() };
    // Keep the raw joint-offset update inside the trust region now that a
    // noise-scale rebound no longer supplies an incidental rollback.
    analysis.SetAtomLocalGaussianResult(
        FittingStage::Second,
        *atom_list.at(0),
        MakeGaussianResult(rg::GaussianModel3D{ 5.8, 0.47, 0.15 }));
    analysis.SetAtomLocalGaussianResult(
        FittingStage::Second,
        *atom_list.at(1),
        MakeGaussianResult(rg::GaussianModel3D{ 4.7, 0.68, -0.08 }));
    return model;
}

std::unique_ptr<rg::ModelObject> BuildIndependentOffsetDefenseModel(
    double intensity_scale,
    bool alternate_keys)
{
    auto model{ BuildDefenseModel(
        {
            std::array<double, 3>{ 0.0, 0.0, 0.0 },
            std::array<double, 3>{ 0.8, 0.0, 0.0 }
        },
        alternate_keys ? std::vector<Spot>{ Spot::N, Spot::O } :
            std::vector<Spot>{ Spot::C, Spot::C },
        { Element::CARBON, Element::CARBON },
        {
            rg::GaussianModel3D{
                6.0 * intensity_scale,
                0.48,
                0.05 * intensity_scale
            },
            rg::GaussianModel3D{
                4.8 * intensity_scale,
                0.62,
                -0.08 * intensity_scale
            }
        },
        rg::GaussianModel3D{
            5.3 * intensity_scale,
            0.57,
            0.02 * intensity_scale
        }) };

    const std::array<double, 2> initial_offset_list{
        0.03 * intensity_scale,
        0.07 * intensity_scale
    };
    auto analysis{ model->EditAnalysis() };
    const auto & atom_list{ model->GetSelectedAtoms() };
    for (std::size_t i = 0; i < atom_list.size(); i++)
    {
        const auto current{
            rg::AtomLocalPotentialView::For(*atom_list.at(i))
                .GetGaussianResult(FittingStage::Second).mdpde.GetModel()
        };
        analysis.SetAtomLocalGaussianResult(
            FittingStage::Second,
            *atom_list.at(i),
            MakeGaussianResult(current.WithOffset(initial_offset_list.at(i))));
    }
    return model;
}

std::unique_ptr<rg::ModelObject> BuildSeparatedRollbackDefenseModel()
{
    auto model{
        BuildDefenseModel(
            {
                std::array<double, 3>{ 0.0, 0.0, 0.0 },
                std::array<double, 3>{ 1.0e-4, 0.0, 0.0 },
                std::array<double, 3>{ 10.0, 0.0, 0.0 },
                std::array<double, 3>{ 10.0001, 0.0, 0.0 }
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
        const auto local_view{ rg::AtomLocalPotentialView::For(*atom) };
        for (const auto & sample : local_view.GetRawSamplingEntries(false))
        {
            double fitted_response{ 0.0 };
            for (const auto * fitted_atom : selected_atoms)
            {
                const auto fitted_view{
                    rg::AtomLocalPotentialView::For(*fitted_atom)
                };
                fitted_response += fitted_view.GetEstimateMDPDE(
                    FittingStage::Second).ResponseAtDistance(
                    Distance(sample.point.position, fitted_atom->GetPosition()));
            }
            const auto residual{ sample.response - fitted_response };
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
    return rg::AtomLocalPotentialView::For(atom).GetEstimateMDPDE(
        FittingStage::Second);
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
            rg::AtomLocalPotentialView::For(*atom).GetEstimateMDPDE(
                FittingStage::Second)
        };
        EXPECT_TRUE(std::isfinite(estimate.GetAmplitude()));
        EXPECT_TRUE(std::isfinite(estimate.GetWidth()));
        EXPECT_TRUE(std::isfinite(estimate.GetOffset()));
    }
}

} // namespace second_stage_test
