#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

namespace {
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

std::unique_ptr<rg::ModelObject> BuildNearCollinearDefenseModel()
{
    return BuildDefenseModel(
        {
            std::array<float, 3>{ 0.0F, 0.0F, 0.0F },
            std::array<float, 3>{ 1.0e-4F, 0.0F, 0.0F }
        },
        { Spot::C, Spot::O },
        { Element::CARBON, Element::OXYGEN },
        {
            rg::GaussianModel3D{ 6.0, 0.55, 0.20 },
            rg::GaussianModel3D{ 5.5, 0.55, -0.15 }
        },
        rg::GaussianModel3D{ 5.75, 0.55, 0.0 });
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

TEST(EstimatorSecondStageDefenseTest, RunSecondStageLocalFittingFallsBackWhenJointOffsetSamplesAreNonFinite)
{
    auto model{ BuildNonFiniteJointOffsetDefenseModel() };
    auto * atom{ model->GetSelectedAtoms().front() };
    const auto previous_model{ GetEstimateModel(*atom) };

    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Warning);
    testing::internal::CaptureStderr();
    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions(false));
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(previous_log_level);

    ExpectGaussianModelsNear(GetEstimateModel(*atom), previous_model, 1.0e-12);
    EXPECT_EQ(
        error_output.find("Second-stage local fitting fallback summary:"),
        std::string::npos);
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
    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions(false));
    const std::string output{ testing::internal::GetCapturedStdout() };
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
    EXPECT_TRUE(
        output.find("acceleration = aa") != std::string::npos ||
        output.find("acceleration = damped-aa") != std::string::npos);
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
        output.find("acceleration = damped-fixed-point"),
        std::string::npos);
    EXPECT_EQ(output.find("objective ="), std::string::npos);
    EXPECT_EQ(output.find("d_amplitude ="), std::string::npos);
    EXPECT_EQ(output.find("active/frozen/thawed atoms"), std::string::npos);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, RunSecondStageLocalFittingReportsClusterLocalRidgeRetry)
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

    EXPECT_NE(
        output.find("Objective backtracking rejected all attempts; retrying after"),
        std::string::npos);
    EXPECT_NE(
        output.find("increased cluster-local objective ridge"),
        std::string::npos);
    EXPECT_EQ(
        error_output.find("best fixed-point candidate"),
        std::string::npos);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, RunSecondStageLocalFittingRetriesAfterIncreasingGlobalRidgeToMaximum)
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

    EXPECT_NE(
        output.find("next attempt uses increased global ridge ratio = 1.00000"),
        std::string::npos);
    const auto stop_warning_position{
        error_output.find("Stopped local fitting because objective backtracking rejected all")
    };
    if (stop_warning_position != std::string::npos)
    {
        EXPECT_NE(
            error_output.find("maximum joint-offset ridge ratio", stop_warning_position),
            std::string::npos);
    }
    ExpectSelectedAtomEstimatesAreFinite(*model);
}
