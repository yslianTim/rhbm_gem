#include "core/estimator/GaussianEstimator.hpp"
#include "core/estimator/InitialStage.hpp"
#include "core/estimator/joint_component/Numerics.hpp"
#include "core/estimator/gaussian_fit/PreparedLocalGaussianFit.hpp"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/hrl/RHBMTrainer.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

namespace rhbm_gem::core {
namespace {
constexpr double kLocalGaussianFittingRadius{1.0};

constexpr std::size_t kMinimumAlphaRTrainingSampleCount{ 10 };
rhbm_trainer::RHBMTrainingOptions MakeTrainingOptions(const FitOptions & options)
{
    rhbm_trainer::RHBMTrainingOptions training_options;
    training_options.execution_options = RHBMExecutionOptions{
        .thread_size = options.thread_size
    };
    return training_options;
}

std::vector<double> CollectSampleResponses(const LocalPotentialSampleList & sample_entries)
{
    std::vector<double> response_list;
    response_list.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        response_list.emplace_back(sample.response);
    }
    return response_list;
}

GaussianModel3DWithUncertainty WithPreservedUncertaintyOffset(
    const GaussianModel3DWithUncertainty & gaussian,
    double offset)
{
    return {
        gaussian.GetModel().WithOffset(offset),
        gaussian.GetStandardDeviationModel()
    };
}

std::vector<GroupGaussianMemberResult> DecodeMemberGaussianResults(
    const RHBMGroupEstimationResult & result,
    const std::vector<double> & member_offset_list)
{
    const auto member_count{ static_cast<std::size_t>(result.beta_posterior_matrix.cols()) };
    if (member_offset_list.size() != member_count ||
        result.capital_sigma_posterior_list.size() != member_count)
    {
        throw std::invalid_argument("Group Gaussian member result count is inconsistent.");
    }
    eigen_validation::RequireVectorSize(
        result.outlier_flag_array, result.beta_posterior_matrix.cols(),
        "outlier_flag_array", "Group Gaussian member result count is inconsistent.");
    eigen_validation::RequireVectorSize(
        result.statistical_distance_array, result.beta_posterior_matrix.cols(),
        "statistical_distance_array", "Group Gaussian member result count is inconsistent.");

    std::vector<GroupGaussianMemberResult> member_results;
    member_results.reserve(member_count);
    for (Eigen::Index i = 0; i < result.beta_posterior_matrix.cols(); i++)
    {
        const auto member_index{ static_cast<std::size_t>(i) };
        const auto offset{ member_offset_list.at(member_index) };
        const auto gaussian{
            linearization_service::DecodeParameterVector(
                result.beta_posterior_matrix.col(i),
                result.capital_sigma_posterior_list.at(member_index))
        };
        const auto gaussian_with_offset{
            WithPreservedUncertaintyOffset(gaussian, offset)
        };
        member_results.emplace_back(GroupGaussianMemberResult{
            gaussian_with_offset,
            static_cast<bool>(result.outlier_flag_array(i)),
            result.statistical_distance_array(i),
            {},
            true,
            std::nullopt,
            std::nullopt
        });
    }
    return member_results;
}

GroupGaussianResult DecodeGroupGaussianResult(
    double alpha_g,
    const RHBMGroupEstimationResult & result,
    const std::vector<double> & member_offset_list)
{
    const auto group_offset{ array_helper::ComputeMedian(member_offset_list) };
    const auto prior{
        linearization_service::DecodeParameterVector(result.mu_prior, result.capital_lambda)
    };
    const auto mean{
        linearization_service::DecodeParameterVector(result.mu_mean).WithOffset(group_offset)
    };
    const auto mdpde{
        linearization_service::DecodeParameterVector(result.mu_mdpde).WithOffset(group_offset)
    };
    return GroupGaussianResult{
        alpha_g,
        mean,
        mdpde,
        WithPreservedUncertaintyOffset(prior, group_offset),
        DecodeMemberGaussianResults(result, member_offset_list)
    };
}

} // namespace

double EvaluateGroupMemberResponseAtDistance(const GroupGaussianMemberResult & member, double distance)
{
    const auto & gaussian = member.posterior.GetModel();
    if (!member.fixed_charge) return gaussian.ResponseAtDistance(distance);
    const auto charge = joint_component::EvaluateKernel(
        distance * distance, member.fixed_charge->width, 2.5).charge;
    return gaussian.SignalAtDistance(distance) + member.fixed_charge->coefficient * charge;
}

void RunFirstStageLocalFitting(ModelObject & model_object, const FitOptions & options)
{
    model_object.EditAnalysis().EnsureSelectedAtomLocalPotentials();
    RunFirstStageLocalFitting(model_object, options, model_object.GetSelectedAtoms());
}

void RunFirstStageLocalFitting(
    ModelObject & model_object,
    const FitOptions & options,
    const std::vector<AtomObject *> & atom_list)
{
    const auto selected_atom_size{ atom_list.size() };
    auto analysis{ model_object.EditAnalysis() };
    std::vector<LocalGaussianResult> local_results(selected_atom_size);
    size_t atom_count{ 0 };
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info,
            "Run local atom fitting for " + std::to_string(selected_atom_size) + " atoms.");
    }

#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(options.thread_size)
#endif
    for (size_t i = 0; i < selected_atom_size; i++)
    {
        auto & atom{ *atom_list[i] };
        local_results[i] = detail::FitFirstStageAtom(atom, options);

        if (!options.quiet_mode)
        {
#ifdef USE_OPENMP
            #pragma omp critical
#endif
            {
                atom_count++;
                Logger::ProgressPercent(atom_count, selected_atom_size);
            }
        }
    }

    for (size_t i = 0; i < selected_atom_size; i++)
    {
        analysis.SetAtomFirstStageGaussianResult(*atom_list[i], std::move(local_results[i]));
    }
}

double TrainAlphaR(
    const LocalPotentialSampleList & sample_entries,
    const FitOptions & options)
{
    return detail::TrainLocalGaussianAlpha(
        rhbm_helper::BuildMemberDataset(sample_entries,0.0,kLocalGaussianFittingRadius),options.thread_size);
}

double TrainAlphaG(
    const std::vector<std::vector<GaussianModel3D>> & model_group_list,
    const FitOptions & options)
{
    std::vector<std::vector<RHBMParameterVector>> beta_group_list;
    beta_group_list.reserve(model_group_list.size());
    for (const auto & model_group : model_group_list)
    {
        std::vector<RHBMParameterVector> beta_list;
        beta_list.reserve(model_group.size());
        for (const auto & model : model_group)
        {
            beta_list.emplace_back(
                linearization_service::EncodeGaussianToParameterVector(model));
        }
        beta_group_list.emplace_back(std::move(beta_list));
    }

    const auto training_options{ MakeTrainingOptions(options) };
    if (beta_group_list.empty())
    {
        return training_options.alpha_min;
    }

    return rhbm_trainer::CrossValidationAlphaG(beta_group_list, training_options).best_alpha;
}

LocalGaussianResult EstimateLocalGaussian(
    const LocalPotentialSampleList & sample_entries,
    double alpha_r,
    const FitOptions & options,
    const GaussianModel3D & offset_model)
{
    const detail::PreparedLocalGaussianDesign design{
        sample_entries,
        0.0,
        kLocalGaussianFittingRadius
    };
    const auto sample_response_list{ CollectSampleResponses(sample_entries) };
    return design.Estimate(
        sample_response_list,
        alpha_r,
        options.thread_size,
        offset_model);
}

GroupGaussianResult EstimateGroupGaussian(
    const std::vector<GroupGaussianMemberInput> & member_list,
    double alpha_g,
    const FitOptions & options)
{
    numeric_validation::RequireFiniteNonNegative(alpha_g, "alpha_g");

    const RHBMExecutionOptions execution_options{
        .thread_size = options.thread_size
    };
    std::vector<RHBMMemberDataset> dataset_list;
    dataset_list.reserve(member_list.size());
    std::vector<RHBMBetaEstimateResult> fit_result_list;
    fit_result_list.reserve(member_list.size());
    std::vector<double> member_offset_list;
    member_offset_list.reserve(member_list.size());
    for (const auto & member : member_list)
    {
        const detail::PreparedLocalGaussianDesign design{
            member.sample_entries,
            0.0,
            kLocalGaussianFittingRadius
        };
        auto dataset{
            design.BuildDataset(
                CollectSampleResponses(member.sample_entries),
                member.local_model)
        };
        fit_result_list.emplace_back(
            rhbm_helper::EstimateBetaMDPDE(
                member.alpha_r,
                dataset,
                execution_options));
        dataset_list.emplace_back(std::move(dataset));
        member_offset_list.emplace_back(member.local_model.GetOffset());
    }
    const auto group_input{ rhbm_helper::BuildGroupInput(dataset_list, fit_result_list) };
    const auto raw_result{ rhbm_helper::EstimateGroup(alpha_g, group_input, execution_options) };
    return DecodeGroupGaussianResult(alpha_g, raw_result, member_offset_list);
}

void RunFirstStageLocalAlphaTraining(ModelObject & model_object, const FitOptions & options)
{
    model_object.EditAnalysis().EnsureSelectedAtomLocalPotentials();
    RunFirstStageLocalAlphaTraining(model_object, options, model_object.GetSelectedAtoms());
}

double detail::EstimateFirstStageAlpha(
    const LocalPotentialSampleList & samples, const FitOptions & options)
{
    std::size_t count{ 0 };
    for (const auto & sample : samples)
    {
        if (sample.point.distance < 0.0 || sample.point.distance > kLocalGaussianFittingRadius) continue;
        if (++count >= kMinimumAlphaRTrainingSampleCount) return TrainAlphaR(samples, options);
    }
    return MakeTrainingOptions(options).alpha_min;
}

void RunFirstStageLocalAlphaTraining(
    ModelObject & model_object,
    const FitOptions & options,
    const std::vector<AtomObject *> & atom_list)
{
    auto analysis{ model_object.EditAnalysis() };

    size_t count{ 0 };
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info,
            "Run local alpha training for " +
            std::to_string(atom_list.size()) + " atoms.");
    }
    for (auto * atom : atom_list)
    {
        const auto samples = AtomLocalPotentialView::For(*atom).GetRawSamplingEntries();
        analysis.SetAtomFirstStageAlphaR(*atom, detail::EstimateFirstStageAlpha(samples, options));
        count++;
        if (!options.quiet_mode)
        {
            Logger::ProgressPercent(count, atom_list.size());
        }
    }
}

} // namespace rhbm_gem::core
