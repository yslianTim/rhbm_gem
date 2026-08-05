#include <cstddef>
#include <rhbm_gem/core/GaussianEstimator.hpp>

#include "core/detail/GaussianEstimatorStages.hpp"
#include "data/detail/AtomClassifier.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/hrl/RHBMTrainer.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <algorithm>
#include <array>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace rhbm_gem::core {
namespace {
constexpr std::size_t kMinimumAlphaRTrainingSampleCount{ 10 };
constexpr std::size_t kMinimumAlphaGTrainingMemberCount{ 10 };
constexpr double kNeighborContributionDistanceMax{ 2.5 };
constexpr double kNeighborAtomSearchRange{ 2.0 * kNeighborContributionDistanceMax };
constexpr std::array<Spot, 5> kGroupPriorSummarySpotList{
    Spot::C, Spot::CA, Spot::CB, Spot::N, Spot::O
};

struct GaussianModelParameterSamples
{
    std::vector<double> amplitude_list{};
    std::vector<double> width_list{};
    std::vector<double> offset_list{};
};

using GroupMedianModelMap = std::unordered_map<GroupKey, GaussianModel3D>;

std::vector<AtomLocalPotentialEditor> BuildAtomLocalEditors(
    ModelObject & model_object,
    const std::vector<AtomObject *> & atom_list)
{
    auto analysis{ model_object.EditAnalysis() };
    std::vector<AtomLocalPotentialEditor> local_editor_list;
    local_editor_list.reserve(atom_list.size());
    for (auto * atom : atom_list)
    {
        local_editor_list.emplace_back(analysis.EnsureAtomLocalPotential(*atom));
    }
    return local_editor_list;
}

bool HasEnoughSamplesInFitRange(
    const LocalPotentialSampleList & sample_entries,
    double fit_range_min,
    double fit_range_max,
    std::size_t minimum_sample_count)
{
    std::size_t count{ 0 };
    for (const auto & sample : sample_entries)
    {
        if (sample.point.distance < fit_range_min || sample.point.distance > fit_range_max) continue;
        count++;
        if (count >= minimum_sample_count) return true;
    }
    return false;
}

RHBMExecutionOptions MakeExecutionOptions(const FitOptions & options)
{
    RHBMExecutionOptions execution_options;
    execution_options.quiet_mode = false;
    execution_options.thread_size = options.thread_size;
    return execution_options;
}

double CalculateZeroOffsetResponse(
    const LocalPotentialSample & sample,
    const GaussianModel3D & model)
{
    const auto distance{ static_cast<double>(sample.point.distance) };
    const auto model_offset{ model.ResponseAtDistance(distance) - model.SignalAtDistance(distance) };
    return static_cast<double>(sample.response) - model_offset;
}

rhbm_trainer::RHBMTrainingOptions MakeTrainingOptions(const FitOptions & options)
{
    rhbm_trainer::RHBMTrainingOptions training_options;
    training_options.execution_options = MakeExecutionOptions(options);
    return training_options;
}

std::size_t GetMinimumDatasetResponseCount(const std::vector<RHBMMemberDataset> & dataset_list)
{
    std::size_t minimum_response_count{ std::numeric_limits<std::size_t>::max() };
    for (const auto & dataset : dataset_list)
    {
        const auto response_count{ static_cast<std::size_t>(dataset.y.size()) };
        if (response_count < minimum_response_count)
        {
            minimum_response_count = response_count;
        }
    }
    return minimum_response_count;
}

LocalPotentialSampleList BuildSamplesForZeroOffsetGaussianFit(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & model)
{
    LocalPotentialSampleList adjusted_sampling_entries;
    adjusted_sampling_entries.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        const auto response{ static_cast<float>(CalculateZeroOffsetResponse(sample, model)) };
        adjusted_sampling_entries.emplace_back(LocalPotentialSample{ response, sample.point });
    }
    return adjusted_sampling_entries;
}

LocalGaussianResult DecodeLocalGaussianResult(
    double alpha_r,
    const RHBMBetaEstimateResult & fit_result,
    double offset = 0.0)
{
    const auto ols_model{
        linearization_service::DecodeParameterVector(fit_result.beta_ols).WithOffset(offset)
    };
    const auto mdpde_model{
        linearization_service::DecodeParameterVector(fit_result.beta_mdpde).WithOffset(offset)
    };
    return LocalGaussianResult{
        alpha_r,
        GaussianModel3DWithUncertainty{ ols_model, GaussianModel3DUncertainty{} },
        GaussianModel3DWithUncertainty{ mdpde_model, GaussianModel3DUncertainty{} },
        std::nullopt,
        false,
        0.0,
        fit_result
    };
}

GroupGaussianResult DecodeGroupGaussianResult(
    double alpha_g,
    const RHBMGroupEstimationResult & result,
    double offset)
{
    const auto prior{
        linearization_service::DecodeParameterVector(result.mu_prior, result.capital_lambda)
    };
    return GroupGaussianResult{
        alpha_g,
        linearization_service::DecodeParameterVector(result.mu_mean).WithOffset(offset),
        linearization_service::DecodeParameterVector(result.mu_mdpde).WithOffset(offset),
        GaussianModel3DWithUncertainty{
            prior.GetModel().WithOffset(offset),
            prior.GetStandardDeviationModel()
        }
    };
}

std::vector<LocalGaussianResult> DecodeMemberGaussianResults(
    const RHBMGroupEstimationResult & result,
    const std::vector<LocalGaussianResult> & member_result_list)
{
    const auto member_count{ static_cast<std::size_t>(result.beta_posterior_matrix.cols()) };
    if (member_result_list.size() != member_count)
    {
        throw std::invalid_argument("Group Gaussian member result count is inconsistent.");
    }
    if (result.capital_sigma_posterior_list.size() != member_count)
    {
        throw std::invalid_argument("Group Gaussian member result count is inconsistent.");
    }
    eigen_validation::RequireVectorSize(
        result.outlier_flag_array, result.beta_posterior_matrix.cols(),
        "outlier_flag_array", "Group Gaussian member result count is inconsistent.");
    eigen_validation::RequireVectorSize(
        result.statistical_distance_array, result.beta_posterior_matrix.cols(),
        "statistical_distance_array", "Group Gaussian member result count is inconsistent.");

    std::vector<LocalGaussianResult> member_results;
    member_results.reserve(member_count);
    for (Eigen::Index i = 0; i < result.beta_posterior_matrix.cols(); i++)
    {
        const auto member_index{ static_cast<std::size_t>(i) };
        const auto offset{
            member_result_list.at(member_index).mdpde.GetModel().GetOffset()
        };
        const auto gaussian{
            linearization_service::DecodeParameterVector(
                result.beta_posterior_matrix.col(i),
                result.capital_sigma_posterior_list.at(member_index))
        };
        const auto gaussian_with_offset{
            GaussianModel3DWithUncertainty{
                gaussian.GetModel().WithOffset(offset),
                gaussian.GetStandardDeviationModel()
            }
        };
        member_results.emplace_back(LocalGaussianResult{
            0.0,
            gaussian_with_offset,
            gaussian_with_offset,
            gaussian_with_offset,
            static_cast<bool>(result.outlier_flag_array(i)),
            result.statistical_distance_array(i)
        });
    }
    return member_results;
}

GroupMedianModelMap BuildGroupMedianMDPDEModelMap(const std::vector<AtomObject *> & atom_list)
{
    std::unordered_map<GroupKey, GaussianModelParameterSamples> parameter_samples_by_group;
    parameter_samples_by_group.reserve(atom_list.size());
    for (const auto * atom : atom_list)
    {
        const auto local_view{ AtomLocalPotentialView::For(*atom) };
        if (!local_view.IsAvailable()) continue;

        const auto & model{ local_view.GetEstimateMDPDE() };
        auto & parameter_samples{
            parameter_samples_by_group[data_internal::GetGroupKey(atom)]
        };
        parameter_samples.amplitude_list.emplace_back(model.GetAmplitude());
        parameter_samples.width_list.emplace_back(model.GetWidth());
        parameter_samples.offset_list.emplace_back(model.GetOffset());
    }

    GroupMedianModelMap median_model_by_group;
    median_model_by_group.reserve(parameter_samples_by_group.size());
    for (const auto & [group_key, parameter_samples] : parameter_samples_by_group)
    {
        if (parameter_samples.amplitude_list.empty()) continue;

        median_model_by_group.emplace(
            group_key,
            GaussianModel3D{
                array_helper::ComputeMedian(parameter_samples.amplitude_list),
                array_helper::ComputeMedian(parameter_samples.width_list),
                array_helper::ComputeMedian(parameter_samples.offset_list)
            });
    }
    return median_model_by_group;
}

template <typename GaussianLookup>
LocalPotentialSampleList BuildPeelingSamplingEntriesWithGaussianLookup(
    const AtomObject & atom,
    GaussianLookup lookup_gaussian)
{
    const auto local_view{ AtomLocalPotentialView::RequireFor(atom) };
    const auto raw_sampling_entries{ local_view.GetRawSamplingEntries(false) };
    const auto & neighbor_atom_list{ atom.FindNeighborAtoms(kNeighborAtomSearchRange) };
    LocalPotentialSampleList peeling_sampling_entries;
    peeling_sampling_entries.reserve(raw_sampling_entries.size());
    for (const auto & sample : raw_sampling_entries)
    {
        auto sample_position{ sample.point.position };
        auto response_value{ sample.response };
        for (const auto * neighbor_atom : neighbor_atom_list)
        {
            const auto * gaussian{ lookup_gaussian(*neighbor_atom) };
            if (gaussian == nullptr) continue;

            auto neighbor_position{ neighbor_atom->GetPosition() };
            auto distance{
                static_cast<double>(array_helper::ComputeNorm<float>(sample_position, neighbor_position))
            };
            if (distance > kNeighborContributionDistanceMax) continue;
            response_value -= static_cast<float>(gaussian->ResponseAtDistance(distance));
        }
        peeling_sampling_entries.emplace_back(LocalPotentialSample{response_value, sample.point });
    }
    return peeling_sampling_entries;
}

LocalPotentialSampleList BuildPeelingSamplingEntriesFromGroupMedianGaussian(
    const AtomObject & atom,
    const GroupMedianModelMap & median_model_by_group)
{
    return BuildPeelingSamplingEntriesWithGaussianLookup(
        atom,
        [&median_model_by_group](const AtomObject & neighbor_atom) -> const GaussianModel3D *
        {
            const auto median_model_iter{
                median_model_by_group.find(data_internal::GetGroupKey(&neighbor_atom))
            };
            if (median_model_iter != median_model_by_group.end())
            {
                return &median_model_iter->second;
            }

            const auto local_view{ AtomLocalPotentialView::For(neighbor_atom) };
            return local_view.IsAvailable() ? &local_view.GetEstimateMDPDE() : nullptr;
        });
}

LocalPotentialSampleList BuildPeelingSamplingEntriesFromFittedGroupGaussian(
    const AtomObject & atom,
    const ModelAnalysisView & analysis_view)
{
    return BuildPeelingSamplingEntriesWithGaussianLookup(
        atom,
        [&analysis_view](const AtomObject & neighbor_atom) -> const GaussianModel3D *
        {
            const auto group_key{ data_internal::GetGroupKey(&neighbor_atom) };
            if (!analysis_view.HasAtomGroup(group_key)) return nullptr;
            return &analysis_view.GetAtomGroupPrior(group_key);
        });
}

void SetPeelingSamplingEntriesFromGroupMedianGaussian(ModelObject & model_object)
{
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    auto local_editor_list{ BuildAtomLocalEditors(model_object, atom_list) };
    const auto median_model_by_group{ BuildGroupMedianMDPDEModelMap(atom_list) };
    for (size_t i = 0; i < atom_list.size(); i++)
    {
        local_editor_list[i].SetPeelingSamplingEntries(
            BuildPeelingSamplingEntriesFromGroupMedianGaussian(*atom_list[i], median_model_by_group)
        );
    }
}

void SetPeelingSamplingEntriesFromFittedGroupGaussian(ModelObject & model_object)
{
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    auto local_editor_list{ BuildAtomLocalEditors(model_object, atom_list) };
    const auto analysis_view{ model_object.GetAnalysisView() };
    for (size_t i = 0; i < atom_list.size(); i++)
    {
        local_editor_list[i].SetPeelingSamplingEntries(
            BuildPeelingSamplingEntriesFromFittedGroupGaussian(*atom_list[i], analysis_view)
        );
    }
}

std::vector<std::string> BuildGroupPriorSpotSummaryLines(const ModelObject & model_object)
{
    const auto analysis_view{ model_object.GetAnalysisView() };
    std::map<Spot, GaussianModelParameterSamples> spot_sample_map;
    for (const auto group_key : analysis_view.CollectAtomGroupKeys())
    {
        const auto & atom_list{ analysis_view.GetAtomObjectList(group_key) };
        if (atom_list.empty()) continue;

        const auto spot{ atom_list.front()->GetSpot() };
        if (std::find(
                kGroupPriorSummarySpotList.begin(), kGroupPriorSummarySpotList.end(),
                spot) == kGroupPriorSummarySpotList.end()) continue;
        const auto & prior{ analysis_view.GetAtomGroupPrior(group_key) };
        auto & sample_list{ spot_sample_map[spot] };
        sample_list.amplitude_list.emplace_back(prior.GetAmplitude());
        sample_list.width_list.emplace_back(prior.GetWidth());
        sample_list.offset_list.emplace_back(prior.GetOffset());
    }

    std::vector<std::string> summary_lines;
    summary_lines.reserve(spot_sample_map.size());
    for (const auto spot : kGroupPriorSummarySpotList)
    {
        const auto sample_iter{ spot_sample_map.find(spot) };
        if (sample_iter == spot_sample_map.end()) continue;
        const auto & sample_list{ sample_iter->second };

        const auto amplitude_mean{
            array_helper::ComputeMean(
                sample_list.amplitude_list.data(), sample_list.amplitude_list.size())
        };
        const auto width_mean{
            array_helper::ComputeMean(
                sample_list.width_list.data(), sample_list.width_list.size())
        };
        const auto offset_mean{
            array_helper::ComputeMean(
                sample_list.offset_list.data(), sample_list.offset_list.size())
        };

        std::ostringstream stream;
        stream << "| " << std::left << std::setw(8)
            << ChemicalDataHelper::GetLabel(spot)
            << " | " << std::right << std::fixed << std::setprecision(2)
            << std::setw(8) << amplitude_mean
            << " | " << std::setw(8)
            << array_helper::ComputeStandardDeviation(
                sample_list.amplitude_list.data(),
                sample_list.amplitude_list.size(),
                amplitude_mean)
            << " | " << std::setw(8) << width_mean
            << " | " << std::setw(8)
            << array_helper::ComputeStandardDeviation(
                sample_list.width_list.data(),
                sample_list.width_list.size(),
                width_mean)
            << " | " << std::setw(8) << offset_mean
            << " | " << std::setw(8)
            << array_helper::ComputeStandardDeviation(
                sample_list.offset_list.data(),
                sample_list.offset_list.size(),
                offset_mean)
            << " |";
        summary_lines.emplace_back(stream.str());
    }
    return summary_lines;
}

void LogGroupPriorSpotSummary(const ModelObject & model_object)
{
    const auto summary_lines{ BuildGroupPriorSpotSummaryLines(model_object) };
    if (summary_lines.empty())
    {
        Logger::Log(LogLevel::Info, "Group fitting prior summary by Spot: no atom groups available.");
        return;
    }

    Logger::Log(LogLevel::Info, "Group fitting prior summary by Spot:");
    Logger::Log(LogLevel::Info,
        "|---Spot---|------Amplitude------|--------Width--------|-------Offset--------|");
    Logger::Log(LogLevel::Info,
        "|          |   mean   |   s.d.   |   mean   |   s.d.   |   mean   |   s.d.   |");
    for (const auto & line : summary_lines)
    {
        Logger::Log(LogLevel::Info, line);
    }
}

void InitializeLocalFittingSeedModels(ModelObject & model_object)
{
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    auto local_editor_list{ BuildAtomLocalEditors(model_object, atom_list) };
    const auto seed_model{ GaussianModel3D{ 0.0, 1.0, 0.0 } };
    for (size_t i = 0; i < atom_list.size(); i++)
    {
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom_list[i]) };
        auto result{ local_view.GetGaussianResult() };
        result.ols = GaussianModel3DWithUncertainty{
            seed_model,
            GaussianModel3DUncertainty{}
        };
        result.mdpde = GaussianModel3DWithUncertainty{
            seed_model,
            GaussianModel3DUncertainty{}
        };
        result.posterior.reset();
        result.is_outlier = false;
        result.statistical_distance = 0.0;
        result.fit_result.reset();
        local_editor_list[i].SetGaussianResult(result);
    }
}

} // namespace

void RunFixedOffsetLocalFitting(
    ModelObject & model_object,
    const FitOptions & options,
    LocalFittingPass pass)
{
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    const auto selected_atom_size{ atom_list.size() };
    auto local_editor_list{ BuildAtomLocalEditors(model_object, atom_list) };
    const auto analysis_view{ model_object.GetAnalysisView() };
    const auto stage_label{
        pass == LocalFittingPass::FirstStage ? "1st-stage" : "3rd-stage"
    };
    size_t atom_count{ 0 };
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info,
            "Run " + std::string{ stage_label } + " local atom fitting for " +
            std::to_string(selected_atom_size) + " atoms.");
    }

#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(options.thread_size)
#endif
    for (size_t i = 0; i < selected_atom_size; i++)
    {
        auto & atom{ *atom_list[i] };
        const auto local_view{ AtomLocalPotentialView::RequireFor(atom) };
        LocalPotentialSampleList sample_entries;
        GaussianModel3D offset_model;
        switch (pass)
        {
        case LocalFittingPass::FirstStage:
            sample_entries = local_view.GetRawSamplingEntries();
            offset_model = local_view.GetGaussianResult().mdpde.GetModel();
            break;
        case LocalFittingPass::ThirdStage:
            sample_entries = local_view.GetPeelingSamplingEntries(false);
            offset_model = analysis_view.GetAtomGroupPrior(data_internal::GetGroupKey(&atom));
            break;
        }

        auto result{
            EstimateLocalGaussian(sample_entries, local_view.GetAlphaR(), options, offset_model)
        };
        local_editor_list[i].SetGaussianResult(result);

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
}

double TrainAlphaR(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const FitOptions & options)
{
    numeric_validation::RequireFiniteNonNegativeRange(
        options.distance_min, options.distance_max, "fit range");

    std::vector<RHBMMemberDataset> dataset_list;
    dataset_list.reserve(sample_entries_list.size());
    for (const auto & sample_entries : sample_entries_list)
    {
        dataset_list.emplace_back(
            rhbm_helper::BuildMemberDataset(
                sample_entries,
                options.distance_min,
                options.distance_max));
    }
    auto training_options{ MakeTrainingOptions(options) };
    if (!dataset_list.empty())
    {
        const auto minimum_response_count{ GetMinimumDatasetResponseCount(dataset_list) };
        if (minimum_response_count < 2)
        {
            return training_options.alpha_min;
        }
        if (training_options.subset_size > minimum_response_count)
        {
            training_options.subset_size = minimum_response_count;
        }
    }
    return rhbm_trainer::CrossValidationAlphaR(dataset_list, training_options).best_alpha;
}

double TrainAlphaG(
    const std::vector<std::vector<LocalGaussianResult>> & member_result_list,
    const FitOptions & options)
{
    std::vector<std::vector<RHBMParameterVector>> beta_group_list;
    beta_group_list.reserve(member_result_list.size());
    for (const auto & member_results : member_result_list)
    {
        std::vector<RHBMParameterVector> beta_list;
        beta_list.reserve(member_results.size());
        for (const auto & member_result : member_results)
        {
            beta_list.emplace_back(
                linearization_service::EncodeGaussianToParameterVector(member_result.mdpde.GetModel()));
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
    auto range_min{ options.distance_min };
    auto range_max{ options.distance_max };
    numeric_validation::RequireFiniteNonNegativeRange(range_min, range_max, "fit range");
    numeric_validation::RequireFiniteNonNegative(alpha_r, "alpha_r");
    numeric_validation::RequireFinite(offset_model.GetOffset(), "offset");

    auto execution_options{ MakeExecutionOptions(options) };
    const auto adjusted_sampling_entries{
        BuildSamplesForZeroOffsetGaussianFit(sample_entries, offset_model)
    };
    auto dataset{
        rhbm_helper::BuildMemberDataset(adjusted_sampling_entries, range_min, range_max)
    };
    const auto result{ rhbm_helper::EstimateBetaMDPDE(alpha_r, dataset, execution_options) };
    return DecodeLocalGaussianResult(alpha_r, result, offset_model.GetOffset());
}

GroupGaussianResult EstimateGroupGaussian(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const std::vector<LocalGaussianResult> & member_result_list,
    double alpha_g,
    const FitOptions & options)
{
    numeric_validation::RequireFiniteNonNegativeRange(
        options.distance_min, options.distance_max, "fit range");
    numeric_validation::RequireFiniteNonNegative(alpha_g, "alpha_g");

    if (sample_entries_list.size() != member_result_list.size())
    {
        throw std::invalid_argument("sample_entries_list and member_result_list sizes are inconsistent.");
    }

    auto execution_options{ MakeExecutionOptions(options) };
    const auto range_min{ options.distance_min };
    const auto range_max{ options.distance_max };
    std::vector<RHBMMemberDataset> dataset_list;
    dataset_list.reserve(sample_entries_list.size());
    for (std::size_t i = 0; i < sample_entries_list.size(); i++)
    {
        const auto sampling_entries{
            BuildSamplesForZeroOffsetGaussianFit(
                sample_entries_list.at(i),
                member_result_list.at(i).mdpde.GetModel())
        };
        dataset_list.emplace_back(
            rhbm_helper::BuildMemberDataset(sampling_entries, range_min, range_max));
    }
    std::vector<RHBMBetaEstimateResult> fit_result_list;
    fit_result_list.reserve(member_result_list.size());
    for (std::size_t i = 0; i < member_result_list.size(); i++)
    {
        fit_result_list.emplace_back(
            rhbm_helper::EstimateBetaMDPDE(
                member_result_list.at(i).alpha_r,
                dataset_list.at(i),
                execution_options));
    }
    const auto group_input{ rhbm_helper::BuildGroupInput(dataset_list, fit_result_list) };
    const auto raw_result{ rhbm_helper::EstimateGroup(alpha_g, group_input, execution_options) };
    std::vector<double> member_offset_list;
    member_offset_list.reserve(member_result_list.size());
    for (const auto & member_result : member_result_list)
    {
        member_offset_list.emplace_back(member_result.mdpde.GetModel().GetOffset());
    }
    const auto group_offset{ array_helper::ComputeMedian(member_offset_list) };
    auto result{ DecodeGroupGaussianResult(alpha_g, raw_result, group_offset) };
    result.member_results = DecodeMemberGaussianResults(raw_result, member_result_list);
    return result;
}

void RunLocalAlphaTraining(
    ModelObject & model_object,
    const FitOptions & options,
    LocalFittingPass pass)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };
    const auto group_key_list{ analysis_view.CollectAtomGroupKeys() };

    size_t count{ 0 };
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run local alpha training for " + std::to_string(group_key_list.size()) + " groups.");
    }
    for (const auto group_key : group_key_list)
    {
        const auto & group_atom_list{ analysis_view.GetAtomObjectList(group_key) };
        std::vector<LocalPotentialSampleList> sample_entries_list;
        sample_entries_list.reserve(group_atom_list.size());
        for (auto * atom : group_atom_list)
        {
            analysis.EnsureAtomLocalPotential(*atom);
            const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
            auto sample_entries{
                pass == LocalFittingPass::FirstStage
                    ? local_view.GetRawSamplingEntries()
                    : local_view.GetPeelingSamplingEntries(false)
            };
            if (!HasEnoughSamplesInFitRange(
                    sample_entries,
                    options.distance_min,
                    options.distance_max,
                    kMinimumAlphaRTrainingSampleCount)) continue;
            sample_entries_list.emplace_back(std::move(sample_entries));
        }
        if (!sample_entries_list.empty())
        {
            const auto alpha_r{ TrainAlphaR(sample_entries_list, options) };
            for (auto * atom : group_atom_list)
            {
                analysis.EnsureAtomLocalPotential(*atom).SetAlphaR(alpha_r);
            }
        }
        count++;
        if (!options.quiet_mode)
        {
            Logger::ProgressPercent(count, group_key_list.size());
        }
    }
}

namespace {

void RunGroupAlphaTraining(ModelObject & model_object, const FitOptions & options)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };
    const auto group_key_list{ analysis_view.CollectAtomGroupKeys() };

    std::vector<std::vector<LocalGaussianResult>> member_result_list;
    member_result_list.reserve(group_key_list.size());
    for (const auto group_key : group_key_list)
    {
        const auto & group_atom_list{ analysis_view.GetAtomObjectList(group_key) };
        if (group_atom_list.size() < kMinimumAlphaGTrainingMemberCount) continue;
        if (group_atom_list.front()->IsMainChainAtom() == false) continue;

        std::vector<LocalGaussianResult> group_member_results;
        group_member_results.reserve(group_atom_list.size());
        for (auto * atom : group_atom_list)
        {
            analysis.EnsureAtomLocalPotential(*atom);
            const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
            group_member_results.emplace_back(local_view.GetGaussianResult());
        }
        member_result_list.emplace_back(std::move(group_member_results));
    }

    const auto alpha_g{ TrainAlphaG(member_result_list, options) };
    for (const auto group_key : group_key_list)
    {
        analysis.SetAtomGroupAlphaG(group_key, alpha_g);
    }
}

void RunGroupPotentialFitting(
    ModelObject & model_object,
    const FitOptions & options,
    bool apply_selection = false,
    bool use_peeling_sampling_entries = true)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };
    const auto & selected_atom_list{ model_object.GetSelectedAtoms() };
    for (auto * atom : selected_atom_list)
    {
        analysis.EnsureAtomLocalPotential(*atom);
    }
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run atom group fitting.");
    }

    auto group_key_list{ analysis_view.CollectAtomGroupKeys() };
    auto group_key_size{ group_key_list.size() };
    size_t key_count{ 0 };

#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(options.thread_size)
#endif
    for (size_t k = 0; k < group_key_size; k++)
    {
        auto group_key{ group_key_list[k] };
        const auto & atom_list{ analysis_view.GetAtomObjectList(group_key) };
        const auto alpha_g{ analysis_view.GetAtomAlphaG(group_key) };
        std::vector<LocalPotentialSampleList> sample_entries_list;
        std::vector<LocalGaussianResult> member_result_list;
        sample_entries_list.reserve(atom_list.size());
        member_result_list.reserve(atom_list.size());
        for (const auto & atom : atom_list)
        {
            const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
            sample_entries_list.emplace_back(
                use_peeling_sampling_entries
                    ? local_view.GetPeelingSamplingEntries(apply_selection)
                    : local_view.GetRawSamplingEntries(apply_selection));
            member_result_list.emplace_back(local_view.GetGaussianResult());
        }
        const auto result{
            EstimateGroupGaussian(sample_entries_list, member_result_list, alpha_g, options)
        };

#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            analysis.ApplyAtomGroupGaussianResult(group_key, result);
            key_count++;
            if (!options.quiet_mode)
            {
                Logger::ProgressBar(key_count, group_key_size);
            }
        }
    }
}

} // namespace

void RunPotentialFittingWorkflow(ModelObject & model_object, const FitOptions & options)
{
    RunLocalAlphaTraining(model_object, options, LocalFittingPass::FirstStage);

    InitializeLocalFittingSeedModels(model_object);
    RunFixedOffsetLocalFitting(model_object, options, LocalFittingPass::FirstStage);
    RunGroupAlphaTraining(model_object, options);
    RunGroupPotentialFitting(model_object, options, true, false);

    RunSecondStageLocalFitting(model_object, options);

    RunGroupAlphaTraining(model_object, options);
    SetPeelingSamplingEntriesFromGroupMedianGaussian(model_object);
    RunGroupPotentialFitting(model_object, options);

    SetPeelingSamplingEntriesFromFittedGroupGaussian(model_object);
    RunLocalAlphaTraining(model_object, options, LocalFittingPass::ThirdStage);
    RunFixedOffsetLocalFitting(model_object, options, LocalFittingPass::ThirdStage);

    RunGroupAlphaTraining(model_object, options);
    RunGroupPotentialFitting(model_object, options);
    if (!options.quiet_mode)
    {
        LogGroupPriorSpotSummary(model_object);
    }
}

} // namespace rhbm_gem::core
