#include <cstddef>
#include <rhbm_gem/core/GaussianEstimator.hpp>

#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/SampleFilter.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/hrl/RHBMTrainer.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <algorithm>
#include <atomic>
#include <iomanip>
#include <limits>
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
constexpr double kResidualInterceptRangeMin{ 1.0 };
constexpr double kResidualInterceptRangeMax{ 1.5 };
constexpr double kInterceptDampingWeight{ 0.5 };
constexpr std::size_t kInterceptCycleHistorySize{ 64 };

std::vector<AtomLocalPotentialEditor> BuildSelectedAtomLocalEditors(ModelObject & model_object)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto & atom_list{ model_object.GetSelectedAtoms() };
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
    execution_options.quiet_mode = true;
    execution_options.thread_size = options.thread_size;
    return execution_options;
}

double EstimateInitialIntercept(const LocalPotentialSampleList & sample_entries)
{
    float maximum_distance{ 0.0f };
    for (const auto & sample : sample_entries)
    {
        if (sample.point.distance > maximum_distance) maximum_distance = sample.point.distance;
    }
    //Logger::Log(LogLevel::Info, "Estimated maximum distance: "+ std::to_string(maximum_distance));

    std::vector<double> response_list;
    response_list.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        if (sample.point.distance != maximum_distance) continue;
        const auto response{ static_cast<double>(sample.response) };
        numeric_validation::RequireFinite(response, "maximum distance shell response");
        response_list.emplace_back(response);
    }
    if (response_list.empty()) return 0.0;

    const auto lowest_response_list{
        array_helper::ComputeSmallestProportionValues(response_list, 0.10)
    };
    return array_helper::ComputeMedian(lowest_response_list);
}

LocalPotentialSampleList BuildResidualSampleEntries(
    const LocalPotentialSampleList & sample_entries,
    const RHBMParameterVector & beta,
    double range_min,
    double range_max)
{
    const auto signal_model{ linearization_service::DecodeParameterVector(beta) };
    const auto median_sample_entries{
        sample_filter::BuildMedianResponseSampleEntriesByRadius(sample_entries)
    };
    LocalPotentialSampleList residual_sample_entries;
    residual_sample_entries.reserve(median_sample_entries.size());
    for (const auto & sample : median_sample_entries)
    {
        const auto distance{ sample.point.distance };
        if (distance < static_cast<float>(range_min)) continue;
        if (distance > static_cast<float>(range_max)) continue;

        const auto residual{
            static_cast<double>(sample.response) -
                signal_model.ResponseAtDistance(static_cast<double>(distance))
        };
        residual_sample_entries.emplace_back(
            LocalPotentialSample{ static_cast<float>(residual), sample.point }
        );
    }
    return residual_sample_entries;
}

double EstimateResidualIntercept(
    const LocalPotentialSampleList & sample_entries,
    const RHBMBetaEstimateResult & fit_result)
{
    const auto residual_sample_entries{
        BuildResidualSampleEntries(
            sample_entries,
            fit_result.beta_mdpde,
            kResidualInterceptRangeMin,
            kResidualInterceptRangeMax)
    };
    std::vector<double> residual_list;
    residual_list.reserve(residual_sample_entries.size());
    for (const auto & sample : residual_sample_entries)
    {
        residual_list.emplace_back(static_cast<double>(sample.response));
    }
    if (residual_list.empty()) return 0.0;
    //return array_helper::ComputeMean(residual_list.data(), residual_list.size());
    return array_helper::ComputeMedian(residual_list);
}

rhbm_trainer::RHBMTrainingOptions MakeTrainingOptions(const FitOptions & options)
{
    rhbm_trainer::RHBMTrainingOptions training_options;
    training_options.execution_options = MakeExecutionOptions(options);
    return training_options;
}

std::vector<RHBMMemberDataset> BuildMemberDatasetList(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const FitOptions & options)
{
    std::vector<RHBMMemberDataset> dataset_list;
    dataset_list.reserve(sample_entries_list.size());
    for (const auto & sample_entries : sample_entries_list)
    {
        dataset_list.emplace_back(
            rhbm_helper::BuildMemberDataset(
                sample_entries,
                options.distance_min,
                options.distance_max)
        );
    }
    return dataset_list;
}

std::vector<RHBMMemberDataset> BuildMemberDatasetList(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const std::vector<LocalGaussianResult> & member_result_list,
    const FitOptions & options)
{
    if (sample_entries_list.size() != member_result_list.size())
    {
        throw std::invalid_argument("sample_entries_list and member_result_list sizes are inconsistent.");
    }

    std::vector<RHBMMemberDataset> dataset_list;
    dataset_list.reserve(sample_entries_list.size());
    for (std::size_t i = 0; i < sample_entries_list.size(); i++)
    {
        const auto intercept{ member_result_list.at(i).mdpde.GetModel().GetIntercept() };
        const auto shifted_sample_entries{
            sample_filter::BuildResponseShiftedSampleEntries(sample_entries_list.at(i), intercept)
        };
        dataset_list.emplace_back(
            rhbm_helper::BuildMemberDataset(
                shifted_sample_entries,
                options.distance_min,
                options.distance_max)
        );
    }
    return dataset_list;
}

LocalGaussianResult DecodeLocalGaussianResult(
    double alpha_r,
    const RHBMBetaEstimateResult & fit_result,
    double intercept = 0.0)
{
    const auto ols_model{
        linearization_service::DecodeParameterVector(fit_result.beta_ols)
            .WithIntercept(intercept)
    };
    const auto mdpde_model{
        linearization_service::DecodeParameterVector(fit_result.beta_mdpde)
            .WithIntercept(intercept)
    };
    return LocalGaussianResult{
        alpha_r,
        GaussianModel3DWithUncertainty{ ols_model, GaussianModel3DUncertainty{} },
        GaussianModel3DWithUncertainty{ mdpde_model, GaussianModel3DUncertainty{} },
        false,
        0.0,
        fit_result
    };
}

GaussianModel3DWithUncertainty WithModelIntercept(
    const GaussianModel3DWithUncertainty & gaussian,
    double intercept)
{
    return GaussianModel3DWithUncertainty{
        gaussian.GetModel().WithIntercept(intercept),
        gaussian.GetStandardDeviationModel()
    };
}

std::vector<double> ExtractMemberInterceptList(
    const std::vector<LocalGaussianResult> & member_result_list)
{
    std::vector<double> intercept_list;
    intercept_list.reserve(member_result_list.size());
    for (const auto & member_result : member_result_list)
    {
        intercept_list.emplace_back(member_result.mdpde.GetModel().GetIntercept());
    }
    return intercept_list;
}

GroupGaussianResult DecodeGroupGaussianResult(
    double alpha_g,
    const RHBMGroupEstimationResult & result,
    double intercept)
{
    return GroupGaussianResult{
        alpha_g,
        linearization_service::DecodeParameterVector(result.mu_mean).WithIntercept(intercept),
        linearization_service::DecodeParameterVector(result.mu_mdpde).WithIntercept(intercept),
        WithModelIntercept(
            linearization_service::DecodeParameterVector(result.mu_prior, result.capital_lambda),
            intercept)
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
        const auto intercept{
            member_result_list.at(member_index).mdpde.GetModel().GetIntercept()
        };
        const auto gaussian{
            WithModelIntercept(
                linearization_service::DecodeParameterVector(
                    result.beta_posterior_matrix.col(i),
                    result.capital_sigma_posterior_list.at(member_index)),
                intercept)
        };
        member_results.emplace_back(LocalGaussianResult{
            0.0,
            gaussian,
            gaussian,
            static_cast<bool>(result.outlier_flag_array(i)),
            result.statistical_distance_array(i)
        });
    }
    return member_results;
}

std::vector<RHBMBetaEstimateResult> BuildMemberFitResultList(
    const std::vector<LocalGaussianResult> & member_result_list)
{
    std::vector<RHBMBetaEstimateResult> fit_result_list;
    fit_result_list.reserve(member_result_list.size());
    for (const auto & member_result : member_result_list)
    {
        if (!member_result.fit_result.has_value())
        {
            throw std::invalid_argument(
                "member_result_list contains a result without transient fit state.");
        }
        fit_result_list.emplace_back(*member_result.fit_result);
    }
    return fit_result_list;
}

using FittedGaussianSnapshot = std::unordered_map<const AtomObject *, GaussianModel3D>;

FittedGaussianSnapshot BuildFittedGaussianSnapshot(const std::vector<AtomObject *> & atom_list)
{
    FittedGaussianSnapshot snapshot;
    snapshot.reserve(atom_list.size());
    for (const auto * atom : atom_list)
    {
        const auto local_view{ AtomLocalPotentialView::For(*atom) };
        if (!local_view.IsAvailable()) continue;

        const auto & gaussian_result{ local_view.GetGaussianResult() };
        if (!gaussian_result.fit_result.has_value()) continue;

        snapshot.emplace(atom, gaussian_result.mdpde.GetModel());
    }
    return snapshot;
}

int GetLocalFittingSpotSortGroup(Spot spot)
{
    switch (spot)
    {
    case Spot::O:
        return 0;
    case Spot::C:
        return 1;
    case Spot::CA:
        return 2;
    case Spot::N:
        return 3;
    default:
        return 4;
    }
}

int GetLocalFittingSpotSortValue(Spot spot)
{
    if (GetLocalFittingSpotSortGroup(spot) < 4)
    {
        return 0;
    }
    return static_cast<int>(spot);
}

std::vector<size_t> BuildLocalFittingUpdateOrder(const std::vector<AtomObject *> & atom_list)
{
    std::vector<size_t> update_order;
    update_order.reserve(atom_list.size());
    for (size_t i = 0; i < atom_list.size(); i++)
    {
        update_order.emplace_back(i);
    }

    std::stable_sort(
        update_order.begin(),
        update_order.end(),
        [&](size_t lhs, size_t rhs)
        {
            const auto * lhs_atom{ atom_list.at(lhs) };
            const auto * rhs_atom{ atom_list.at(rhs) };
            if (lhs_atom->GetSequenceID() != rhs_atom->GetSequenceID())
            {
                return lhs_atom->GetSequenceID() < rhs_atom->GetSequenceID();
            }

            const auto lhs_spot{ lhs_atom->GetSpot() };
            const auto rhs_spot{ rhs_atom->GetSpot() };
            const auto lhs_group{ GetLocalFittingSpotSortGroup(lhs_spot) };
            const auto rhs_group{ GetLocalFittingSpotSortGroup(rhs_spot) };
            if (lhs_group != rhs_group)
            {
                return lhs_group < rhs_group;
            }
            return GetLocalFittingSpotSortValue(lhs_spot) <
                GetLocalFittingSpotSortValue(rhs_spot);
        });
    return update_order;
}

LocalPotentialSampleList UpdateSampleListWithFittedGaussian(
    const AtomObject & atom,
    const FittedGaussianSnapshot & snapshot)
{
    const auto local_view{ AtomLocalPotentialView::RequireFor(atom) };
    const auto sample_entries{ local_view.GetSamplingEntries(false) };
    const auto & neighbor_atom_list{ atom.FindNeighborAtoms() };
    LocalPotentialSampleList updated_list;
    updated_list.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        auto sample_position{ sample.point.position };
        auto response_value{ sample.response };
        for (const auto * neighbor_atom : neighbor_atom_list)
        {
            const auto gaussian_iter{ snapshot.find(neighbor_atom) };
            if (gaussian_iter == snapshot.end()) continue;

            auto neighbor_position{ neighbor_atom->GetPosition() };
            auto distance{
                static_cast<double>(
                    array_helper::ComputeNorm<float>(sample_position, neighbor_position))
            };
            response_value -= static_cast<float>(gaussian_iter->second.SignalAtDistance(distance));
        }
        updated_list.emplace_back(LocalPotentialSample{response_value, sample.point });
    }
    return updated_list;
}

} // namespace

double TrainAlphaR(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const FitOptions & options)
{
    numeric_validation::RequireFiniteNonNegativeRange(
        options.distance_min, options.distance_max, "fit range");

    const auto dataset_list{ BuildMemberDatasetList(sample_entries_list, options) };
    const auto training_options{ MakeTrainingOptions(options) };
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
    double intercept)
{
    auto range_min{ options.distance_min };
    auto range_max{ options.distance_max };
    numeric_validation::RequireFiniteNonNegativeRange(range_min, range_max, "fit range");
    numeric_validation::RequireFiniteNonNegative(alpha_r, "alpha_r");
    numeric_validation::RequireFinite(intercept, "intercept");
    auto execution_options{ MakeExecutionOptions(options) };
    auto shifted_sample_entries{
        sample_filter::BuildResponseShiftedSampleEntries(sample_entries, intercept)
    };
    auto dataset{
        rhbm_helper::BuildMemberDataset(shifted_sample_entries, range_min, range_max)
    };
    const auto result{ rhbm_helper::EstimateBetaMDPDE(alpha_r, dataset, execution_options) };
    return DecodeLocalGaussianResult(alpha_r, result, intercept);
}

LocalGaussianResult EstimateLocalGaussianWithIntercept(
    const LocalPotentialSampleList & sample_entries,
    double alpha_r,
    const FitOptions & options,
    double intercept_initial)
{
    auto range_min{ options.distance_min };
    auto range_max{ options.distance_max };
    numeric_validation::RequireFiniteNonNegativeRange(range_min, range_max, "fit range");
    numeric_validation::RequireFiniteNonNegative(alpha_r, "alpha_r");

    auto execution_options{ MakeExecutionOptions(options) };
    double intercept{ intercept_initial };
    LocalGaussianResult result;
    std::vector<double> intercept_history;
    intercept_history.reserve(kInterceptCycleHistorySize);
    double best_intercept{ intercept };
    double best_defect{ std::numeric_limits<double>::infinity() };
    auto max_iterations{ execution_options.max_iterations };
    auto tolerance{ execution_options.tolerance };
    for (int t = 0; t < max_iterations; t++)
    {
        result = EstimateLocalGaussian(sample_entries, alpha_r, options, intercept);
        if (!intercept_history.empty() &&
            std::abs(intercept - intercept_history.back()) < tolerance)
        {
            break;
        }

        bool cycle_detected{ false };
        for (std::size_t period = 2; period <= intercept_history.size(); period++)
        {
            const auto cycle_begin{ intercept_history.size() - period };
            if (std::abs(intercept - intercept_history[cycle_begin]) >= tolerance) continue;

            double final_intercept{ 0.0 };
            for (std::size_t i = cycle_begin; i < intercept_history.size(); i++)
            {
                final_intercept += intercept_history[i];
            }
            final_intercept /= static_cast<double>(period);
            result = EstimateLocalGaussian(sample_entries, alpha_r, options, final_intercept);
            Logger::Log(LogLevel::Debug,
                "Cycle detected in local Gaussian intercept estimation with period " +
                std::to_string(period) + "; refitting at cycle mean.");
            cycle_detected = true;
            break;
        }
        if (cycle_detected) break;

        const auto raw_intercept{
            EstimateResidualIntercept(sample_entries, *result.fit_result)
        };
        const auto defect{ std::abs(raw_intercept - intercept) };
        if (defect < best_defect)
        {
            best_intercept = intercept;
            best_defect = defect;
        }
        if (t + 1 == max_iterations)
        {
            result = EstimateLocalGaussian(sample_entries, alpha_r, options, best_intercept);
            Logger::Log(LogLevel::Warning,
                "Maximum iterations reached in local Gaussian estimation with intercept; "
                "refitting at best fixed-point candidate with defect = " +
                std::to_string(best_defect) + ".");
            break;
        }

        if (intercept_history.size() == kInterceptCycleHistorySize)
        {
            intercept_history.erase(intercept_history.begin());
        }
        intercept_history.emplace_back(intercept);
        intercept += kInterceptDampingWeight * (raw_intercept - intercept);
    }
    //Logger::Log(LogLevel::Info, "Estimated intercept: " + std::to_string(intercept));
    return result;
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
    const auto dataset_list{ BuildMemberDatasetList(sample_entries_list, member_result_list, options) };
    const auto fit_result_list{ BuildMemberFitResultList(member_result_list) };
    const auto group_input{ rhbm_helper::BuildGroupInput(dataset_list, fit_result_list) };
    const auto raw_result{ rhbm_helper::EstimateGroup(alpha_g, group_input, execution_options) };
    const auto member_intercept_list{ ExtractMemberInterceptList(member_result_list) };
    const auto group_intercept{ array_helper::ComputeMedian(member_intercept_list) };
    auto result{ DecodeGroupGaussianResult(alpha_g, raw_result, group_intercept) };
    result.member_results = DecodeMemberGaussianResults(raw_result, member_result_list);
    return result;
}

void RunLocalAlphaTraining(ModelObject & model_object, const FitOptions & options)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };
    const auto class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };
    const auto group_key_list{ analysis_view.CollectAtomGroupKeys(class_key) };

    size_t count{ 0 };
    Logger::Log(LogLevel::Info, "Run local alpha training for " + std::to_string(group_key_list.size()) + " groups.");
    for (const auto group_key : group_key_list)
    {
        const auto & group_atom_list{
            analysis_view.GetAtomObjectList(group_key, class_key)
        };
        std::vector<LocalPotentialSampleList> sample_entries_list;
        sample_entries_list.reserve(group_atom_list.size());
        for (auto * atom : group_atom_list)
        {
            analysis.EnsureAtomLocalPotential(*atom);
            const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
            const auto & sample_entries{ local_view.GetSamplingEntries() };
            if (!HasEnoughSamplesInFitRange(
                    sample_entries,
                    options.distance_min,
                    options.distance_max,
                    kMinimumAlphaRTrainingSampleCount)) continue;
            sample_entries_list.emplace_back(sample_entries);
        }
        sample_entries_list.shrink_to_fit();
        if (!sample_entries_list.empty())
        {
            const auto alpha_r{ TrainAlphaR(sample_entries_list, options) };
            for (auto * atom : group_atom_list)
            {
                analysis.EnsureAtomLocalPotential(*atom).SetAlphaR(alpha_r);
            }
        }
        Logger::ProgressPercent(++count, group_key_list.size());
    }
}

void RunGroupAlphaTraining(ModelObject & model_object, const FitOptions & options)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };
    const auto class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };
    const auto group_key_list{ analysis_view.CollectAtomGroupKeys(class_key) };

    std::vector<std::vector<LocalGaussianResult>> member_result_list;
    member_result_list.reserve(group_key_list.size());
    for (const auto group_key : group_key_list)
    {
        const auto & group_atom_list{
            analysis_view.GetAtomObjectList(group_key, class_key)
        };
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
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto & class_key_tmp{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        for (const auto group_key : analysis_view.CollectAtomGroupKeys(class_key_tmp))
        {
            analysis.SetAtomGroupAlphaG(group_key, class_key_tmp, alpha_g);
        }
    }
}

void RunLocalPotentialFitting(ModelObject & model_object, const FitOptions & options)
{
    const auto selected_atom_size{ model_object.GetSelectedAtomCount() };
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    auto local_editor_list{ BuildSelectedAtomLocalEditors(model_object) };
    std::atomic<size_t> atom_count{ 0 };
    Logger::Log(LogLevel::Info,
        "Run local atom fitting for " + std::to_string(selected_atom_size) + " atoms.");

#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(options.thread_size)
#endif
    for (size_t i = 0; i < selected_atom_size; i++)
    {
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom_list[i]) };
        auto sample_entries{ local_view.GetSamplingEntries() };
        auto intercept_initial{ EstimateInitialIntercept(sample_entries) };
        const auto result{
            EstimateLocalGaussianWithIntercept(
                sample_entries, local_view.GetAlphaR(), options, intercept_initial)
        };
        local_editor_list[i].SetGaussianResult(result);

#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            atom_count++;
            Logger::ProgressPercent(atom_count, selected_atom_size);
        }
    }

    const size_t maximum_iter_size{ 100 };
    constexpr double convergence_tolerance{ 1.0e-5 };
    constexpr double convergence_percentile{ 0.90 };
    constexpr std::size_t required_consecutive_convergence_count{ 2 };
    const auto update_order{ BuildLocalFittingUpdateOrder(atom_list) };
    std::vector<LocalPotentialSampleList> sample_entries_list(selected_atom_size);
    std::vector<Eigen::VectorXd> previous_estimation_list(selected_atom_size);
    std::size_t consecutive_convergence_count{ 0 };
    for (size_t i = 0; i < selected_atom_size; i++)
    {
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom_list[i]) };
        previous_estimation_list[i] = local_view.GetGaussianResult().mdpde.GetModel().ToVector();
    }

    Logger::Log(LogLevel::Info, "Run updated local atom fitting with iterations...");
    auto snapshot{ BuildFittedGaussianSnapshot(atom_list) };
    for (size_t iter = 0; iter < maximum_iter_size; iter++)
    {
        std::vector<Eigen::VectorXd> current_estimation_list(selected_atom_size);
        for (const auto i : update_order)
        {
            const auto & atom{ *atom_list[i] };
            const auto local_view{ AtomLocalPotentialView::RequireFor(atom) };
            auto sample_entries{
                UpdateSampleListWithFittedGaussian(atom, snapshot)
            };
            const auto intercept{ snapshot.at(&atom).GetIntercept() };
            const auto result{
                EstimateLocalGaussianWithIntercept(sample_entries, local_view.GetAlphaR(), options, intercept)
            };
            sample_entries_list[i] = std::move(sample_entries);
            current_estimation_list[i] = result.mdpde.GetModel().ToVector();
            local_editor_list[i].SetGaussianResult(result);
            snapshot.at(&atom) = result.mdpde.GetModel();
        }
        std::vector<double> parameter_change_list(selected_atom_size);
        double maximum_parameter_change{ 0.0 };
        for (size_t i = 0; i < selected_atom_size; i++)
        {
            const auto parameter_change{
                (current_estimation_list[i] - previous_estimation_list[i]).norm()
            };
            parameter_change_list[i] = parameter_change;
            if (parameter_change > maximum_parameter_change)
            {
                maximum_parameter_change = parameter_change;
            }
        }

        const auto parameter_change_percentile90{
            array_helper::ComputePercentile(parameter_change_list, convergence_percentile)
        };
        if (maximum_parameter_change * maximum_parameter_change < convergence_tolerance)
        {
            consecutive_convergence_count++;
        }
        else
        {
            consecutive_convergence_count = 0;
        }

        std::ostringstream progress_message;
        progress_message << "Local fitting iteration " << iter + 1 << '/'
            << maximum_iter_size << ", max parameter change = "
            << std::fixed << std::setprecision(6) << maximum_parameter_change
            << ", 90th percentile parameter change = "
            << parameter_change_percentile90
            << ", stable streak = " << consecutive_convergence_count << '/'
            << required_consecutive_convergence_count;
        Logger::ProgressLine(progress_message.str());

        if (consecutive_convergence_count >= required_consecutive_convergence_count)
        {
            Logger::FinishProgressLine();
            Logger::Log(LogLevel::Info,
                "Converged after " + std::to_string(iter + 1) +
                " iterations with max parameter change = " +
                std::to_string(maximum_parameter_change) +
                " and 90th percentile parameter change = " +
                std::to_string(parameter_change_percentile90) + ".");
            break;
        }
        previous_estimation_list = std::move(current_estimation_list);
        if (iter == maximum_iter_size - 1)
        {
            Logger::FinishProgressLine();
            Logger::Log(LogLevel::Warning,
                "Reached maximum iteration size with max parameter change = " +
                std::to_string(maximum_parameter_change) +
                " and 90th percentile parameter change = " +
                std::to_string(parameter_change_percentile90));
        }
    }

    for (size_t i = 0; i < selected_atom_size; i++)
    {
        local_editor_list[i].SetSamplingEntries(std::move(sample_entries_list[i]));
    }
}

void RunGroupPotentialFitting(ModelObject & model_object, const FitOptions & options)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };
    for (auto * atom : model_object.GetSelectedAtoms())
    {
        analysis.EnsureAtomLocalPotential(*atom);
    }
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        Logger::Log(LogLevel::Info, "Class type: " + class_key);

        auto group_key_list{ analysis_view.CollectAtomGroupKeys(class_key) };
        auto group_key_size{ group_key_list.size() };
        std::atomic<size_t> key_count{ 0 };

#ifdef USE_OPENMP
        #pragma omp parallel for num_threads(options.thread_size)
#endif
        for (size_t k = 0; k < group_key_size; k++)
        {
            auto group_key{ group_key_list[k] };
            const auto & atom_list{ analysis_view.GetAtomObjectList(group_key, class_key) };
            std::vector<LocalPotentialSampleList> member_sample_list;
            std::vector<LocalGaussianResult> member_result_list;
            member_sample_list.reserve(atom_list.size());
            member_result_list.reserve(atom_list.size());
            for (const auto & atom : atom_list)
            {
                const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
                member_sample_list.emplace_back(local_view.GetSamplingEntries(false));
                member_result_list.emplace_back(local_view.GetGaussianResult());
            }
            const auto result{
                EstimateGroupGaussian(
                    member_sample_list,
                    member_result_list,
                    analysis_view.GetAtomAlphaG(group_key, class_key),
                    options)
            };

#ifdef USE_OPENMP
            #pragma omp critical
#endif
            {
                analysis.ApplyAtomGroupGaussianResult(group_key, class_key, result);
                key_count++;
                Logger::ProgressBar(key_count, group_key_size);
            }
        }
    }
}

} // namespace rhbm_gem::core
