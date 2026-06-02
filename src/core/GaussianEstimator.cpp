#include <rhbm_gem/core/GaussianEstimator.hpp>

#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/hrl/RHBMTrainer.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <atomic>
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

double EstimateResidualIntercept(const LocalPotentialSampleList & residual_sample_entries)
{
    std::vector<double> response_list;
    response_list.reserve(residual_sample_entries.size());
    for (const auto & sample : residual_sample_entries)
    {
        response_list.emplace_back(static_cast<double>(sample.response));
    }
    if (response_list.empty()) return 0.0;
    return array_helper::ComputeMean(response_list.data(), response_list.size());
}

LocalPotentialSampleList BuildShiftedSampleEntries(
    const LocalPotentialSampleList & sample_entries,
    double intercept)
{
    LocalPotentialSampleList shifted_sample_entries;
    shifted_sample_entries.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        shifted_sample_entries.emplace_back(
            LocalPotentialSample{
                static_cast<float>(static_cast<double>(sample.response) - intercept),
                sample.point
            }
        );
    }
    return shifted_sample_entries;
}

LocalPotentialSampleList BuildResidualSampleEntries(
    const LocalPotentialSampleList & sample_entries,
    const RHBMParameterVector & beta,
    LocalGaussianFitModel fit_model,
    double range_min,
    double range_max)
{
    const auto signal_model{ linearization_service::DecodeParameterVector(beta, fit_model) };
    LocalPotentialSampleList residual_sample_entries;
    residual_sample_entries.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
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
                options.distance_max,
                options.local_fit_model)
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
            BuildShiftedSampleEntries(sample_entries_list.at(i), intercept)
        };
        dataset_list.emplace_back(
            rhbm_helper::BuildMemberDataset(
                shifted_sample_entries,
                options.distance_min,
                options.distance_max,
                options.local_fit_model)
        );
    }
    return dataset_list;
}

LocalGaussianResult DecodeLocalGaussianResult(
    double alpha_r,
    const RHBMBetaEstimateResult & fit_result,
    LocalGaussianFitModel fit_model,
    double intercept = 0.0)
{
    const auto ols_model{
        linearization_service::DecodeParameterVector(fit_result.beta_ols, fit_model)
            .WithIntercept(intercept)
    };
    const auto mdpde_model{
        linearization_service::DecodeParameterVector(fit_result.beta_mdpde, fit_model)
            .WithIntercept(intercept)
    };
    return LocalGaussianResult{
        alpha_r,
        GaussianModel3DWithUncertainty{ ols_model, GaussianModel3DUncertainty{} },
        GaussianModel3DWithUncertainty{ mdpde_model, GaussianModel3DUncertainty{} },
        false,
        0.0,
        fit_model,
        fit_result
    };
}

GroupGaussianResult DecodeGroupGaussianResult(
    double alpha_g,
    const RHBMGroupEstimationResult & result,
    LocalGaussianFitModel fit_model)
{
    return GroupGaussianResult{
        alpha_g,
        linearization_service::DecodeParameterVector(result.mu_mean, fit_model),
        linearization_service::DecodeParameterVector(result.mu_mdpde, fit_model),
        linearization_service::DecodeParameterVector(
            result.mu_prior, result.capital_lambda, fit_model)
    };
}

std::vector<LocalGaussianResult> DecodeMemberGaussianResults(
    const RHBMGroupEstimationResult & result,
    LocalGaussianFitModel fit_model)
{
    const auto member_count{ static_cast<std::size_t>(result.beta_posterior_matrix.cols()) };
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
        const auto gaussian{
            linearization_service::DecodeParameterVector(
                result.beta_posterior_matrix.col(i),
                result.capital_sigma_posterior_list.at(static_cast<std::size_t>(i)),
                fit_model)
        };
        member_results.emplace_back(LocalGaussianResult{
            0.0,
            gaussian,
            gaussian,
            static_cast<bool>(result.outlier_flag_array(i)),
            result.statistical_distance_array(i),
            fit_model
        });
    }
    return member_results;
}

std::vector<RHBMBetaEstimateResult> BuildMemberFitResultList(
    const std::vector<LocalGaussianResult> & member_result_list,
    LocalGaussianFitModel fit_model)
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
        if (member_result.fit_model != fit_model)
        {
            throw std::invalid_argument(
                "group Gaussian estimation requires member results to match the requested local fit model.");
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

LocalPotentialSampleList UpdateSampleListWithFittedGaussian(
    const AtomObject & atom,
    const FittedGaussianSnapshot & fitted_gaussian_snapshot)
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
            const auto gaussian_iter{ fitted_gaussian_snapshot.find(neighbor_atom) };
            if (gaussian_iter == fitted_gaussian_snapshot.end()) continue;

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
            if (member_result.fit_model != options.local_fit_model)
            {
                throw std::invalid_argument(
                    "Alpha_G training requires member results to match the requested local fit model.");
            }
            beta_list.emplace_back(
                linearization_service::EncodeGaussianToParameterVector(
                    member_result.mdpde.GetModel(),
                    member_result.fit_model));
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
    const FitOptions & options)
{
    numeric_validation::RequireFiniteNonNegativeRange(
        options.distance_min, options.distance_max, "fit range");
    numeric_validation::RequireFiniteNonNegative(alpha_r, "alpha_r");

    auto execution_options{ MakeExecutionOptions(options) };
    auto dataset{
        rhbm_helper::BuildMemberDataset(
            sample_entries,
            options.distance_min,
            options.distance_max,
            options.local_fit_model)
    };
    const auto result{ rhbm_helper::EstimateBetaMDPDE(alpha_r, dataset, execution_options) };
    return DecodeLocalGaussianResult(alpha_r, result, options.local_fit_model);
}

LocalGaussianResult EstimateLocalGaussianWithIntercept(
    const LocalPotentialSampleList & sample_entries,
    double alpha_r,
    const FitOptions & options)
{
    numeric_validation::RequireFiniteNonNegativeRange(
        options.distance_min, options.distance_max, "fit range");
    numeric_validation::RequireFiniteNonNegative(alpha_r, "alpha_r");

    auto execution_options{ MakeExecutionOptions(options) };
    double intercept{ EstimateInitialIntercept(sample_entries) };
    RHBMBetaEstimateResult result;
    bool has_previous_intercept{ false };

    for (int t = 0; t < execution_options.max_iterations; t++)
    {
        auto shifted_sample_entries{ BuildShiftedSampleEntries(sample_entries, intercept) };
        auto dataset{
            rhbm_helper::BuildMemberDataset(
                shifted_sample_entries,
                options.distance_min,
                options.distance_max,
                options.local_fit_model)
        };
        result = rhbm_helper::EstimateBetaMDPDE(alpha_r, dataset, execution_options);
        const auto residual_sample_entries{
            BuildResidualSampleEntries(
                sample_entries,
                result.beta_mdpde,
                options.local_fit_model,
                options.distance_min,
                options.distance_max)
        };
        const auto next_intercept{ EstimateResidualIntercept(residual_sample_entries) };
        if (has_previous_intercept && std::fabs(next_intercept - intercept) < execution_options.tolerance)
        {
            break;
        }
        has_previous_intercept = true;
        if (t + 1 < execution_options.max_iterations)
        {
            intercept = next_intercept;
        }
    }
    return DecodeLocalGaussianResult(alpha_r, result, options.local_fit_model, intercept);
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
    const auto fit_result_list{ BuildMemberFitResultList(member_result_list, options.local_fit_model) };
    const auto group_input{ rhbm_helper::BuildGroupInput(dataset_list, fit_result_list) };
    const auto raw_result{ rhbm_helper::EstimateGroup(alpha_g, group_input, execution_options) };
    auto result{ DecodeGroupGaussianResult(alpha_g, raw_result, options.local_fit_model) };
    result.member_results = DecodeMemberGaussianResults(raw_result, options.local_fit_model);
    return result;
}

void RunLocalAlphaTraining(ModelObject & model_object, const FitOptions & options)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };
    const auto component_class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };
    const auto component_group_keys{ analysis_view.CollectAtomGroupKeys(component_class_key) };

    size_t count{ 0 };
    for (const auto group_key : component_group_keys)
    {
        const auto & group_atom_list{
            analysis_view.GetAtomObjectList(group_key, component_class_key)
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
        Logger::ProgressPercent(++count, component_group_keys.size());
    }
}

void RunGroupAlphaTraining(ModelObject & model_object, const FitOptions & options)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };
    const auto component_class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };
    const auto component_group_keys{ analysis_view.CollectAtomGroupKeys(component_class_key) };

    RunLocalPotentialFitting(model_object, options);

    std::vector<std::vector<LocalGaussianResult>> member_result_list;
    member_result_list.reserve(component_group_keys.size());
    for (const auto group_key : component_group_keys)
    {
        const auto & group_atom_list{
            analysis_view.GetAtomObjectList(group_key, component_class_key)
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
        const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        for (const auto group_key : analysis_view.CollectAtomGroupKeys(class_key))
        {
            analysis.SetAtomGroupAlphaG(group_key, class_key, alpha_g);
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
        "Run Local atom fitting for " + std::to_string(selected_atom_size) + " atoms.");

#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(options.thread_size)
#endif
    for (size_t i = 0; i < selected_atom_size; i++)
    {
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom_list[i]) };
        auto sample_entries{ local_view.GetSamplingEntries() };
        const auto result{
            EstimateLocalGaussianWithIntercept(sample_entries, local_view.GetAlphaR(), options)
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

    const size_t iter_size{ 10 };
    std::vector<LocalPotentialSampleList> updated_sample_entries_list(selected_atom_size);
    Logger::Log(LogLevel::Info, "Run updated local atom fitting with iterations...");
    for (size_t iter = 0; iter < iter_size; iter++)
    {
        const auto fitted_gaussian_snapshot{ BuildFittedGaussianSnapshot(atom_list) };
#ifdef USE_OPENMP
        #pragma omp parallel for num_threads(options.thread_size)
#endif
        for (size_t i = 0; i < selected_atom_size; i++)
        {
            const auto & atom{ *atom_list[i] };
            const auto local_view{ AtomLocalPotentialView::RequireFor(atom) };
            auto updated_sample_entries{
                UpdateSampleListWithFittedGaussian(atom, fitted_gaussian_snapshot)
            };
            const auto result{
                EstimateLocalGaussian(
                    updated_sample_entries,
                    local_view.GetAlphaR(),
                    options)
            };
            updated_sample_entries_list[i] = std::move(updated_sample_entries);
            local_editor_list[i].SetGaussianResult(result);
        }
        Logger::ProgressBar(iter+1, iter_size);
    }

    for (size_t i = 0; i < selected_atom_size; i++)
    {
        local_editor_list[i].SetSamplingEntries(std::move(updated_sample_entries_list[i]));
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
