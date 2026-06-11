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

#include <atomic>
#include <cmath>
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
constexpr std::size_t kInterceptCycleHistorySize{ 64 };
constexpr std::size_t kLocalFittingMaximumIterations{ 100 };
constexpr double kLocalFittingParameterChangeTolerance{ 1.0e-5 };
constexpr double kLocalFittingChangePercentile{ 0.95 };

struct LocalFittingParameterChangeStats
{
    double max_amplitude_change{ 0.0 };
    double max_width_change{ 0.0 };
    double max_intercept_change{ 0.0 };
    double amplitude_change_percentile{ 0.0 };
    double width_change_percentile{ 0.0 };
    double intercept_change_percentile{ 0.0 };
};

struct LocalFittingAtomChangeStats
{
    double max_atom_change{ 0.0 };
    double atom_change_percentile{ 0.0 };
};

struct LocalFittingIterationResult
{
    std::vector<LocalPotentialSampleList> sample_entries_list;
    std::vector<LocalGaussianResult> result_list;
    std::vector<Eigen::VectorXd> estimation_list;
    std::vector<double> atom_change_list;
};

bool IsSquaredChangeBelowTolerance(double change, double tolerance)
{
    return change * change < tolerance;
}

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

LocalPotentialSampleList BuildSamplesForZeroInterceptGaussianFit(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & model)
{
    const auto zero_intercept_model{ model.WithIntercept(0.0) };
    LocalPotentialSampleList zero_intercept_sample_entries;
    zero_intercept_sample_entries.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        const auto distance{ static_cast<double>(sample.point.distance) };
        const auto model_offset{
            model.ResponseAtDistance(distance) -
                zero_intercept_model.ResponseAtDistance(distance)
        };
        zero_intercept_sample_entries.emplace_back(
            LocalPotentialSample{
                static_cast<float>(
                    static_cast<double>(sample.response) - model_offset),
                sample.point
            }
        );
    }
    return zero_intercept_sample_entries;
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
        const auto zero_intercept_sample_entries{
            BuildSamplesForZeroInterceptGaussianFit(
                sample_entries_list.at(i),
                member_result_list.at(i).mdpde.GetModel())
        };
        dataset_list.emplace_back(
            rhbm_helper::BuildMemberDataset(
                zero_intercept_sample_entries,
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

FittedGaussianSnapshot BuildFittedGaussianSnapshot(
    const std::vector<AtomObject *> & atom_list,
    const std::vector<Eigen::VectorXd> & estimation_list)
{
    if (atom_list.size() != estimation_list.size())
    {
        throw std::invalid_argument("atom_list and estimation_list sizes are inconsistent.");
    }

    FittedGaussianSnapshot snapshot;
    snapshot.reserve(atom_list.size());
    for (std::size_t i = 0; i < atom_list.size(); i++)
    {
        snapshot.emplace(atom_list.at(i), GaussianModel3D::FromVector(estimation_list.at(i)));
    }
    return snapshot;
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
            //if (distance > 2.0) continue; // TEST : Skip long distance contribution
            if (neighbor_atom->GetElement() == Element::OXYGEN)
            {
                response_value -= static_cast<float>(gaussian_iter->second.ResponseAtDistance(distance));
            }
            else
            {
                response_value -= static_cast<float>(gaussian_iter->second.SignalAtDistance(distance));
            }
        }
        updated_list.emplace_back(LocalPotentialSample{response_value, sample.point });
    }
    return updated_list;
}

double CalculateLocalFittingAtomChange(
    const Eigen::VectorXd & current_estimation,
    const Eigen::VectorXd & previous_estimation)
{
    const auto parameter_delta{ current_estimation - previous_estimation };
    const auto amplitude_change{
        std::abs(parameter_delta(GaussianModel3D::AmplitudeIndex()))
    };
    const auto width_change{
        std::abs(parameter_delta(GaussianModel3D::WidthIndex()))
    };
    const auto shape_change{
        amplitude_change > width_change ? amplitude_change : width_change
    };
    const auto intercept_change{
        std::abs(parameter_delta(GaussianModel3D::InterceptIndex()))
    };
    return shape_change > intercept_change ? shape_change : intercept_change;
}

LocalFittingAtomChangeStats CalculateLocalFittingAtomChangeStats(
    const std::vector<double> & atom_change_list)
{
    LocalFittingAtomChangeStats stats;
    for (const auto atom_change : atom_change_list)
    {
        if (atom_change > stats.max_atom_change)
        {
            stats.max_atom_change = atom_change;
        }
    }
    stats.atom_change_percentile = array_helper::ComputePercentile(
        atom_change_list,
        kLocalFittingChangePercentile);
    return stats;
}

std::size_t CountFrozenLocalFittingAtoms(const std::vector<bool> & frozen_atom_list)
{
    std::size_t frozen_count{ 0 };
    for (const auto frozen : frozen_atom_list)
    {
        if (frozen) frozen_count++;
    }
    return frozen_count;
}

LocalFittingParameterChangeStats CalculateLocalFittingParameterChangeStats(
    const std::vector<Eigen::VectorXd> & current_estimation_list,
    const std::vector<Eigen::VectorXd> & previous_estimation_list)
{
    LocalFittingParameterChangeStats stats;
    std::vector<double> amplitude_change_list(current_estimation_list.size());
    std::vector<double> width_change_list(current_estimation_list.size());
    std::vector<double> intercept_change_list(current_estimation_list.size());
    for (size_t i = 0; i < current_estimation_list.size(); i++)
    {
        const auto parameter_delta{ current_estimation_list[i] - previous_estimation_list[i] };
        const auto amplitude_change{
            std::abs(parameter_delta(GaussianModel3D::AmplitudeIndex()))
        };
        const auto width_change{
            std::abs(parameter_delta(GaussianModel3D::WidthIndex()))
        };
        const auto intercept_change{
            std::abs(parameter_delta(GaussianModel3D::InterceptIndex()))
        };
        amplitude_change_list[i] = amplitude_change;
        width_change_list[i] = width_change;
        intercept_change_list[i] = intercept_change;
        if (amplitude_change > stats.max_amplitude_change)
        {
            stats.max_amplitude_change = amplitude_change;
        }
        if (width_change > stats.max_width_change)
        {
            stats.max_width_change = width_change;
        }
        if (intercept_change > stats.max_intercept_change)
        {
            stats.max_intercept_change = intercept_change;
        }
    }

    stats.amplitude_change_percentile = array_helper::ComputePercentile(
        amplitude_change_list,
        kLocalFittingChangePercentile);
    stats.width_change_percentile = array_helper::ComputePercentile(
        width_change_list,
        kLocalFittingChangePercentile);
    stats.intercept_change_percentile = array_helper::ComputePercentile(
        intercept_change_list,
        kLocalFittingChangePercentile);
    return stats;
}

bool IsLocalFittingParameterChangeConverged(const LocalFittingParameterChangeStats & stats)
{
    return
        IsSquaredChangeBelowTolerance(
            stats.max_amplitude_change,
            kLocalFittingParameterChangeTolerance) &&
        IsSquaredChangeBelowTolerance(
            stats.max_width_change,
            kLocalFittingParameterChangeTolerance) &&
        IsSquaredChangeBelowTolerance(
            stats.max_intercept_change,
            kLocalFittingParameterChangeTolerance);
}

double GetLocalFittingParameterChange(const LocalFittingParameterChangeStats & stats)
{
    const auto shape_change{
        stats.amplitude_change_percentile > stats.width_change_percentile ?
            stats.amplitude_change_percentile :
            stats.width_change_percentile
    };
    return shape_change > stats.intercept_change_percentile ?
        shape_change :
        stats.intercept_change_percentile;
}

bool IsBetterLocalFittingCandidate(
    const LocalFittingParameterChangeStats & stats,
    const LocalFittingParameterChangeStats & best_stats)
{
    return GetLocalFittingParameterChange(stats) < GetLocalFittingParameterChange(best_stats);
}

LocalFittingIterationResult RunLocalFittingIteration(
    const std::vector<AtomObject *> & atom_list,
    const std::vector<Eigen::VectorXd> & input_estimation_list,
    const std::vector<LocalGaussianResult> & input_result_list,
    const std::vector<bool> & frozen_atom_list,
    const FitOptions & options)
{
    const auto selected_atom_size{ atom_list.size() };
    if (input_result_list.size() != selected_atom_size ||
        frozen_atom_list.size() != selected_atom_size)
    {
        throw std::invalid_argument("Local fitting iteration input sizes are inconsistent.");
    }
    const auto snapshot{ BuildFittedGaussianSnapshot(atom_list, input_estimation_list) };
    LocalFittingIterationResult iteration_result{
        std::vector<LocalPotentialSampleList>(selected_atom_size),
        std::vector<LocalGaussianResult>(selected_atom_size),
        std::vector<Eigen::VectorXd>(selected_atom_size),
        std::vector<double>(selected_atom_size)
    };

#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(options.thread_size)
#endif
    for (size_t i = 0; i < selected_atom_size; i++)
    {
        const auto & atom{ *atom_list[i] };
        const auto local_view{ AtomLocalPotentialView::RequireFor(atom) };
        auto sample_entries{
            UpdateSampleListWithFittedGaussian(atom, snapshot)
        };
        auto result{ input_result_list.at(i) };
        if (!frozen_atom_list.at(i))
        {
            const auto intercept{ snapshot.at(&atom).GetIntercept() };
            result = EstimateLocalGaussianWithIntercept(
                sample_entries,
                local_view.GetAlphaR(),
                options,
                intercept);
        }
        const auto fitted_model{ result.mdpde.GetModel() };
        iteration_result.sample_entries_list[i] = std::move(sample_entries);
        iteration_result.estimation_list[i] = fitted_model.ToVector();
        iteration_result.atom_change_list[i] = CalculateLocalFittingAtomChange(
            iteration_result.estimation_list[i],
            input_estimation_list.at(i));
        iteration_result.result_list[i] = std::move(result);
    }
    return iteration_result;
}

void RefreshLocalFittingIterationSampleEntries(
    LocalFittingIterationResult & iteration_result,
    const std::vector<AtomObject *> & atom_list)
{
    const auto snapshot{ BuildFittedGaussianSnapshot(atom_list, iteration_result.estimation_list) };
    for (std::size_t i = 0; i < atom_list.size(); i++)
    {
        iteration_result.sample_entries_list.at(i) =
            UpdateSampleListWithFittedGaussian(*atom_list.at(i), snapshot);
    }
}

void ApplyLocalFittingIterationResult(
    const LocalFittingIterationResult & iteration_result,
    std::vector<AtomLocalPotentialEditor> & local_editor_list,
    std::vector<LocalPotentialSampleList> & sample_entries_list)
{
    if (local_editor_list.size() != iteration_result.result_list.size())
    {
        throw std::invalid_argument(
            "local_editor_list and iteration_result sizes are inconsistent.");
    }

    sample_entries_list = iteration_result.sample_entries_list;
    for (std::size_t i = 0; i < local_editor_list.size(); i++)
    {
        local_editor_list.at(i).SetGaussianResult(iteration_result.result_list.at(i));
    }
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
    const auto zero_intercept_sample_entries{
        BuildSamplesForZeroInterceptGaussianFit(
            sample_entries,
            GaussianModel3D{ 0.0, 1.0, intercept })
    };
    auto dataset{
        rhbm_helper::BuildMemberDataset(zero_intercept_sample_entries, range_min, range_max)
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
        const auto raw_intercept{
            EstimateResidualIntercept(sample_entries, *result.fit_result)
        };
        const auto defect{ std::abs(raw_intercept - intercept) };
        if (defect < best_defect)
        {
            best_intercept = intercept;
            best_defect = defect;
        }
        if (defect < tolerance)
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
            auto cycle_result{
                EstimateLocalGaussian(sample_entries, alpha_r, options, final_intercept)
            };
            const auto cycle_raw_intercept{
                EstimateResidualIntercept(sample_entries, *cycle_result.fit_result)
            };
            const auto cycle_defect{ std::abs(cycle_raw_intercept - final_intercept) };
            if (cycle_defect < best_defect)
            {
                best_intercept = final_intercept;
                best_defect = cycle_defect;
            }
            if (cycle_defect < tolerance)
            {
                result = std::move(cycle_result);
                if (!options.quiet_mode)
                {
                    Logger::Log(LogLevel::Debug,
                        "Cycle detected in local Gaussian intercept estimation with period " +
                        std::to_string(period) + "; refitting at cycle mean with defect = " +
                        std::to_string(cycle_defect) + ".");
                }
            }
            else
            {
                result = EstimateLocalGaussian(sample_entries, alpha_r, options, best_intercept);
                if (!options.quiet_mode)
                {
                    Logger::Log(LogLevel::Debug,
                        "Cycle mean fixed-point defect = " + std::to_string(cycle_defect) +
                        " did not satisfy local Gaussian intercept tolerance; "
                        "refitting at best fixed-point candidate with defect = " +
                        std::to_string(best_defect) + ".");
                }
            }
            cycle_detected = true;
            break;
        }
        if (cycle_detected) break;

        if (t + 1 == max_iterations)
        {
            result = EstimateLocalGaussian(sample_entries, alpha_r, options, best_intercept);
            if (!options.quiet_mode)
            {
                Logger::Log(LogLevel::Warning,
                    "Maximum iterations reached in local Gaussian estimation with intercept; "
                    "refitting at best fixed-point candidate with defect = " +
                    std::to_string(best_defect) + ".");
            }
            break;
        }

        if (intercept_history.size() == kInterceptCycleHistorySize)
        {
            intercept_history.erase(intercept_history.begin());
        }
        intercept_history.emplace_back(intercept);
        intercept = raw_intercept;
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
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run local alpha training for " + std::to_string(group_key_list.size()) + " groups.");
    }
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
        count++;
        if (!options.quiet_mode)
        {
            Logger::ProgressPercent(count, group_key_list.size());
        }
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

void RunFirstStageLocalFitting(ModelObject & model_object, const FitOptions & options)
{
    const auto selected_atom_size{ model_object.GetSelectedAtomCount() };
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    auto local_editor_list{ BuildSelectedAtomLocalEditors(model_object) };
    std::atomic<size_t> atom_count{ 0 };

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
            if (!options.quiet_mode)
            {
                Logger::ProgressPercent(atom_count, selected_atom_size);
            }
        }
    }
}

void RunSecondStageLocalFitting(
    ModelObject & model_object,
    const std::vector<AtomObject *> & atom_list,
    const FitOptions & options)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto atom_size{ atom_list.size() };
    std::vector<AtomLocalPotentialEditor> local_editor_list;
    local_editor_list.reserve(atom_list.size());
    for (auto * atom : atom_list)
    {
        local_editor_list.emplace_back(analysis.EnsureAtomLocalPotential(*atom));
    }

    std::vector<LocalPotentialSampleList> sample_entries_list(atom_size);
    std::vector<Eigen::VectorXd> previous_estimation_list(atom_size);
    std::vector<LocalGaussianResult> previous_result_list(atom_size);
    for (size_t i = 0; i < atom_size; i++)
    {
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom_list[i]) };
        previous_result_list[i] = local_view.GetGaussianResult();
        previous_estimation_list[i] = previous_result_list[i].mdpde.GetModel().ToVector();
    }

    LocalFittingIterationResult best_iteration_result;
    LocalFittingParameterChangeStats best_change_stats;
    bool has_best_iteration_result{ false };
    std::vector<bool> frozen_atom_list(atom_size, false);
    for (size_t iter = 0; iter < kLocalFittingMaximumIterations; iter++)
    {
        auto iteration_result{
            RunLocalFittingIteration(
                atom_list,
                previous_estimation_list,
                previous_result_list,
                frozen_atom_list,
                options)
        };
        const auto change_stats{
            CalculateLocalFittingParameterChangeStats(
                iteration_result.estimation_list,
                previous_estimation_list)
        };
        const auto atom_change_stats{
            CalculateLocalFittingAtomChangeStats(iteration_result.atom_change_list)
        };
        for (std::size_t i = 0; i < atom_size; i++)
        {
            if (!frozen_atom_list.at(i) &&
                iteration_result.atom_change_list.at(i) < kLocalFittingParameterChangeTolerance)
            {
                frozen_atom_list.at(i) = true;
            }
        }
        const auto frozen_atom_count{ CountFrozenLocalFittingAtoms(frozen_atom_list) };
        if (!has_best_iteration_result ||
            IsBetterLocalFittingCandidate(change_stats, best_change_stats))
        {
            best_iteration_result = iteration_result;
            best_change_stats = change_stats;
            has_best_iteration_result = true;
        }

        if (!options.quiet_mode)
        {
            std::ostringstream progress_message;
            progress_message << "Local fitting iteration " << iter + 1 << '/'
                << kLocalFittingMaximumIterations
                << std::fixed << std::setprecision(5)
                //<< ", max amplitude change = "
                //<< change_stats.max_amplitude_change
                << ", percentile amplitude change = "
                << change_stats.amplitude_change_percentile
                //<< ", max width change = "
                //<< change_stats.max_width_change
                << ", percentile width change = "
                << change_stats.width_change_percentile
                //<< ", max intercept change = "
                //<< change_stats.max_intercept_change
                << ", percentile intercept change = "
                << change_stats.intercept_change_percentile
                << ", percentile atom change = "
                << atom_change_stats.atom_change_percentile
                << ", frozen atoms = "
                << frozen_atom_count << '/' << atom_size;
            Logger::ProgressLine(progress_message.str());
        }

        if (frozen_atom_count == atom_size)
        {
            RefreshLocalFittingIterationSampleEntries(iteration_result, atom_list);
            ApplyLocalFittingIterationResult(
                iteration_result,
                local_editor_list,
                sample_entries_list);
            if (!options.quiet_mode)
            {
                Logger::FinishProgressLine();
                Logger::Log(LogLevel::Info,
                    "All local fitting atoms frozen after " + std::to_string(iter + 1) +
                    " iterations with max atom change = " +
                    std::to_string(atom_change_stats.max_atom_change) +
                    " and percentile atom change = " +
                    std::to_string(atom_change_stats.atom_change_percentile) + ".");
            }
            break;
        }

        const auto converged{ IsLocalFittingParameterChangeConverged(change_stats) };
        if (converged)
        {
            RefreshLocalFittingIterationSampleEntries(iteration_result, atom_list);
            ApplyLocalFittingIterationResult(
                iteration_result,
                local_editor_list,
                sample_entries_list);
            if (!options.quiet_mode)
            {
                Logger::FinishProgressLine();
                Logger::Log(LogLevel::Info,
                    "Converged after " + std::to_string(iter + 1) +
                    " iterations with max amplitude change = " +
                    std::to_string(change_stats.max_amplitude_change) +
                    ", percentile amplitude change = " +
                    std::to_string(change_stats.amplitude_change_percentile) +
                    ", max width change = " +
                    std::to_string(change_stats.max_width_change) +
                    ", percentile width change = " +
                    std::to_string(change_stats.width_change_percentile) +
                    ", max intercept change = " +
                    std::to_string(change_stats.max_intercept_change) +
                    ", and percentile intercept change = " +
                    std::to_string(change_stats.intercept_change_percentile) + ".");
            }
            break;
        }

        if (iter + 1 == kLocalFittingMaximumIterations)
        {
            RefreshLocalFittingIterationSampleEntries(best_iteration_result, atom_list);
            ApplyLocalFittingIterationResult(
                best_iteration_result,
                local_editor_list,
                sample_entries_list);
            if (!options.quiet_mode)
            {
                Logger::FinishProgressLine();
                Logger::Log(LogLevel::Warning,
                    "Reached maximum iteration size; refitting at best fixed-point candidate "
                    "with max amplitude change = " +
                    std::to_string(best_change_stats.max_amplitude_change) +
                    ", percentile amplitude change = " +
                    std::to_string(best_change_stats.amplitude_change_percentile) +
                    ", max width change = " +
                    std::to_string(best_change_stats.max_width_change) +
                    ", percentile width change = " +
                    std::to_string(best_change_stats.width_change_percentile) +
                    ", max intercept change = " +
                    std::to_string(best_change_stats.max_intercept_change) +
                    ", and percentile intercept change = " +
                    std::to_string(best_change_stats.intercept_change_percentile));
            }
        }
        previous_estimation_list = std::move(iteration_result.estimation_list);
        previous_result_list = std::move(iteration_result.result_list);
    }

    for (size_t i = 0; i < atom_size; i++)
    {
        local_editor_list[i].SetSamplingEntries(std::move(sample_entries_list[i]));
    }
}

void RunLocalPotentialFitting(ModelObject & model_object, const FitOptions & options)
{
    RunFirstStageLocalFitting(model_object, options);

    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run updated local atom fitting with iterations...");
    }

    const auto & atom_list{ model_object.GetSelectedAtoms() };
    RunSecondStageLocalFitting(model_object, atom_list, options);
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
        if (!options.quiet_mode)
        {
            Logger::Log(LogLevel::Info, "Class type: " + class_key);
        }

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
                if (!options.quiet_mode)
                {
                    Logger::ProgressBar(key_count, group_key_size);
                }
            }
        }
    }
}

} // namespace rhbm_gem::core
