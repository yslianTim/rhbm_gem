#include <cstddef>
#include <rhbm_gem/core/GaussianEstimator.hpp>

#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>
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
#include <cmath>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace rhbm_gem::core {
namespace {
constexpr std::size_t kMinimumAlphaRTrainingSampleCount{ 10 };
constexpr std::size_t kMinimumAlphaGTrainingMemberCount{ 10 };
constexpr double kResidualOffsetRangeMin{ 1.0 };
constexpr double kResidualOffsetRangeMax{ 2.0 };
constexpr double kOffsetDampingFactor{ 0.5 };
constexpr double kNeighborContributionDistanceMax{ 2.5 };
constexpr std::size_t kLocalFittingMaximumIterations{ 100 };
constexpr double kLocalFittingParameterChangeTolerance{ 1.0e-6 };
constexpr double kLocalFittingChangePercentile{ 0.95 };
constexpr int kHuberSlopeMaximumIterations{ 50 };
constexpr double kHuberSlopeTolerance{ 1.0e-8 };
constexpr double kHuberScaleMultiplier{ 1.4826 };
constexpr double kHuberScaleMin{ 1.0e-12 };
constexpr double kHuberCutoffMultiplier{ 1.345 };
constexpr double kOffsetRegularizationAmplitudeRatio{ 0.1 };
constexpr double kOffsetRegularizationPriorScaleMin{ 1.0e-12 };
constexpr double kJointOffsetRidgeRatio{ 1.0e-3 };

struct ResidualOffsetSample
{
    double basis{ 0.0 };
    double residual{ 0.0 };
};

struct LocalFittingParameterChangeStats
{
    double amplitude_change_percentile{ 0.0 };
    double width_change_percentile{ 0.0 };
    double offset_change_percentile{ 0.0 };
};

struct JointOffsetSystem
{
    Eigen::SparseMatrix<double> design_matrix;
    Eigen::VectorXd response;
    Eigen::VectorXd previous_offset;
    Eigen::VectorXd ridge_diagonal;
};

struct JointOffsetRow
{
    std::vector<std::pair<Eigen::Index, double>> basis_entries;
    double response{ 0.0 };
};

struct LocalFittingIterationResult
{
    std::vector<LocalGaussianResult> result_list;
    std::vector<Eigen::VectorXd> estimation_list;
};

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

bool CanBuildFiniteZeroOffsetSamples(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & model)
{
    for (const auto & sample : sample_entries)
    {
        const auto distance{ static_cast<double>(sample.point.distance) };
        const auto model_offset{ model.ResponseAtDistance(distance) - model.SignalAtDistance(distance) };
        const auto response{ static_cast<double>(sample.response) - model_offset };
        if (!std::isfinite(response)) return false;
        if (std::abs(response) > static_cast<double>(std::numeric_limits<float>::max()))
        {
            return false;
        }
    }
    return true;
}

bool EstimateOrdinarySlopeThroughOrigin(
    const std::vector<ResidualOffsetSample> & sample_list,
    double & slope)
{
    double numerator{ 0.0 };
    double denominator{ 0.0 };
    for (const auto & sample : sample_list)
    {
        numerator += sample.basis * sample.residual;
        denominator += sample.basis * sample.basis;
    }
    if (!std::isfinite(numerator) ||
        !std::isfinite(denominator) ||
        denominator <= std::numeric_limits<double>::epsilon())
    {
        return false;
    }
    slope = numerator / denominator;
    return std::isfinite(slope);
}

double ComputeHuberResidualScale(
    const std::vector<ResidualOffsetSample> & sample_list,
    double slope)
{
    std::vector<double> residual_list;
    residual_list.reserve(sample_list.size());
    for (const auto & sample : sample_list)
    {
        residual_list.emplace_back(sample.residual - slope * sample.basis);
    }

    const auto median_residual{ array_helper::ComputeMedian(residual_list) };
    std::vector<double> deviation_list;
    deviation_list.reserve(residual_list.size());
    for (const auto residual : residual_list)
    {
        deviation_list.emplace_back(std::abs(residual - median_residual));
    }

    return std::max(
        kHuberScaleMultiplier * array_helper::ComputeMedian(deviation_list),
        kHuberScaleMin);
}

double ComputeOffsetRegularizationLambda(double residual_scale, double amplitude)
{
    const auto prior_scale{
        std::max(
            std::abs(amplitude) * kOffsetRegularizationAmplitudeRatio,
            kOffsetRegularizationPriorScaleMin)
    };
    const auto lambda{ residual_scale / prior_scale };
    return lambda * lambda;
}

bool EstimateHuberSlopeThroughOrigin(
    const std::vector<ResidualOffsetSample> & sample_list,
    double amplitude,
    double & slope)
{
    if (sample_list.empty())
    {
        return false;
    }
    if (!EstimateOrdinarySlopeThroughOrigin(sample_list, slope))
    {
        return false;
    }

    for (int t = 0; t < kHuberSlopeMaximumIterations; t++)
    {
        const auto scale{ ComputeHuberResidualScale(sample_list, slope) };
        const auto cutoff{ kHuberCutoffMultiplier * scale };
        const auto lambda{ ComputeOffsetRegularizationLambda(scale, amplitude) };
        double numerator{ 0.0 };
        double denominator{ 0.0 };
        for (const auto & sample : sample_list)
        {
            const auto error{ sample.residual - slope * sample.basis };
            const auto abs_error{ std::abs(error) };
            const auto weight{ abs_error <= cutoff ? 1.0 : cutoff / abs_error };
            numerator += weight * sample.basis * sample.residual;
            denominator += weight * sample.basis * sample.basis;
        }
        denominator += lambda;
        if (!std::isfinite(numerator) ||
            !std::isfinite(denominator) ||
            denominator <= std::numeric_limits<double>::epsilon())
        {
            return false;
        }

        const auto updated_slope{ numerator / denominator };
        if (!std::isfinite(updated_slope))
        {
            return false;
        }
        if (std::abs(updated_slope - slope) < kHuberSlopeTolerance)
        {
            slope = updated_slope;
            return true;
        }
        slope = updated_slope;
    }
    return true;
}

double EstimateResidualOffsetParameter(
    const LocalPotentialSampleList & sample_entries,
    const RHBMBetaEstimateResult & fit_result,
    double current_offset)
{
    const auto signal_model{ linearization_service::DecodeParameterVector(fit_result.beta_mdpde) };
    const auto width{ signal_model.GetWidth() };
    if (!std::isfinite(width) || width <= 0.0) return current_offset;

    std::vector<ResidualOffsetSample> residual_sample_list;
    residual_sample_list.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        const auto distance{ static_cast<double>(sample.point.distance) };
        if (distance < kResidualOffsetRangeMin) continue;
        if (distance > kResidualOffsetRangeMax) continue;

        const auto basis{ signal_model.OffsetBasisAtDistance(distance) };
        if (!std::isfinite(basis) || std::abs(basis) <= std::numeric_limits<double>::epsilon())
        {
            continue;
        }
        const auto residual{
            static_cast<double>(sample.response) - signal_model.SignalAtDistance(distance)
        };
        if (!std::isfinite(residual))
        {
            continue;
        }
        residual_sample_list.emplace_back(ResidualOffsetSample{ basis, residual });
    }
    double candidate_offset{ current_offset };
    if (!EstimateHuberSlopeThroughOrigin(
            residual_sample_list,
            signal_model.GetAmplitude(),
            candidate_offset))
    {
        return current_offset;
    }
    const auto candidate_model{ signal_model.WithOffset(candidate_offset) };
    if (!CanBuildFiniteZeroOffsetSamples(sample_entries, candidate_model))
    {
        return current_offset;
    }
    return candidate_offset;
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
    LocalPotentialSampleList updated_sample_entries;
    updated_sample_entries.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        const auto distance{ static_cast<double>(sample.point.distance) };
        const auto model_offset{ model.ResponseAtDistance(distance) - model.SignalAtDistance(distance) };
        const auto response{ static_cast<float>(static_cast<double>(sample.response) - model_offset)};
        updated_sample_entries.emplace_back(LocalPotentialSample{ response, sample.point });
    }
    return updated_sample_entries;
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
    auto range_min{ options.distance_min };
    auto range_max{ options.distance_max };
    std::vector<RHBMMemberDataset> dataset_list;
    dataset_list.reserve(sample_entries_list.size());
    for (std::size_t i = 0; i < sample_entries_list.size(); i++)
    {
        const auto sampling_entries{
            BuildSamplesForZeroOffsetGaussianFit(
                sample_entries_list.at(i),
                member_result_list.at(i).mdpde.GetModel())
        };
        dataset_list.emplace_back(rhbm_helper::BuildMemberDataset(sampling_entries, range_min, range_max));
    }
    return dataset_list;
}

LocalGaussianResult DecodeLocalGaussianResult(
    double alpha_r,
    const RHBMBetaEstimateResult & fit_result,
    double offset = 0.0)
{
    const auto ols_model{
        linearization_service::DecodeParameterVector(fit_result.beta_ols)
            .WithOffset(offset)
    };
    const auto mdpde_model{
        linearization_service::DecodeParameterVector(fit_result.beta_mdpde)
            .WithOffset(offset)
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

GaussianModel3DWithUncertainty WithModelOffset(
    const GaussianModel3DWithUncertainty & gaussian,
    double offset)
{
    return GaussianModel3DWithUncertainty{
        gaussian.GetModel().WithOffset(offset),
        gaussian.GetStandardDeviationModel()
    };
}

GroupGaussianResult DecodeGroupGaussianResult(
    double alpha_g,
    const RHBMGroupEstimationResult & result,
    double offset)
{
    return GroupGaussianResult{
        alpha_g,
        linearization_service::DecodeParameterVector(result.mu_mean).WithOffset(offset),
        linearization_service::DecodeParameterVector(result.mu_mdpde).WithOffset(offset),
        WithModelOffset(
            linearization_service::DecodeParameterVector(result.mu_prior, result.capital_lambda),
            offset)
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
            WithModelOffset(
                linearization_service::DecodeParameterVector(
                    result.beta_posterior_matrix.col(i),
                    result.capital_sigma_posterior_list.at(member_index)),
                offset)
        };
        member_results.emplace_back(LocalGaussianResult{
            0.0,
            gaussian,
            gaussian,
            gaussian,
            static_cast<bool>(result.outlier_flag_array(i)),
            result.statistical_distance_array(i)
        });
    }
    return member_results;
}

std::vector<RHBMBetaEstimateResult> BuildMemberFitResultList(
    const std::vector<RHBMMemberDataset> & dataset_list,
    const std::vector<LocalGaussianResult> & member_result_list,
    const FitOptions & options)
{
    if (dataset_list.size() != member_result_list.size())
    {
        throw std::invalid_argument("dataset_list and member_result_list sizes are inconsistent.");
    }
    const auto execution_options{ MakeExecutionOptions(options) };
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

JointOffsetSystem BuildJointOffsetSystem(
    const std::vector<AtomObject *> & atom_list,
    const FittedGaussianSnapshot & snapshot)
{
    std::unordered_map<const AtomObject *, Eigen::Index> atom_column_map;
    atom_column_map.reserve(atom_list.size());
    for (std::size_t i = 0; i < atom_list.size(); i++)
    {
        atom_column_map.emplace(atom_list.at(i), static_cast<Eigen::Index>(i));
    }

    std::vector<JointOffsetRow> row_list;
    for (const auto * atom : atom_list)
    {
        const auto model_iter{ snapshot.find(atom) };
        if (model_iter == snapshot.end())
        {
            throw std::invalid_argument("Joint offset snapshot is missing an atom.");
        }
        const auto target_column{ atom_column_map.at(atom) };
        const auto & target_model{ model_iter->second };
        const auto sample_entries{
            AtomLocalPotentialView::RequireFor(*atom).GetSamplingEntries(false)
        };
        const auto neighbor_atom_list{ atom->FindNeighborAtoms() };
        for (const auto & sample : sample_entries)
        {
            if (!std::isfinite(static_cast<double>(sample.response)))
            {
                throw std::runtime_error("Joint offset sample response is not finite.");
            }
            const auto target_distance{ static_cast<double>(sample.point.distance) };
            const auto target_signal{ target_model.SignalAtDistance(target_distance) };
            const auto target_basis{ target_model.OffsetBasisAtDistance(target_distance) };
            if (!std::isfinite(target_signal) || !std::isfinite(target_basis))
            {
                throw std::runtime_error("Joint offset target model evaluation is not finite.");
            }
            auto residual{ static_cast<double>(sample.response) - target_signal };
            JointOffsetRow row;
            if (std::abs(target_basis) > std::numeric_limits<double>::epsilon())
            {
                row.basis_entries.emplace_back(target_column, target_basis);
            }

            for (const auto * neighbor_atom : neighbor_atom_list)
            {
                const auto column_iter{ atom_column_map.find(neighbor_atom) };
                if (column_iter == atom_column_map.end()) continue;
                const auto neighbor_model_iter{ snapshot.find(neighbor_atom) };
                if (neighbor_model_iter == snapshot.end()) continue;

                const auto distance{
                    static_cast<double>(
                        array_helper::ComputeNorm<float>(sample.point.position, neighbor_atom->GetPosition()))
                };
                if (distance > kNeighborContributionDistanceMax) continue;
                const auto & neighbor_model{ neighbor_model_iter->second };
                const auto signal{ neighbor_model.SignalAtDistance(distance) };
                const auto basis{ neighbor_model.OffsetBasisAtDistance(distance) };
                if (!std::isfinite(signal) || !std::isfinite(basis))
                {
                    throw std::runtime_error(
                        "Joint offset neighbor model evaluation is not finite.");
                }
                residual -= signal;
                if (std::abs(basis) > std::numeric_limits<double>::epsilon())
                {
                    row.basis_entries.emplace_back(column_iter->second, basis);
                }
            }
            if (!std::isfinite(residual))
            {
                throw std::runtime_error("Joint offset residual is not finite.");
            }
            if (row.basis_entries.empty()) continue;
            row.response = residual;
            row_list.emplace_back(std::move(row));
        }
    }

    const auto row_count{ static_cast<Eigen::Index>(row_list.size()) };
    const auto column_count{ static_cast<Eigen::Index>(atom_list.size()) };
    std::vector<Eigen::Triplet<double>> triplet_list;
    Eigen::VectorXd response{ Eigen::VectorXd::Zero(row_count) };
    Eigen::VectorXd column_square_sum{ Eigen::VectorXd::Zero(column_count) };
    for (Eigen::Index row_index = 0; row_index < row_count; row_index++)
    {
        const auto & row{ row_list.at(static_cast<std::size_t>(row_index)) };
        response(row_index) = row.response;
        for (const auto & [column_index, basis] : row.basis_entries)
        {
            triplet_list.emplace_back(row_index, column_index, basis);
            column_square_sum(column_index) += basis * basis;
        }
    }

    JointOffsetSystem system;
    system.design_matrix.resize(row_count, column_count);
    system.design_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
    system.response = std::move(response);
    system.previous_offset = Eigen::VectorXd::Zero(column_count);
    system.ridge_diagonal = Eigen::VectorXd::Zero(column_count);
    for (Eigen::Index column_index = 0; column_index < column_count; column_index++)
    {
        const auto & model{ snapshot.at(atom_list.at(static_cast<std::size_t>(column_index))) };
        system.previous_offset(column_index) = model.GetOffset();
        const auto square_sum{ column_square_sum(column_index) };
        system.ridge_diagonal(column_index) =
            square_sum > std::numeric_limits<double>::epsilon() ? kJointOffsetRidgeRatio * square_sum : 1.0;
    }
    return system;
}

bool SolveJointOffsetWeightedRidge(
    const JointOffsetSystem & system,
    const Eigen::VectorXd & weight,
    Eigen::VectorXd & offset)
{
    auto weighted_design{ system.design_matrix };
    for (Eigen::Index column = 0; column < weighted_design.outerSize(); column++)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator iter(weighted_design, column);
             iter;
             ++iter)
        {
            iter.valueRef() *= std::sqrt(weight(iter.row()));
        }
    }
    const auto weighted_response{ weight.array().sqrt().matrix().cwiseProduct(system.response) };
    Eigen::SparseMatrix<double> normal_matrix{
        weighted_design.transpose() * weighted_design
    };
    const auto column_count{ normal_matrix.cols() };
    for (Eigen::Index i = 0; i < column_count; i++)
    {
        normal_matrix.coeffRef(i, i) += system.ridge_diagonal(i);
    }
    normal_matrix.makeCompressed();
    const Eigen::VectorXd right_hand_side{
        weighted_design.transpose() * weighted_response +
        system.ridge_diagonal.cwiseProduct(system.previous_offset)
    };

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(normal_matrix);
    if (solver.info() != Eigen::Success) return false;
    offset = solver.solve(right_hand_side);
    return solver.info() == Eigen::Success && offset.allFinite();
}

Eigen::VectorXd EstimateJointOffsets(
    const std::vector<AtomObject *> & atom_list,
    const FittedGaussianSnapshot & snapshot)
{
    Eigen::VectorXd previous_offset{
        Eigen::VectorXd::Zero(static_cast<Eigen::Index>(atom_list.size()))
    };
    for (std::size_t i = 0; i < atom_list.size(); i++)
    {
        previous_offset(static_cast<Eigen::Index>(i)) = snapshot.at(atom_list.at(i)).GetOffset();
    }
    JointOffsetSystem system;
    try
    {
        system = BuildJointOffsetSystem(atom_list, snapshot);
    }
    catch (const std::exception &)
    {
        return previous_offset;
    }
    if (system.response.size() == 0 || system.previous_offset.size() == 0)
    {
        return previous_offset;
    }

    Eigen::VectorXd weight{ Eigen::VectorXd::Ones(system.response.size()) };
    Eigen::VectorXd offset;
    if (!SolveJointOffsetWeightedRidge(system, weight, offset))
    {
        return previous_offset;
    }

    for (int iteration = 0; iteration < kHuberSlopeMaximumIterations; iteration++)
    {
        const Eigen::VectorXd residual{ system.response - system.design_matrix * offset };
        std::vector<double> residual_list(residual.data(), residual.data() + residual.size());
        const auto median_residual{ array_helper::ComputeMedian(residual_list) };
        std::vector<double> deviation_list;
        deviation_list.reserve(residual_list.size());
        for (const auto value : residual_list)
        {
            deviation_list.emplace_back(std::abs(value - median_residual));
        }
        const auto residual_scale{ std::max(
            kHuberScaleMultiplier * array_helper::ComputeMedian(deviation_list),
            kHuberScaleMin)
        };
        const auto cutoff{ kHuberCutoffMultiplier * residual_scale };
        for (Eigen::Index i = 0; i < residual.size(); i++)
        {
            const auto absolute_residual{ std::abs(residual(i)) };
            weight(i) = absolute_residual <= cutoff ? 1.0 : cutoff / absolute_residual;
        }

        Eigen::VectorXd updated_offset;
        if (!SolveJointOffsetWeightedRidge(system, weight, updated_offset))
        {
            return system.previous_offset;
        }
        const auto maximum_change{ (updated_offset - offset).cwiseAbs().maxCoeff() };
        offset = std::move(updated_offset);
        if (maximum_change < kHuberSlopeTolerance) break;
    }

    return offset;
}

FittedGaussianSnapshot WithJointOffsets(
    const std::vector<AtomObject *> & atom_list,
    const FittedGaussianSnapshot & snapshot,
    const Eigen::VectorXd & offset)
{
    if (atom_list.size() != static_cast<std::size_t>(offset.size()))
    {
        throw std::invalid_argument("Joint offset result size is inconsistent.");
    }
    auto updated_snapshot{ snapshot };
    for (std::size_t i = 0; i < atom_list.size(); i++)
    {
        updated_snapshot.at(atom_list.at(i)) =
            updated_snapshot.at(atom_list.at(i)).WithOffset(
                offset(static_cast<Eigen::Index>(i)));
    }
    return updated_snapshot;
}

template <typename GaussianLookup>
LocalPotentialSampleList UpdateSampleListWithGaussianLookup(
    const AtomObject & atom,
    GaussianLookup lookup_gaussian)
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
            const auto * gaussian{ lookup_gaussian(*neighbor_atom) };
            if (gaussian == nullptr) continue;

            auto neighbor_position{ neighbor_atom->GetPosition() };
            auto distance{
                static_cast<double>(
                    array_helper::ComputeNorm<float>(sample_position, neighbor_position))
            };
            if (distance > kNeighborContributionDistanceMax) continue;
            response_value -= static_cast<float>(gaussian->ResponseAtDistance(distance));
        }
        updated_list.emplace_back(LocalPotentialSample{response_value, sample.point });
    }
    return updated_list;
}

LocalPotentialSampleList UpdateSampleListWithFittedGaussian(
    const AtomObject & atom,
    const FittedGaussianSnapshot & snapshot)
{
    return UpdateSampleListWithGaussianLookup(
        atom,
        [&snapshot](const AtomObject & neighbor_atom) -> const GaussianModel3D *
        {
            const auto gaussian_iter{ snapshot.find(&neighbor_atom) };
            return gaussian_iter == snapshot.end() ? nullptr : &gaussian_iter->second;
        });
}

LocalPotentialSampleList UpdateSampleListWithFittedGaussian(const AtomObject & atom)
{
    return UpdateSampleListWithGaussianLookup(
        atom,
        [](const AtomObject & neighbor_atom) -> const GaussianModel3D *
        {
            const auto local_view{ AtomLocalPotentialView::For(neighbor_atom) };
            return local_view.IsAvailable() ? &local_view.GetEstimateMDPDE() : nullptr;
        });
}

LocalFittingParameterChangeStats CalculateLocalFittingParameterChangeStats(
    const std::vector<Eigen::VectorXd> & current_estimation_list,
    const std::vector<Eigen::VectorXd> & previous_estimation_list)
{
    LocalFittingParameterChangeStats stats;
    std::vector<double> amplitude_change_list(current_estimation_list.size());
    std::vector<double> width_change_list(current_estimation_list.size());
    std::vector<double> offset_change_list(current_estimation_list.size());
    for (size_t i = 0; i < current_estimation_list.size(); i++)
    {
        const auto parameter_delta{ current_estimation_list[i] - previous_estimation_list[i] };
        const auto amplitude_change{
            std::abs(parameter_delta(GaussianModel3D::AmplitudeIndex()))
        };
        const auto width_change{
            std::abs(parameter_delta(GaussianModel3D::WidthIndex()))
        };
        const auto offset_change{
            std::abs(parameter_delta(GaussianModel3D::OffsetIndex()))
        };
        amplitude_change_list[i] = amplitude_change;
        width_change_list[i] = width_change;
        offset_change_list[i] = offset_change;
    }

    stats.amplitude_change_percentile = array_helper::ComputePercentile(
        amplitude_change_list,
        kLocalFittingChangePercentile);
    stats.width_change_percentile = array_helper::ComputePercentile(
        width_change_list,
        kLocalFittingChangePercentile);
    stats.offset_change_percentile = array_helper::ComputePercentile(
        offset_change_list,
        kLocalFittingChangePercentile);
    return stats;
}

bool IsLocalFittingParameterChangeConverged(const LocalFittingParameterChangeStats & stats)
{
    return
        (std::pow(stats.amplitude_change_percentile, 2) < kLocalFittingParameterChangeTolerance) &&
        (std::pow(stats.width_change_percentile, 2) < kLocalFittingParameterChangeTolerance) &&
        (std::pow(stats.offset_change_percentile, 2) < kLocalFittingParameterChangeTolerance);
}

double GetLocalFittingParameterChange(const LocalFittingParameterChangeStats & stats)
{
    const auto shape_change{
        stats.amplitude_change_percentile > stats.width_change_percentile ?
            stats.amplitude_change_percentile :
            stats.width_change_percentile
    };
    return shape_change > stats.offset_change_percentile ?
        shape_change :
        stats.offset_change_percentile;
}

bool IsBetterLocalFittingCandidate(
    const LocalFittingParameterChangeStats & stats,
    const LocalFittingParameterChangeStats & best_stats)
{
    return GetLocalFittingParameterChange(stats) < GetLocalFittingParameterChange(best_stats);
}

void ApplyLocalFittingUnderRelaxation(
    LocalFittingIterationResult & iteration_result,
    const std::vector<Eigen::VectorXd> & previous_estimation_list,
    double beta)
{
    if (iteration_result.estimation_list.size() != previous_estimation_list.size() ||
        iteration_result.result_list.size() != previous_estimation_list.size())
    {
        throw std::invalid_argument("Local fitting relaxation input sizes are inconsistent.");
    }
    for (std::size_t i = 0; i < iteration_result.estimation_list.size(); i++)
    {
        auto relaxed_estimation{
            (beta * iteration_result.estimation_list.at(i) +
            (1.0 - beta) * previous_estimation_list.at(i)).eval()
        };
        const auto relaxed_model{ GaussianModel3D::FromVector(relaxed_estimation) };
        auto & result{ iteration_result.result_list.at(i) };
        result.mdpde = GaussianModel3DWithUncertainty{
            relaxed_model,
            result.mdpde.GetStandardDeviationModel()
        };
        iteration_result.estimation_list.at(i) = relaxed_estimation;
    }
}

LocalGaussianResult EstimateLocalGaussianWithOffsetModel(
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
    const auto updated_sample_entries{
        BuildSamplesForZeroOffsetGaussianFit(sample_entries, offset_model)
    };
    auto dataset{
        rhbm_helper::BuildMemberDataset(updated_sample_entries, range_min, range_max)
    };
    const auto result{ rhbm_helper::EstimateBetaMDPDE(alpha_r, dataset, execution_options) };
    return DecodeLocalGaussianResult(alpha_r, result, offset_model.GetOffset());
}

LocalFittingIterationResult RunLocalFittingIteration(
    const std::vector<AtomObject *> & atom_list,
    const std::vector<Eigen::VectorXd> & input_estimation_list,
    const std::vector<LocalGaussianResult> & input_result_list,
    const FitOptions & options)
{
    const auto selected_atom_size{ atom_list.size() };
    if (input_result_list.size() != selected_atom_size)
    {
        throw std::invalid_argument("Local fitting iteration input sizes are inconsistent.");
    }
    const auto current_snapshot{
        BuildFittedGaussianSnapshot(atom_list, input_estimation_list)
    };
    const auto joint_offset{
        EstimateJointOffsets(atom_list, current_snapshot)
    };
    const auto offset_snapshot{
        WithJointOffsets(atom_list, current_snapshot, joint_offset)
    };
    LocalFittingIterationResult iteration_result{
        std::vector<LocalGaussianResult>(selected_atom_size),
        std::vector<Eigen::VectorXd>(selected_atom_size)
    };

#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(options.thread_size)
#endif
    for (size_t i = 0; i < selected_atom_size; i++)
    {
        const auto & atom{ *atom_list[i] };
        const auto local_view{ AtomLocalPotentialView::RequireFor(atom) };
        auto sample_entries{
            UpdateSampleListWithFittedGaussian(atom, offset_snapshot)
        };
        auto result{ input_result_list.at(i) };
        const auto & offset_model{ offset_snapshot.at(&atom) };
        bool candidate_accepted{ false };
        try
        {
            auto candidate_result{
                EstimateLocalGaussianWithOffsetModel(
                    sample_entries,
                    local_view.GetAlphaR(),
                    options,
                    offset_model)
            };
            if (CanBuildFiniteZeroOffsetSamples(sample_entries, candidate_result.mdpde.GetModel()))
            {
                result = std::move(candidate_result);
                candidate_accepted = true;
            }
        }
        catch (const std::exception &)
        {
        }
        if (!candidate_accepted)
        {
            result.ols = WithModelOffset(result.ols, offset_model.GetOffset());
            result.mdpde = WithModelOffset(result.mdpde, offset_model.GetOffset());
        }
        const auto fitted_model{ result.mdpde.GetModel() };
        iteration_result.estimation_list[i] = fitted_model.ToVector();
        iteration_result.result_list[i] = std::move(result);
    }
    return iteration_result;
}

void ApplyLocalFittingIterationResult(
    const LocalFittingIterationResult & iteration_result,
    std::vector<AtomLocalPotentialEditor> & local_editor_list)
{
    if (local_editor_list.size() != iteration_result.result_list.size())
    {
        throw std::invalid_argument(
            "local_editor_list and iteration_result sizes are inconsistent.");
    }

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
    double offset)
{
    numeric_validation::RequireFinite(offset, "offset");
    const auto zero_offset_result{
        EstimateLocalGaussianWithOffsetModel(
            sample_entries,
            alpha_r,
            options,
            GaussianModel3D{ 0.0, 1.0, 0.0 })
    };
    if (offset == 0.0)
    {
        return zero_offset_result;
    }
    const auto offset_model{ zero_offset_result.mdpde.GetModel().WithOffset(offset) };
    return EstimateLocalGaussianWithOffsetModel(sample_entries, alpha_r, options, offset_model);
}

LocalGaussianResult EstimateLocalGaussianWithOffset(
    const LocalPotentialSampleList & sample_entries,
    double alpha_r,
    const FitOptions & options,
    double offset_initial)
{
    numeric_validation::RequireFiniteNonNegativeRange(
        options.distance_min, options.distance_max, "fit range");
    numeric_validation::RequireFiniteNonNegative(alpha_r, "alpha_r");
    numeric_validation::RequireFinite(offset_initial, "offset_initial");

    auto execution_options{ MakeExecutionOptions(options) };
    auto result{
        EstimateLocalGaussianWithOffsetModel(
            sample_entries,
            alpha_r,
            options,
            GaussianModel3D{ 0.0, 1.0, 0.0 })
    };
    auto current_model{ result.mdpde.GetModel().WithOffset(offset_initial) };
    double best_offset{ current_model.GetOffset() };
    double best_error{ std::numeric_limits<double>::infinity() };
    auto best_result{ result };
    bool has_best_result{ false };
    auto max_iterations{ execution_options.max_iterations };
    auto tolerance{ execution_options.tolerance };
    for (int t = 0; t < max_iterations; t++)
    {
        const auto offset{ current_model.GetOffset() };
        result = EstimateLocalGaussianWithOffsetModel(sample_entries, alpha_r, options, current_model);
        const auto raw_offset{
            EstimateResidualOffsetParameter(sample_entries, *result.fit_result, offset)
        };
        const auto candidate_model{ result.mdpde.GetModel().WithOffset(raw_offset) };
        const auto error{
            (candidate_model.ToVector() - current_model.ToVector()).norm()
        };
        if (error < best_error)
        {
            best_offset = offset;
            best_error = error;
            best_result = result;
            has_best_result = true;
        }
        if (error < tolerance)
        {
            break;
        }

        if (t + 1 == max_iterations)
        {
            result = has_best_result ?
                best_result :
                EstimateLocalGaussian(sample_entries, alpha_r, options, best_offset);
            if (!options.quiet_mode)
            {
                Logger::Log(LogLevel::Debug,
                    "Maximum iterations reached in local Gaussian estimation with offset; "
                    "refitting at best fixed-point candidate with error = " +
                    std::to_string(best_error) + ".");
            }
            break;
        }

        const auto damped_offset{ offset + kOffsetDampingFactor * (raw_offset - offset) };
        current_model = result.mdpde.GetModel().WithOffset(damped_offset);
    }
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
    const auto fit_result_list{ BuildMemberFitResultList(dataset_list, member_result_list, options) };
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

void RunLocalAlphaTraining(ModelObject & model_object, const FitOptions & options)
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
        const auto & group_atom_list{
            analysis_view.GetAtomObjectList(group_key)
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
    const auto group_key_list{ analysis_view.CollectAtomGroupKeys() };

    std::vector<std::vector<LocalGaussianResult>> member_result_list;
    member_result_list.reserve(group_key_list.size());
    for (const auto group_key : group_key_list)
    {
        const auto & group_atom_list{
            analysis_view.GetAtomObjectList(group_key)
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
    for (const auto group_key : analysis_view.CollectAtomGroupKeys())
    {
        analysis.SetAtomGroupAlphaG(group_key, alpha_g);
    }
}

void RunFirstStageLocalFitting(ModelObject & model_object, const FitOptions & options)
{
    const auto selected_atom_size{ model_object.GetSelectedAtomCount() };
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    auto local_editor_list{ BuildSelectedAtomLocalEditors(model_object) };
    std::atomic<size_t> atom_count{ 0 };
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info,
            "Run first-stage local atom fitting for " +
            std::to_string(selected_atom_size) + " atoms.");
    }

#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(options.thread_size)
#endif
    for (size_t i = 0; i < selected_atom_size; i++)
    {
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom_list[i]) };
        auto sample_entries{ local_view.GetSamplingEntries() };
        const auto result{
            EstimateLocalGaussianWithOffset(
                sample_entries, local_view.GetAlphaR(), options, 0.0)
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
    for (size_t iter = 0; iter < kLocalFittingMaximumIterations; iter++)
    {
        auto iteration_result{
            RunLocalFittingIteration(
                atom_list,
                previous_estimation_list,
                previous_result_list,
                options)
        };
        ApplyLocalFittingUnderRelaxation(iteration_result, previous_estimation_list, options.relaxation_factor);
        const auto change_stats{
            CalculateLocalFittingParameterChangeStats(
                iteration_result.estimation_list, previous_estimation_list)
        };
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
                << ", percentile amplitude change = "
                << change_stats.amplitude_change_percentile
                << ", percentile width change = "
                << change_stats.width_change_percentile
                << ", percentile offset change = "
                << change_stats.offset_change_percentile;
            Logger::ProgressLine(progress_message.str());
        }

        const auto converged{ IsLocalFittingParameterChangeConverged(change_stats) };
        if (converged)
        {
            ApplyLocalFittingIterationResult(iteration_result, local_editor_list);
            if (!options.quiet_mode)
            {
                Logger::FinishProgressLine();
                Logger::Log(LogLevel::Info,
                    "Converged after " + std::to_string(iter + 1) +
                    " iterations with percentile amplitude change = " +
                    std::to_string(change_stats.amplitude_change_percentile) +
                    ", percentile width change = " +
                    std::to_string(change_stats.width_change_percentile) +
                    ", and percentile offset change = " +
                    std::to_string(change_stats.offset_change_percentile) + ".");
            }
            break;
        }

        if (iter + 1 == kLocalFittingMaximumIterations)
        {
            ApplyLocalFittingIterationResult(best_iteration_result, local_editor_list);
            if (!options.quiet_mode)
            {
                Logger::FinishProgressLine();
                Logger::Log(LogLevel::Warning,
                    "Reached maximum iteration size; refitting at best fixed-point candidate "
                    "with percentile amplitude change = " +
                    std::to_string(best_change_stats.amplitude_change_percentile) +
                    ", percentile width change = " +
                    std::to_string(best_change_stats.width_change_percentile) +
                    ", and percentile offset change = " +
                    std::to_string(best_change_stats.offset_change_percentile));
            }
        }
        previous_estimation_list = std::move(iteration_result.estimation_list);
        previous_result_list = std::move(iteration_result.result_list);
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
    std::atomic<size_t> key_count{ 0 };

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
            sample_entries_list.emplace_back(UpdateSampleListWithFittedGaussian(*atom));
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

void RunPotentialFittingWorkflow(ModelObject & model_object, const FitOptions & options)
{
    RunLocalAlphaTraining(model_object, options);
    RunLocalPotentialFitting(model_object, options);
    RunGroupAlphaTraining(model_object, options);
    RunGroupPotentialFitting(model_object, options);
}

} // namespace rhbm_gem::core
