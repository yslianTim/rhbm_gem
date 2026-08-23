#pragma once

#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

namespace rhbm_gem::core::detail {

struct LocalGaussianDesignTemplate
{
    std::size_t source_sample_count{ 0 };
    std::vector<std::size_t> source_sample_index_list{};
    std::vector<double> distance_list{};
    RHBMDesignMatrix design_matrix{};
};

inline RHBMExecutionOptions MakeExecutionOptions(const FitOptions & options)
{
    RHBMExecutionOptions execution_options;
    execution_options.quiet_mode = false;
    execution_options.thread_size = options.thread_size;
    return execution_options;
}

inline float CalculateAdjustedResponse(
    double sample_response,
    double distance,
    const GaussianModel3D & offset_model)
{
    const auto evaluation{ offset_model.EvaluateAtDistance(distance) };
    return static_cast<float>(sample_response - (evaluation.response - evaluation.signal));
}

inline LocalGaussianResult DecodeLocalGaussianResult(
    double alpha_r,
    const RHBMBetaEstimateResult & fit_result,
    double offset)
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

inline LocalGaussianDesignTemplate BuildLocalGaussianDesignTemplate(
    const LocalPotentialSampleList & sample_entries,
    double range_min,
    double range_max)
{
    numeric_validation::RequireFiniteNonNegativeRange(range_min, range_max, "fit range");

    LocalGaussianDesignTemplate design_template;
    design_template.source_sample_count = sample_entries.size();
    design_template.source_sample_index_list.reserve(sample_entries.size());
    design_template.distance_list.reserve(sample_entries.size());
    for (std::size_t sample_index = 0; sample_index < sample_entries.size(); sample_index++)
    {
        const auto distance{
            static_cast<double>(sample_entries.at(sample_index).point.distance)
        };
        if (distance < range_min || distance > range_max) continue;
        design_template.source_sample_index_list.emplace_back(sample_index);
        design_template.distance_list.emplace_back(distance);
    }

    const auto row_count{
        static_cast<Eigen::Index>(design_template.distance_list.size())
    };
    design_template.design_matrix = RHBMDesignMatrix::Zero(row_count, 2);
    for (Eigen::Index row = 0; row < row_count; row++)
    {
        const auto distance{
            design_template.distance_list.at(static_cast<std::size_t>(row))
        };
        design_template.design_matrix(row, 0) = 1.0;
        design_template.design_matrix(row, 1) = -0.5 * distance * distance;
    }
    return design_template;
}

inline RHBMMemberDataset BuildLocalGaussianPreparedDataset(
    const LocalGaussianDesignTemplate & design_template,
    const std::vector<double> & sample_response_list,
    const GaussianModel3D & offset_model)
{
    numeric_validation::RequireFinite(offset_model.GetOffset(), "offset");
    if (sample_response_list.size() != design_template.source_sample_count ||
        design_template.source_sample_index_list.size() != design_template.distance_list.size() ||
        design_template.design_matrix.rows() != static_cast<Eigen::Index>(design_template.distance_list.size()) ||
        design_template.design_matrix.cols() != 2)
    {
        throw std::invalid_argument("Prepared local Gaussian design inputs are inconsistent.");
    }

    const auto candidate_count{
        static_cast<Eigen::Index>(design_template.distance_list.size())
    };
    RHBMMemberDataset dataset;
    dataset.X = RHBMDesignMatrix::Zero(candidate_count, 2);
    dataset.y = RHBMResponseVector::Zero(candidate_count);
    Eigen::Index retained_count{ 0 };
    for (std::size_t row = 0; row < design_template.distance_list.size(); row++)
    {
        const auto sample_index{ design_template.source_sample_index_list.at(row) };
        if (sample_index >= sample_response_list.size())
        {
            throw std::invalid_argument("Prepared local Gaussian sample index is out of range.");
        }
        const auto adjusted_response{
            static_cast<double>(CalculateAdjustedResponse(
                sample_response_list.at(sample_index),
                design_template.distance_list.at(row),
                offset_model))
        };
        if (adjusted_response <= 0.0) continue;
        numeric_validation::RequireFinite(
            adjusted_response,
            "response",
            "Member dataset contains non-finite value.");
        dataset.X.row(retained_count) = design_template.design_matrix.row(
            static_cast<Eigen::Index>(row));
        dataset.y(retained_count) = std::log(adjusted_response);
        retained_count++;
    }

    if (retained_count == 0)
    {
        dataset.X = RHBMDesignMatrix::Zero(1, 2);
        dataset.y = RHBMResponseVector::Zero(1);
    }
    else
    {
        dataset.X.conservativeResize(retained_count, 2);
        dataset.y.conservativeResize(retained_count);
    }

    return dataset;
}

inline LocalGaussianResult EstimateLocalGaussianPrepared(
    const LocalGaussianDesignTemplate & design_template,
    const std::vector<double> & sample_response_list,
    double alpha_r,
    const FitOptions & options,
    const GaussianModel3D & offset_model)
{
    numeric_validation::RequireFiniteNonNegative(alpha_r, "alpha_r");
    auto dataset{
        BuildLocalGaussianPreparedDataset(design_template, sample_response_list, offset_model)
    };
    const auto result{
        rhbm_helper::EstimateBetaMDPDE(alpha_r, dataset, MakeExecutionOptions(options))
    };
    return DecodeLocalGaussianResult(alpha_r, result, offset_model.GetOffset());
}

} // namespace rhbm_gem::core::detail
