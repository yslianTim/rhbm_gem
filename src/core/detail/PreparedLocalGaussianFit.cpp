#include "core/detail/PreparedLocalGaussianFit.hpp"

#include <cmath>
#include <stdexcept>

#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

namespace rhbm_gem::core::detail {

namespace {

double CalculateAdjustedResponse(
    double sample_response,
    double distance,
    const GaussianModel3D & offset_model)
{
    const auto evaluation{ offset_model.EvaluateAtDistance(distance) };
    return sample_response - (evaluation.response - evaluation.signal);
}

LocalGaussianResult DecodeLocalGaussianResult(
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
        GaussianModel3DWithUncertainty{
            ols_model,
            GaussianModel3DUncertainty{}
        },
        GaussianModel3DWithUncertainty{
            mdpde_model,
            GaussianModel3DUncertainty{}
        },
        std::nullopt,
        false,
        0.0,
        fit_result
    };
}

} // namespace

PreparedLocalGaussianDesign::PreparedLocalGaussianDesign(
    const LocalPotentialSampleList & sample_entries,
    double range_min,
    double range_max)
    : m_source_sample_count{ sample_entries.size() }
{
    numeric_validation::RequireFiniteNonNegativeRange(range_min, range_max, "fit range");
    m_row_list.reserve(sample_entries.size());
    for (std::size_t sample_index = 0; sample_index < sample_entries.size(); sample_index++)
    {
        const auto distance{ sample_entries.at(sample_index).point.distance };
        if (distance < range_min || distance > range_max)
        {
            continue;
        }
        m_row_list.emplace_back(Row{
            sample_index,
            distance
        });
    }

    const auto row_count{ static_cast<Eigen::Index>(m_row_list.size()) };
    m_design_matrix = RHBMDesignMatrix::Zero(row_count, 2);
    for (Eigen::Index row = 0; row < row_count; row++)
    {
        const auto distance{ m_row_list.at(static_cast<std::size_t>(row)).distance };
        m_design_matrix(row, 0) = 1.0;
        m_design_matrix(row, 1) = -0.5 * distance * distance;
    }
}

RHBMMemberDataset PreparedLocalGaussianDesign::BuildDataset(
    const std::vector<double> & sample_response_list,
    const GaussianModel3D & offset_model) const
{
    numeric_validation::RequireFinite(offset_model.GetOffset(), "offset");
    if (sample_response_list.size() != m_source_sample_count ||
        m_design_matrix.rows() != static_cast<Eigen::Index>(m_row_list.size()) ||
        m_design_matrix.cols() != 2)
    {
        throw std::invalid_argument("Prepared local Gaussian design inputs are inconsistent.");
    }

    const auto candidate_count{ static_cast<Eigen::Index>(m_row_list.size()) };
    RHBMMemberDataset dataset;
    dataset.X = RHBMDesignMatrix::Zero(candidate_count, 2);
    dataset.y = RHBMResponseVector::Zero(candidate_count);
    Eigen::Index retained_count{ 0 };
    for (std::size_t row = 0; row < m_row_list.size(); row++)
    {
        const auto & design_row{ m_row_list.at(row) };
        if (design_row.source_sample_index >= sample_response_list.size())
        {
            throw std::invalid_argument("Prepared local Gaussian sample index is out of range.");
        }
        const auto adjusted_response{
            CalculateAdjustedResponse(
                sample_response_list.at(design_row.source_sample_index),
                design_row.distance,
                offset_model)
        };
        if (adjusted_response <= 0.0) continue;
        numeric_validation::RequireFinite(adjusted_response, "response", "Member dataset contains non-finite value.");
        dataset.X.row(retained_count) = m_design_matrix.row(static_cast<Eigen::Index>(row));
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

LocalGaussianResult PreparedLocalGaussianDesign::Estimate(
    const std::vector<double> & sample_response_list,
    double alpha_r,
    int thread_size,
    const GaussianModel3D & offset_model) const
{
    numeric_validation::RequireFiniteNonNegative(alpha_r, "alpha_r");
    auto dataset{ BuildDataset(sample_response_list, offset_model) };
    const auto result{
        rhbm_helper::EstimateBetaMDPDE(
            alpha_r,
            dataset,
            RHBMExecutionOptions{ .thread_size = thread_size })
    };
    return DecodeLocalGaussianResult(alpha_r, result, offset_model.GetOffset());
}

} // namespace rhbm_gem::core::detail
