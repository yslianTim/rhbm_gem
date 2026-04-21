#include <rhbm_gem/utils/hrl/HRLDataTransform.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace
{
int CheckedCastToInt(std::size_t value)
{
    if (value > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::overflow_error("data_array size exceeds maximum int");
    }
    return static_cast<int>(value);
}
} // namespace

HRLMemberDataset HRLDataTransform::BuildMemberDataset(
    const SeriesPointList & series_point_list,
    bool quiet_mode)
{
    (void)quiet_mode;
    if (series_point_list.empty())
    {
        throw std::invalid_argument("series_point_list must not be empty.");
    }

    const auto expected_basis_size{ series_point_list.front().GetBasisSize() };
    const auto basis_size{ CheckedCastToInt(expected_basis_size) };
    rhbm_gem::NumericValidation::RequirePositive(basis_size, "basis_size");
    const auto data_size{ CheckedCastToInt(series_point_list.size()) };

    HRLMemberDataset dataset;
    dataset.X = HRLDesignMatrix::Zero(data_size, basis_size);
    dataset.y = HRLResponseVector::Zero(data_size);
    dataset.score = HRLScoreVector::Zero(data_size);
    for (int i = 0; i < data_size; i++)
    {
        const auto & point{ series_point_list.at(static_cast<std::size_t>(i)) };
        if (point.GetBasisSize() != expected_basis_size)
        {
            throw std::invalid_argument("All data entries must share the same basis size.");
        }
        if (!std::isfinite(point.response) || !std::isfinite(point.score))
        {
            throw std::invalid_argument("Member dataset contains non-finite value.");
        }
        for (int j = 0; j < basis_size; j++)
        {
            const auto value{ point.GetBasisValue(static_cast<std::size_t>(j)) };
            if (!std::isfinite(value))
            {
                throw std::invalid_argument("Member dataset contains non-finite value.");
            }
            dataset.X(i, j) = value;
        }
        dataset.y(i) = point.response;
        dataset.score(i) = point.score;
    }

    return dataset;
}

HRLBetaMatrix HRLDataTransform::BuildBetaMatrix(
    const std::vector<HRLBetaVector> & beta_list,
    bool quiet_mode)
{
    (void)quiet_mode;
    if (beta_list.empty())
    {
        throw std::invalid_argument("beta_list must not be empty.");
    }

    const auto basis_size{ static_cast<int>(beta_list.front().rows()) };
    rhbm_gem::NumericValidation::RequirePositive(basis_size, "basis_size");
    const auto member_size{ CheckedCastToInt(beta_list.size()) };
    HRLBetaMatrix beta_matrix{ HRLBetaMatrix::Zero(basis_size, member_size) };
    for (int i = 0; i < member_size; i++)
    {
        const auto & beta_vector{ beta_list.at(static_cast<std::size_t>(i)) };
        if (beta_vector.rows() != basis_size)
        {
            throw std::invalid_argument("All beta vectors must share the same basis size.");
        }
        try
        {
            rhbm_gem::EigenValidation::RequireFinite(beta_vector, "beta");
        }
        catch (const std::invalid_argument &)
        {
            throw std::invalid_argument("beta_list contains non-finite value.");
        }
        beta_matrix.col(i) = beta_vector;
    }
    return beta_matrix;
}

HRLGroupEstimationInput HRLDataTransform::BuildGroupInput(
    const std::vector<HRLMemberDataset> & member_datasets,
    const std::vector<HRLMemberLocalEstimate> & member_estimates)
{
    const auto member_size{ member_datasets.size() };
    if (member_size == 0)
    {
        throw std::invalid_argument("member_datasets must not be empty.");
    }
    if (member_estimates.size() != member_size)
    {
        throw std::invalid_argument("Group estimation inputs must have consistent member counts.");
    }

    const auto basis_size{ static_cast<int>(member_datasets.front().X.cols()) };
    rhbm_gem::NumericValidation::RequirePositive(basis_size, "basis_size");

    HRLGroupEstimationInput input;
    input.basis_size = basis_size;
    input.member_datasets.reserve(member_size);
    input.member_estimates.reserve(member_size);

    for (std::size_t i = 0; i < member_size; i++)
    {
        const auto & dataset{ member_datasets.at(i) };
        if (dataset.X.cols() != basis_size)
        {
            throw std::invalid_argument("Member dataset basis size is inconsistent.");
        }
        if (dataset.X.rows() != dataset.y.size())
        {
            throw std::invalid_argument("Member dataset shape is inconsistent.");
        }

        const auto & estimate{ member_estimates.at(i) };
        if (estimate.beta_mdpde.rows() != basis_size)
        {
            throw std::invalid_argument("Member beta basis size is inconsistent.");
        }
        try
        {
            rhbm_gem::EigenValidation::RequireVectorSize(
                estimate.data_weight.diagonal(),
                dataset.y.size(),
                "Member data weight");
            rhbm_gem::EigenValidation::RequireVectorSize(
                estimate.data_covariance.diagonal(),
                dataset.y.size(),
                "Member data covariance");
        }
        catch (const std::invalid_argument &)
        {
            throw std::invalid_argument("Member covariance or weight size is inconsistent.");
        }

        input.member_datasets.emplace_back(dataset);
        input.member_estimates.emplace_back(estimate);
    }

    return input;
}
