#include <rhbm_gem/utils/hrl/HRLDataTransform.hpp>

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

void ValidateBasisSize(int basis_size)
{
    if (basis_size <= 0)
    {
        throw std::invalid_argument("basis_size must be positive.");
    }
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
    ValidateBasisSize(basis_size);
    const auto data_size{ CheckedCastToInt(series_point_list.size()) };

    HRLMemberDataset dataset;
    dataset.X = Eigen::MatrixXd::Zero(data_size, basis_size);
    dataset.y = Eigen::VectorXd::Zero(data_size);
    dataset.w = Eigen::VectorXd::Zero(data_size);
    for (int i = 0; i < data_size; i++)
    {
        const auto & point{ series_point_list.at(static_cast<std::size_t>(i)) };
        if (point.GetBasisSize() != expected_basis_size)
        {
            throw std::invalid_argument("All data entries must share the same basis size.");
        }
        if (!std::isfinite(point.response) || !std::isfinite(point.weight))
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
        dataset.w(i) = point.weight;
    }

    return dataset;
}

Eigen::MatrixXd HRLDataTransform::BuildBetaMatrix(
    const std::vector<Eigen::VectorXd> & beta_list,
    bool quiet_mode)
{
    (void)quiet_mode;
    if (beta_list.empty())
    {
        throw std::invalid_argument("beta_list must not be empty.");
    }

    const auto basis_size{ static_cast<int>(beta_list.front().rows()) };
    ValidateBasisSize(basis_size);
    const auto member_size{ CheckedCastToInt(beta_list.size()) };
    Eigen::MatrixXd beta_array{ Eigen::MatrixXd::Zero(basis_size, member_size) };
    for (int i = 0; i < member_size; i++)
    {
        const auto & beta{ beta_list.at(static_cast<std::size_t>(i)) };
        if (beta.rows() != basis_size)
        {
            throw std::invalid_argument("All beta vectors must share the same basis size.");
        }
        if (!beta.array().allFinite())
        {
            throw std::invalid_argument("beta_list contains non-finite value.");
        }
        beta_array.col(i) = beta;
    }
    return beta_array;
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
    ValidateBasisSize(basis_size);

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
        if (estimate.data_weight.diagonal().size() != dataset.y.size() ||
            estimate.data_covariance.diagonal().size() != dataset.y.size())
        {
            throw std::invalid_argument("Member covariance or weight size is inconsistent.");
        }

        input.member_datasets.emplace_back(dataset);
        input.member_estimates.emplace_back(estimate);
    }

    return input;
}
