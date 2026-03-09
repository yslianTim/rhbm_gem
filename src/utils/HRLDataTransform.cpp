#include <rhbm_gem/utils/hrl/HRLDataTransform.hpp>

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
    const std::vector<Eigen::VectorXd> & data_vector,
    bool quiet_mode)
{
    (void)quiet_mode;
    if (data_vector.empty())
    {
        throw std::invalid_argument("data_vector must not be empty.");
    }

    const auto expected_rows{ static_cast<int>(data_vector.front().rows()) };
    const auto basis_size{ expected_rows - 1 };
    ValidateBasisSize(basis_size);
    const auto data_size{ CheckedCastToInt(data_vector.size()) };

    HRLMemberDataset dataset;
    dataset.X = Eigen::MatrixXd::Zero(data_size, basis_size);
    dataset.y = Eigen::VectorXd::Zero(data_size);
    for (int i = 0; i < data_size; i++)
    {
        const auto & data{ data_vector.at(static_cast<std::size_t>(i)) };
        if (data.rows() != expected_rows)
        {
            throw std::invalid_argument("All data entries must share the same basis size.");
        }
        if (!data.array().allFinite())
        {
            throw std::invalid_argument("Member dataset contains non-finite value.");
        }
        dataset.X.row(i) = data.head(basis_size);
        dataset.y(i) = data(basis_size);
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
    int basis_size,
    const std::vector<std::vector<Eigen::VectorXd>> & data_entries,
    const std::vector<Eigen::VectorXd> & beta_list,
    const std::vector<double> & sigma_square_list,
    const std::vector<HRLDiagonalMatrix> & weight_list,
    const std::vector<HRLDiagonalMatrix> & covariance_list)
{
    ValidateBasisSize(basis_size);
    const auto member_size{ data_entries.size() };
    if (member_size == 0)
    {
        throw std::invalid_argument("member_datasets must not be empty.");
    }
    if (beta_list.size() != member_size ||
        sigma_square_list.size() != member_size ||
        weight_list.size() != member_size ||
        covariance_list.size() != member_size)
    {
        throw std::invalid_argument("Group estimation inputs must have consistent member counts.");
    }

    HRLGroupEstimationInput input;
    input.basis_size = basis_size;
    input.member_datasets.reserve(member_size);
    input.member_estimates.reserve(member_size);

    for (std::size_t i = 0; i < member_size; i++)
    {
        auto dataset{ BuildMemberDataset(data_entries.at(i), true) };
        if (dataset.X.cols() != basis_size)
        {
            throw std::invalid_argument("Member dataset basis size is inconsistent.");
        }

        const auto & beta{ beta_list.at(i) };
        const auto & weight{ weight_list.at(i) };
        const auto & covariance{ covariance_list.at(i) };
        if (beta.rows() != basis_size)
        {
            throw std::invalid_argument("Member beta basis size is inconsistent.");
        }
        if (weight.diagonal().size() != dataset.y.size() ||
            covariance.diagonal().size() != dataset.y.size())
        {
            throw std::invalid_argument("Member covariance or weight size is inconsistent.");
        }

        HRLMemberLocalEstimate estimate;
        estimate.beta_mdpde = beta;
        estimate.sigma_square = sigma_square_list.at(i);
        estimate.data_weight = weight;
        estimate.data_covariance = covariance;

        input.member_datasets.emplace_back(std::move(dataset));
        input.member_estimates.emplace_back(std::move(estimate));
    }

    return input;
}
