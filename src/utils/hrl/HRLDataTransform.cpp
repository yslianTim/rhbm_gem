#include <rhbm_gem/utils/hrl/HRLDataTransform.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <Eigen/Core>

#include <limits>
#include <stdexcept>

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
    rhbm_gem::NumericValidation::RequireAtMost(
        expected_basis_size,
        static_cast<std::size_t>(std::numeric_limits<int>::max()),
        "basis_size");
    const auto basis_size{ static_cast<int>(expected_basis_size) };
    rhbm_gem::NumericValidation::RequirePositive(basis_size, "basis_size");
    rhbm_gem::NumericValidation::RequireAtMost(
        series_point_list.size(),
        static_cast<std::size_t>(std::numeric_limits<int>::max()),
        "series_point_list size");
    const auto data_size{ static_cast<int>(series_point_list.size()) };

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
        try
        {
            rhbm_gem::NumericValidation::RequireFinite(point.response, "response");
            rhbm_gem::NumericValidation::RequireFinite(point.score, "score");
            rhbm_gem::NumericValidation::RequireAllFinite(point.basis_list, "basis_list");
        }
        catch (const std::invalid_argument &)
        {
            throw std::invalid_argument("Member dataset contains non-finite value.");
        }

        dataset.X.row(i) = Eigen::Map<const Eigen::RowVectorXd>(
            point.basis_list.data(),
            basis_size);
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
    rhbm_gem::NumericValidation::RequireAtMost(
        beta_list.size(),
        static_cast<std::size_t>(std::numeric_limits<int>::max()),
        "beta_list size");
    const auto member_size{ static_cast<int>(beta_list.size()) };
    HRLBetaMatrix beta_matrix{ HRLBetaMatrix::Zero(basis_size, member_size) };
    for (int i = 0; i < member_size; i++)
    {
        const auto & beta_vector{ beta_list.at(static_cast<std::size_t>(i)) };
        try
        {
            rhbm_gem::EigenValidation::RequireVectorSize(beta_vector, basis_size, "beta");
        }
        catch (const std::invalid_argument &)
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
        const auto & estimate{ member_estimates.at(i) };

        if (dataset.X.cols() != basis_size)
        {
            throw std::invalid_argument("Member dataset basis size is inconsistent.");
        }

        try
        {
            rhbm_gem::EigenValidation::RequireVectorSize(
                dataset.y,
                dataset.X.rows(),
                "Member dataset response");
            rhbm_gem::EigenValidation::RequireSameSize(
                dataset.score,
                dataset.y,
                "Member dataset shape");
        }
        catch (const std::invalid_argument &)
        {
            throw std::invalid_argument("Member dataset shape is inconsistent.");
        }

        try
        {
            rhbm_gem::EigenValidation::RequireVectorSize(
                estimate.beta_mdpde,
                basis_size,
                "Member beta");
        }
        catch (const std::invalid_argument &)
        {
            throw std::invalid_argument("Member beta basis size is inconsistent.");
        }

        try
        {
            rhbm_gem::EigenValidation::RequireSameSize(
                estimate.data_weight.diagonal(),
                dataset.y,
                "Member data weight");
            rhbm_gem::EigenValidation::RequireSameSize(
                estimate.data_covariance.diagonal(),
                dataset.y,
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
