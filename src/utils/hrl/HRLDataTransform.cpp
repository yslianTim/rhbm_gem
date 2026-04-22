#include <rhbm_gem/utils/hrl/HRLDataTransform.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <Eigen/Core>

#include <limits>
#include <stdexcept>

HRLMemberDataset HRLDataTransform::BuildMemberDataset(const SeriesPointList & series_point_list)
{
    if (series_point_list.empty())
    {
        throw std::invalid_argument("series_point_list must not be empty.");
    }

    const auto basis_size{ static_cast<int>(series_point_list.front().GetBasisSize()) };
    rhbm_gem::numeric_validation::RequirePositive(basis_size, "basis_size");
    rhbm_gem::numeric_validation::RequireAtMost(
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
        if (static_cast<int>(point.GetBasisSize()) != basis_size)
        {
            throw std::invalid_argument("All data entries must share the same basis size.");
        }
        rhbm_gem::numeric_validation::RequireFinite(
            point.response,
            "response",
            "Member dataset contains non-finite value.");
        rhbm_gem::numeric_validation::RequireFinite(
            point.score,
            "score",
            "Member dataset contains non-finite value.");
        rhbm_gem::numeric_validation::RequireAllFinite(
            point.basis_list,
            "basis_list",
            "Member dataset contains non-finite value.");

        dataset.X.row(i) = Eigen::Map<const Eigen::RowVectorXd>(
            point.basis_list.data(), basis_size);
        dataset.y(i) = point.response;
        dataset.score(i) = point.score;
    }

    return dataset;
}

HRLBetaMatrix HRLDataTransform::BuildBetaMatrix(const std::vector<HRLBetaVector> & beta_list)
{
    if (beta_list.empty())
    {
        throw std::invalid_argument("beta_list must not be empty.");
    }

    const auto basis_size{ static_cast<int>(beta_list.front().rows()) };
    rhbm_gem::numeric_validation::RequirePositive(basis_size, "basis_size");
    rhbm_gem::numeric_validation::RequireAtMost(
        beta_list.size(),
        static_cast<std::size_t>(std::numeric_limits<int>::max()),
        "beta_list size");
    const auto member_size{ static_cast<int>(beta_list.size()) };
    HRLBetaMatrix beta_matrix{ HRLBetaMatrix::Zero(basis_size, member_size) };
    for (int i = 0; i < member_size; i++)
    {
        const auto & beta_vector{ beta_list.at(static_cast<std::size_t>(i)) };
        rhbm_gem::eigen_validation::RequireVectorSize(
            beta_vector,
            basis_size,
            "beta",
            "All beta vectors must share the same basis size.");
        rhbm_gem::eigen_validation::RequireFinite(
            beta_vector,
            "beta",
            "beta_list contains non-finite value.");
        beta_matrix.col(i) = beta_vector;
    }
    return beta_matrix;
}

HRLGroupEstimationInput HRLDataTransform::BuildGroupInput(
    const std::vector<HRLMemberDataset> & member_datasets,
    const std::vector<HRLBetaEstimateResult> & member_fit_results)
{
    const auto member_size{ member_datasets.size() };
    if (member_size == 0)
    {
        throw std::invalid_argument("member_datasets must not be empty.");
    }
    if (member_fit_results.size() != member_size)
    {
        throw std::invalid_argument("Group estimation inputs must have consistent member counts.");
    }

    const auto basis_size{ static_cast<int>(member_datasets.front().X.cols()) };
    rhbm_gem::numeric_validation::RequirePositive(basis_size, "basis_size");

    HRLGroupEstimationInput input;
    input.basis_size = basis_size;
    input.member_datasets.reserve(member_size);
    input.member_fit_results.reserve(member_size);

    for (std::size_t i = 0; i < member_size; i++)
    {
        const auto & dataset{ member_datasets.at(i) };
        const auto & fit_result{ member_fit_results.at(i) };

        if (dataset.X.cols() != basis_size)
        {
            throw std::invalid_argument("Member dataset basis size is inconsistent.");
        }
        rhbm_gem::eigen_validation::RequireVectorSize(
            dataset.y,
            dataset.X.rows(),
            "Member dataset response",
            "Member dataset shape is inconsistent.");
        rhbm_gem::eigen_validation::RequireSameSize(
            dataset.score,
            dataset.y,
            "Member dataset shape",
            "Member dataset shape is inconsistent.");
        rhbm_gem::eigen_validation::RequireVectorSize(
            fit_result.beta_mdpde,
            basis_size,
            "Member fit beta",
            "Member beta basis size is inconsistent.");
        rhbm_gem::eigen_validation::RequireSameSize(
            fit_result.data_weight.diagonal(),
            dataset.y,
            "Member fit data weight",
            "Member covariance or weight size is inconsistent.");
        rhbm_gem::eigen_validation::RequireSameSize(
            fit_result.data_covariance.diagonal(),
            dataset.y,
            "Member fit data covariance",
            "Member covariance or weight size is inconsistent.");

        input.member_datasets.emplace_back(dataset);
        input.member_fit_results.emplace_back(fit_result);
    }

    return input;
}
