#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>

namespace rg = rhbm_gem;

namespace
{
template <typename Callable>
void ExpectInvalidArgumentMessage(
    Callable && callable,
    const std::string & expected_message)
{
    try
    {
        callable();
        FAIL() << "Expected std::invalid_argument";
    }
    catch (const std::invalid_argument & ex)
    {
        EXPECT_EQ(expected_message, ex.what());
    }
}

Eigen::VectorXd MakeVector(std::initializer_list<double> values)
{
    Eigen::VectorXd result(static_cast<Eigen::Index>(values.size()));
    Eigen::Index index{ 0 };
    for (double value : values)
    {
        result(index++) = value;
    }
    return result;
}

rg::RHBMDiagonalMatrix MakeDiagonal(std::initializer_list<double> values)
{
    rg::RHBMDiagonalMatrix result(static_cast<Eigen::Index>(values.size()));
    result.diagonal() = MakeVector(values);
    return result;
}

rg::RHBMMemberDataset MakeDataset(
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & y)
{
    return rg::RHBMMemberDataset{ X, y };
}

using DatasetRow = std::vector<double>;

DatasetRow MakeDatasetRow(std::initializer_list<double> values)
{
    return DatasetRow(values);
}

rg::RHBMMemberDataset MakeDataset(std::initializer_list<DatasetRow> rows)
{
    const auto data_size{ static_cast<Eigen::Index>(rows.size()) };
    const auto basis_size{ static_cast<Eigen::Index>(rows.begin()->size() - 1) };
    rg::RHBMMemberDataset dataset;
    dataset.X = rg::RHBMDesignMatrix::Zero(data_size, basis_size);
    dataset.y = rg::RHBMResponseVector::Zero(data_size);

    Eigen::Index row_index{ 0 };
    for (const auto & row : rows)
    {
        for (Eigen::Index column = 0; column < basis_size; column++)
        {
            dataset.X(row_index, column) = row.at(static_cast<std::size_t>(column));
        }
        dataset.y(row_index) = row.back();
        row_index++;
    }
    return dataset;
}

rg::RHBMBetaEstimateResult MakeFitResult(
    std::initializer_list<double> beta_values,
    double sigma_square,
    std::initializer_list<double> weight_values,
    std::initializer_list<double> covariance_values,
    rg::RHBMEstimationStatus status = rg::RHBMEstimationStatus::SUCCESS)
{
    rg::RHBMBetaEstimateResult fit_result;
    fit_result.status = status;
    fit_result.beta_mdpde = MakeVector(beta_values);
    fit_result.sigma_square = sigma_square;
    fit_result.data_weight = MakeDiagonal(weight_values);
    fit_result.data_covariance = MakeDiagonal(covariance_values);
    return fit_result;
}

int AlternateThreadCount(int thread_count)
{
    return (thread_count == 1) ? 2 : 1;
}
} // namespace

TEST(RHBMHelperTest, BuildMemberDatasetTransformsSamplesToLogQuadraticDataset)
{
    const LocalPotentialSampleList sampling_entries{
        LocalPotentialSample{ 2.0F, SamplingPoint{ 1.0F } },
        LocalPotentialSample{ 4.0F, SamplingPoint{ 2.0F } }
    };

    const auto dataset{ rhbm_gem::rhbm_helper::BuildMemberDataset(sampling_entries, 0.0, 3.0) };

    Eigen::MatrixXd expected_X(2, 2);
    expected_X << 1.0, -0.5,
                  1.0, -2.0;
    EXPECT_TRUE(dataset.X.isApprox(expected_X, 1e-12));
    EXPECT_TRUE(dataset.y.isApprox(MakeVector({ std::log(2.0), std::log(4.0) }), 1e-12));
}

TEST(RHBMHelperTest, BuildMemberDatasetTransformsSamplesToDifferentialDataset)
{
    const LocalPotentialSampleList sampling_entries{
        LocalPotentialSample{ 2.0F, SamplingPoint{ 0.5F } },
        LocalPotentialSample{ 4.0F, SamplingPoint{ 1.0F } },
        LocalPotentialSample{ 6.0F, SamplingPoint{ 1.0F } },
        LocalPotentialSample{ 8.0F, SamplingPoint{ 2.0F } }
    };

    const auto dataset{
        rhbm_gem::rhbm_helper::BuildMemberDataset(
            sampling_entries,
            1.0,
            2.0,
            rhbm_gem::LocalGaussianFitModel::DifferentialMethod)
    };

    Eigen::MatrixXd expected_X(3, 3);
    expected_X << 1.0, -1.75, 1.0,
                  1.0, -1.75, 1.0,
                  1.0, -12.25, 4.0;
    EXPECT_TRUE(dataset.X.isApprox(expected_X, 1e-12));
    EXPECT_TRUE(dataset.y.isApprox(MakeVector({ 4.0, 6.0, 8.0 }), 1e-12));
}

TEST(RHBMHelperTest, BuildMemberDatasetTransformsSeriesPoints)
{
    const SeriesPointList series_point_list{
        SeriesPoint{ std::vector<double>{ 1.0, 0.25 }, 2.0 },
        SeriesPoint{ std::vector<double>{ 1.0, 0.50 }, -1.0 }
    };

    const auto dataset{ rhbm_gem::rhbm_helper::BuildMemberDataset(series_point_list) };

    Eigen::MatrixXd expected_X(2, 2);
    expected_X << 1.0, 0.25,
                  1.0, 0.50;
    EXPECT_TRUE(dataset.X.isApprox(expected_X, 1e-12));
    EXPECT_TRUE(dataset.y.isApprox(MakeVector({ 2.0, -1.0 }), 1e-12));
}

TEST(RHBMHelperTest, BuildMemberDatasetFiltersByRangeAndPositiveResponse)
{
    const LocalPotentialSampleList sampling_entries{
        LocalPotentialSample{ 2.0F, SamplingPoint{ 0.5F } },
        LocalPotentialSample{ 3.0F, SamplingPoint{ 1.0F } },
        LocalPotentialSample{ -1.0F, SamplingPoint{ 2.0F } },
        LocalPotentialSample{ 5.0F, SamplingPoint{ 3.0F } },
        LocalPotentialSample{ 7.0F, SamplingPoint{ 3.5F } }
    };

    const auto dataset{ rhbm_gem::rhbm_helper::BuildMemberDataset(sampling_entries, 1.0, 3.0) };

    Eigen::MatrixXd expected_X(2, 2);
    expected_X << 1.0, -0.5,
                  1.0, -4.5;
    EXPECT_TRUE(dataset.X.isApprox(expected_X, 1e-12));
    EXPECT_TRUE(dataset.y.isApprox(MakeVector({ std::log(3.0), std::log(5.0) }), 1e-12));
}

TEST(RHBMHelperTest, BuildMemberDatasetReturnsZeroFallbackWhenNoSampleIsValid)
{
    const LocalPotentialSampleList sampling_entries{
        LocalPotentialSample{ 2.0F, SamplingPoint{ 0.5F } },
        LocalPotentialSample{ 0.0F, SamplingPoint{ 1.5F } }
    };

    const auto dataset{ rhbm_gem::rhbm_helper::BuildMemberDataset(sampling_entries, 1.0, 2.0) };

    Eigen::MatrixXd expected_X{ Eigen::MatrixXd::Zero(1, 2) };
    EXPECT_TRUE(dataset.X.isApprox(expected_X, 1e-12));
    EXPECT_TRUE(dataset.y.isApprox(MakeVector({ 0.0 }), 1e-12));
}

TEST(RHBMHelperTest, BuildMemberDatasetRejectsInvalidRange)
{
    const LocalPotentialSampleList sampling_entries{
        LocalPotentialSample{ 2.0F, SamplingPoint{ 1.0F } }
    };

    EXPECT_THROW(
        rhbm_gem::rhbm_helper::BuildMemberDataset(sampling_entries, 2.0, 1.0),
        std::invalid_argument);
}

TEST(RHBMHelperTest, BuildMemberDatasetRejectsNonFiniteValues)
{
    const LocalPotentialSampleList distance_nan{
        LocalPotentialSample{ 2.0F, SamplingPoint{ std::numeric_limits<float>::quiet_NaN() } }
    };
    const LocalPotentialSampleList response_nan{
        LocalPotentialSample{ std::numeric_limits<float>::quiet_NaN(), SamplingPoint{ 1.0F } }
    };

    EXPECT_THROW(
        rhbm_gem::rhbm_helper::BuildMemberDataset(distance_nan, 0.0, 2.0),
        std::invalid_argument);
    EXPECT_THROW(
        rhbm_gem::rhbm_helper::BuildMemberDataset(response_nan, 0.0, 2.0),
        std::invalid_argument);
}

TEST(RHBMHelperTest, BuildBetaMatrixMapsVectorsToColumns)
{
    const std::vector<Eigen::VectorXd> beta_list{
        MakeVector({ 1.0, 2.0 }),
        MakeVector({ 3.0, 4.0 }),
        MakeVector({ 5.0, 6.0 })
    };

    const auto beta_matrix{ rhbm_gem::rhbm_helper::BuildBetaMatrix(beta_list) };

    Eigen::MatrixXd expected(2, 3);
    expected << 1.0, 3.0, 5.0,
                2.0, 4.0, 6.0;
    EXPECT_TRUE(beta_matrix.isApprox(expected, 1e-12));
}

TEST(RHBMHelperTest, BuildBetaMatrixRejectsEmptyInput)
{
    const std::vector<Eigen::VectorXd> beta_list;
    EXPECT_THROW(rhbm_gem::rhbm_helper::BuildBetaMatrix(beta_list), std::invalid_argument);
}

TEST(RHBMHelperTest, BuildBetaMatrixRejectsEmptyBasis)
{
    const std::vector<Eigen::VectorXd> beta_list{
        Eigen::VectorXd{}
    };
    EXPECT_THROW(rhbm_gem::rhbm_helper::BuildBetaMatrix(beta_list), std::invalid_argument);
}

TEST(RHBMHelperTest, BuildBetaMatrixRejectsNonFiniteValue)
{
    const std::vector<Eigen::VectorXd> beta_list{
        MakeVector({ 1.0, std::numeric_limits<double>::quiet_NaN() })
    };

    EXPECT_THROW(rhbm_gem::rhbm_helper::BuildBetaMatrix(beta_list), std::invalid_argument);
}

TEST(RHBMHelperTest, BuildBetaMatrixRejectsInconsistentBasisSize)
{
    const std::vector<Eigen::VectorXd> beta_list{
        MakeVector({ 1.0, 2.0 }),
        MakeVector({ 3.0, 4.0, 5.0 })
    };

    EXPECT_THROW(rhbm_gem::rhbm_helper::BuildBetaMatrix(beta_list), std::invalid_argument);
}

TEST(RHBMHelperTest, BuildGroupInputBuildsStructuredRequest)
{
    const auto input{
        rhbm_gem::rhbm_helper::BuildGroupInput(
            {
                MakeDataset({ MakeDatasetRow({ 1.0, 0.0, 1.0 }), MakeDatasetRow({ 1.0, 1.0, 3.0 }) }),
                MakeDataset({ MakeDatasetRow({ 1.0, 0.0, 2.0 }), MakeDatasetRow({ 1.0, 1.0, 4.0 }) })
            },
            {
                MakeFitResult({ 1.0, 2.0 }, 0.25, { 1.0, 1.0 }, { 0.25, 0.25 }),
                MakeFitResult({ 2.0, 2.0 }, 0.5, { 1.0, 1.0 }, { 0.5, 0.5 })
            }
        )
    };

    ASSERT_EQ(2u, input.member_datasets.size());
    ASSERT_EQ(2u, input.member_fit_results.size());
    EXPECT_EQ(2, input.basis_size);
    EXPECT_TRUE(input.member_fit_results[0].beta_mdpde.isApprox(MakeVector({ 1.0, 2.0 }), 1e-12));
}

TEST(RHBMHelperTest, BuildGroupInputRejectsMismatchedMemberCounts)
{
    EXPECT_THROW(
        rhbm_gem::rhbm_helper::BuildGroupInput(
            { MakeDataset({ MakeDatasetRow({ 1.0, 0.0, 1.0 }) }) },
            {
                MakeFitResult({ 1.0, 2.0 }, 0.25, { 1.0 }, { 1.0 }),
                MakeFitResult({ 3.0, 4.0 }, 0.5, { 1.0 }, { 1.0 })
            }
        ),
        std::invalid_argument
    );
}

TEST(RHBMHelperTest, BuildGroupInputRejectsInconsistentWeightSize)
{
    EXPECT_THROW(
        rhbm_gem::rhbm_helper::BuildGroupInput(
            {
                MakeDataset({ MakeDatasetRow({ 1.0, 0.0, 1.0 }), MakeDatasetRow({ 1.0, 1.0, 3.0 }) })
            },
            {
                MakeFitResult({ 1.0, 2.0 }, 0.25, { 1.0 }, { 0.25, 0.25 })
            }
        ),
        std::invalid_argument
    );
}

TEST(RHBMHelperTest, BuildGroupInputRejectsInconsistentMemberBetaBasisSize)
{
    EXPECT_THROW(
        rhbm_gem::rhbm_helper::BuildGroupInput(
            {
                MakeDataset({
                    MakeDatasetRow({ 1.0, 0.0, 1.0 }),
                    MakeDatasetRow({ 1.0, 1.0, 3.0 })
                })
            },
            {
                MakeFitResult({ 1.0 }, 0.25, { 1.0, 1.0 }, { 0.25, 0.25 })
            }
        ),
        std::invalid_argument
    );
}

TEST(RHBMHelperTest, EstimateGroupSingleMemberUsesFallbackResult)
{
    const auto input{
        rhbm_gem::rhbm_helper::BuildGroupInput(
            {
                MakeDataset({ MakeDatasetRow({ 1.0, 0.0, 1.0 }), MakeDatasetRow({ 1.0, 1.0, 3.0 }) })
            },
            {
                MakeFitResult({ 1.0, 2.0 }, 0.25, { 1.0, 1.0 }, { 0.25, 0.25 })
            }
        )
    };

    const auto result{ rhbm_gem::rhbm_helper::EstimateGroup(0.0, input) };

    EXPECT_EQ(rg::RHBMEstimationStatus::SINGLE_MEMBER, result.status);
    EXPECT_TRUE(result.mu_mean.isApprox(MakeVector({ 1.0, 2.0 }), 1e-12));
    EXPECT_TRUE(result.mu_prior.isApprox(MakeVector({ 1.0, 2.0 }), 1e-12));
    EXPECT_TRUE(result.beta_posterior_matrix.col(0).isApprox(MakeVector({ 1.0, 2.0 }), 1e-12));
    ASSERT_EQ(1u, result.capital_sigma_posterior_list.size());
    EXPECT_TRUE(result.capital_sigma_posterior_list[0].isApprox(Eigen::MatrixXd::Zero(2, 2), 1e-12));
}

TEST(RHBMHelperTest, EstimateGroupRejectsMissingMemberFitResults)
{
    rg::RHBMGroupEstimationInput input;
    input.basis_size = 2;
    input.member_datasets.push_back(
        { Eigen::MatrixXd::Identity(2, 2), MakeVector({ 1.0, 2.0 }) });

    EXPECT_THROW(rhbm_gem::rhbm_helper::EstimateGroup(0.0, input), std::invalid_argument);
}

TEST(RHBMHelperTest, EstimateGroupRejectsMismatchedWeightSize)
{
    rg::RHBMGroupEstimationInput input;
    input.basis_size = 2;
    input.member_datasets.push_back({
        Eigen::MatrixXd::Identity(2, 2),
        MakeVector({ 1.0, 2.0 })
    });

    auto fit_result{ MakeFitResult({ 1.0, 2.0 }, 0.25, { 1.0 }, { 0.25, 0.25 }) };
    input.member_fit_results.push_back(fit_result);

    EXPECT_THROW(rhbm_gem::rhbm_helper::EstimateGroup(0.0, input), std::invalid_argument);
}

TEST(RHBMHelperTest, EstimateGroupRejectsBasisSizeMismatchWithMemberDatasets)
{
    rg::RHBMGroupEstimationInput input;
    input.basis_size = 3;
    input.member_datasets.push_back({
        Eigen::MatrixXd::Identity(2, 2),
        MakeVector({ 1.0, 2.0 })
    });
    input.member_fit_results.push_back(
        MakeFitResult({ 1.0, 2.0 }, 0.25, { 1.0, 1.0 }, { 0.25, 0.25 }));

    EXPECT_THROW(rhbm_gem::rhbm_helper::EstimateGroup(0.0, input), std::invalid_argument);
}

TEST(RHBMHelperTest, EstimateGroupTwoMembersPinsMuPriorToMuMDPDE)
{
    const auto input{
        rhbm_gem::rhbm_helper::BuildGroupInput(
            {
                MakeDataset({ MakeDatasetRow({ 1.0, 0.0, 1.0 }), MakeDatasetRow({ 1.0, 1.0, 3.0 }) }),
                MakeDataset({ MakeDatasetRow({ 1.0, 0.0, 2.0 }), MakeDatasetRow({ 1.0, 1.0, 4.0 }) })
            },
            {
                MakeFitResult({ 1.0, 2.0 }, 0.25, { 1.0, 1.0 }, { 0.25, 0.25 }),
                MakeFitResult({ 2.0, 2.0 }, 0.5, { 1.0, 1.0 }, { 0.5, 0.5 })
            }
        )
    };

    const auto result{ rhbm_gem::rhbm_helper::EstimateGroup(0.0, input) };

    EXPECT_EQ(rg::RHBMEstimationStatus::SUCCESS, result.status);
    EXPECT_TRUE(result.mu_prior.isApprox(result.mu_mdpde, 1e-12));
    EXPECT_EQ(2, result.beta_posterior_matrix.cols());
    ASSERT_EQ(2u, result.capital_sigma_posterior_list.size());
    EXPECT_EQ(2, result.statistical_distance_array.size());
    EXPECT_EQ(2, result.outlier_flag_array.size());
}

TEST(RHBMHelperTest, EstimateGroupAllowsConvergenceFallbackStatuses)
{
    const auto input{
        rhbm_gem::rhbm_helper::BuildGroupInput(
            {
                MakeDataset({ MakeDatasetRow({ 1.0, 0.0, 1.0 }), MakeDatasetRow({ 1.0, 1.0, 3.0 }) }),
                MakeDataset({ MakeDatasetRow({ 1.0, 0.0, 2.0 }), MakeDatasetRow({ 1.0, 1.0, 4.0 }) })
            },
            {
                MakeFitResult(
                    { 1.0, 2.0 },
                    0.25,
                    { 1.0, 1.0 },
                    { 0.25, 0.25 },
                    rg::RHBMEstimationStatus::MAX_ITERATIONS_REACHED),
                MakeFitResult(
                    { 2.0, 2.0 },
                    0.5,
                    { 1.0, 1.0 },
                    { 0.5, 0.5 },
                    rg::RHBMEstimationStatus::NUMERICAL_FALLBACK)
            }
        )
    };

    const auto result{ rhbm_gem::rhbm_helper::EstimateGroup(0.0, input) };

    EXPECT_EQ(rg::RHBMEstimationStatus::SUCCESS, result.status);
    EXPECT_EQ(2, result.beta_posterior_matrix.cols());
}

TEST(RHBMHelperTest, EstimateBetaMDPDEAlphaZeroMatchesOLS)
{
    Eigen::MatrixXd X(4, 2);
    X << 1.0, 0.0,
         1.0, 1.0,
         1.0, 2.0,
         1.0, 3.0;
    const Eigen::VectorXd y{ MakeVector({ 1.0, 2.1, 2.9, 4.2 }) };

    const Eigen::VectorXd expected_beta{
        (X.transpose() * X).inverse() * X.transpose() * y
    };

    const auto result{ rhbm_gem::rhbm_helper::EstimateBetaMDPDE(0.0, MakeDataset(X, y)) };

    EXPECT_EQ(rg::RHBMEstimationStatus::SUCCESS, result.status);
    EXPECT_TRUE(result.beta_ols.isApprox(expected_beta, 1e-12));
    EXPECT_TRUE(result.beta_mdpde.isApprox(expected_beta, 1e-12));
}

TEST(RHBMHelperTest, EstimateBetaMDPDESingleDatumReturnsInsufficientData)
{
    const Eigen::MatrixXd X{ MakeVector({ 1.0, 2.0 }).transpose() };
    const Eigen::VectorXd y{ MakeVector({ 3.0 }) };

    const auto result{ rhbm_gem::rhbm_helper::EstimateBetaMDPDE(0.2, MakeDataset(X, y)) };

    EXPECT_EQ(rg::RHBMEstimationStatus::INSUFFICIENT_DATA, result.status);
    EXPECT_TRUE(result.beta_mdpde.isApprox(Eigen::VectorXd::Zero(2), 1e-12));
    EXPECT_DOUBLE_EQ(std::numeric_limits<double>::max(), result.sigma_square);
}

TEST(RHBMHelperTest, EstimateBetaMDPDERestoresEigenThreadCountAfterEarlyReturn)
{
    const int original_thread_count{ Eigen::nbThreads() };
    const int requested_thread_count{ AlternateThreadCount(original_thread_count) };
    const Eigen::MatrixXd X{ MakeVector({ 1.0, 2.0 }).transpose() };
    const Eigen::VectorXd y{ MakeVector({ 3.0 }) };
    rg::RHBMExecutionOptions options;
    options.thread_size = requested_thread_count;

    const auto result{ rhbm_gem::rhbm_helper::EstimateBetaMDPDE(0.2, MakeDataset(X, y), options) };

    EXPECT_EQ(rg::RHBMEstimationStatus::INSUFFICIENT_DATA, result.status);
    EXPECT_EQ(original_thread_count, Eigen::nbThreads());
    Eigen::setNbThreads(original_thread_count);
}

TEST(RHBMHelperTest, EstimateBetaMDPDERestoresEigenThreadCountAfterNonPositiveRequest)
{
    const int original_thread_count{ Eigen::nbThreads() };
    constexpr int test_thread_count{ 2 };
    Eigen::setNbThreads(test_thread_count);
    const Eigen::MatrixXd X{ MakeVector({ 1.0, 2.0 }).transpose() };
    const Eigen::VectorXd y{ MakeVector({ 3.0 }) };
    rg::RHBMExecutionOptions options;
    options.thread_size = 0;

    const auto result{ rhbm_gem::rhbm_helper::EstimateBetaMDPDE(0.2, MakeDataset(X, y), options) };

    EXPECT_EQ(rg::RHBMEstimationStatus::INSUFFICIENT_DATA, result.status);
    EXPECT_EQ(test_thread_count, Eigen::nbThreads());
    Eigen::setNbThreads(original_thread_count);
}

TEST(RHBMHelperTest, EstimateBetaMDPDERejectsInvalidParameters)
{
    const Eigen::MatrixXd X{ Eigen::MatrixXd::Identity(2, 2) };
    const Eigen::VectorXd y{ MakeVector({ 1.0, 2.0 }) };
    rg::RHBMExecutionOptions options;

    EXPECT_THROW(
        rhbm_gem::rhbm_helper::EstimateBetaMDPDE(-0.1, MakeDataset(X, y), options),
        std::invalid_argument
    );

    options.max_iterations = 0;
    EXPECT_THROW(
        rhbm_gem::rhbm_helper::EstimateBetaMDPDE(0.0, MakeDataset(X, y), options),
        std::invalid_argument
    );

    options = rg::RHBMExecutionOptions{};
    options.tolerance = -1.0;
    EXPECT_THROW(
        rhbm_gem::rhbm_helper::EstimateBetaMDPDE(0.0, MakeDataset(X, y), options),
        std::invalid_argument
    );

    options = rg::RHBMExecutionOptions{};
    options.data_weight_min = 0.0;
    EXPECT_THROW(
        rhbm_gem::rhbm_helper::EstimateBetaMDPDE(0.0, MakeDataset(X, y), options),
        std::invalid_argument
    );
}

TEST(RHBMHelperTest, EstimateBetaMDPDERejectsInvalidDatasetShapeAndValues)
{
    Eigen::MatrixXd X(2, 2);
    X << 1.0, 0.0,
         1.0, 1.0;

    EXPECT_THROW(
        rhbm_gem::rhbm_helper::EstimateBetaMDPDE(0.0, MakeDataset(X, MakeVector({ 1.0 }))),
        std::invalid_argument
    );

    Eigen::MatrixXd X_with_nan{ X };
    X_with_nan(0, 0) = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(
        rhbm_gem::rhbm_helper::EstimateBetaMDPDE(0.0, MakeDataset(X_with_nan, MakeVector({ 1.0, 2.0 }))),
        std::invalid_argument
    );
}

TEST(RHBMHelperTest, EstimateBetaMDPDEDatasetMessageUsesHelperContext)
{
    Eigen::MatrixXd X(2, 2);
    X << 1.0, 0.0,
         1.0, 1.0;

    ExpectInvalidArgumentMessage(
        [&]()
        {
            rhbm_gem::rhbm_helper::EstimateBetaMDPDE(
                0.0,
                MakeDataset(X, MakeVector({ 1.0 })));
        },
        "EstimateBetaMDPDE dataset is invalid. Details: X must have 1 rows."
    );
}

TEST(RHBMHelperTest, EstimateMuMDPDESingleMemberReturnsSingleMember)
{
    Eigen::MatrixXd beta_array(2, 1);
    beta_array << 1.0,
                  2.0;

    const auto result{ rhbm_gem::rhbm_helper::EstimateMuMDPDE(0.2, beta_array) };

    EXPECT_EQ(rg::RHBMEstimationStatus::SINGLE_MEMBER, result.status);
    EXPECT_TRUE(result.mu_mean.isApprox(MakeVector({ 1.0, 2.0 }), 1e-12));
    EXPECT_TRUE(result.mu_mdpde.isApprox(MakeVector({ 1.0, 2.0 }), 1e-12));
}

TEST(RHBMHelperTest, EstimateMuMDPDERejectsInvalidParameters)
{
    Eigen::MatrixXd beta_array(2, 2);
    beta_array << 1.0, 2.0,
                  3.0, 4.0;
    rg::RHBMExecutionOptions options;

    EXPECT_THROW(
        rhbm_gem::rhbm_helper::EstimateMuMDPDE(std::numeric_limits<double>::quiet_NaN(), beta_array, options),
        std::invalid_argument
    );

    options.max_iterations = 0;
    EXPECT_THROW(rhbm_gem::rhbm_helper::EstimateMuMDPDE(0.0, beta_array, options), std::invalid_argument);

    options = rg::RHBMExecutionOptions{};
    options.member_weight_min = 0.0;
    EXPECT_THROW(rhbm_gem::rhbm_helper::EstimateMuMDPDE(0.0, beta_array, options), std::invalid_argument);
}

TEST(RHBMHelperTest, EstimateWEBRejectsInvalidMemberCovarianceShapeWithHelperMessage)
{
    const std::vector<rg::RHBMMemberDataset> member_datasets{
        { Eigen::MatrixXd::Identity(2, 2), MakeVector({ 1.0, 2.0 }) },
        { Eigen::MatrixXd::Identity(2, 2), MakeVector({ 3.0, 4.0 }) }
    };

    ExpectInvalidArgumentMessage(
        [&]()
        {
            rhbm_gem::rhbm_helper::EstimateWEB(
                member_datasets,
                { MakeDiagonal({ 1.0, 1.0 }), MakeDiagonal({ 1.0, 1.0 }) },
                MakeVector({ 1.0, 2.0 }),
                { Eigen::MatrixXd::Identity(1, 1), Eigen::MatrixXd::Identity(2, 2) }
            );
        },
        "EstimateWEB member covariance input is invalid. Details: member_capital_lambda must have shape 2x2."
    );
}

TEST(RHBMHelperTest, EstimateWEBReturnsSingleMemberStatus)
{
    const std::vector<rg::RHBMMemberDataset> member_datasets{
        { Eigen::MatrixXd::Identity(2, 2), MakeVector({ 1.0, 2.0 }) }
    };

    const auto result{
        rhbm_gem::rhbm_helper::EstimateWEB(
            member_datasets,
            { MakeDiagonal({ 1.0, 1.0 }) },
            MakeVector({ 3.0, 4.0 }),
            { Eigen::MatrixXd::Identity(2, 2) }
        )
    };

    EXPECT_EQ(rg::RHBMEstimationStatus::SINGLE_MEMBER, result.status);
}

TEST(RHBMHelperTest, EstimateWEBRestoresEigenThreadCountAfterValidationException)
{
    const int original_thread_count{ Eigen::nbThreads() };
    const int requested_thread_count{ AlternateThreadCount(original_thread_count) };
    const std::vector<rg::RHBMMemberDataset> member_datasets{
        { Eigen::MatrixXd::Identity(2, 2), MakeVector({ 1.0, 2.0 }) },
        { Eigen::MatrixXd::Identity(2, 2), MakeVector({ 3.0, 4.0 }) }
    };
    rg::RHBMExecutionOptions options;
    options.thread_size = requested_thread_count;

    EXPECT_THROW(
        rhbm_gem::rhbm_helper::EstimateWEB(
            member_datasets,
            { MakeDiagonal({ 1.0, 1.0 }), MakeDiagonal({ 1.0, 1.0 }) },
            MakeVector({ 1.0, 2.0 }),
            { Eigen::MatrixXd::Identity(1, 1), Eigen::MatrixXd::Identity(2, 2) },
            options
        ),
        std::invalid_argument
    );
    EXPECT_EQ(original_thread_count, Eigen::nbThreads());
    Eigen::setNbThreads(original_thread_count);
}

TEST(RHBMHelperTest, EstimateWEBForTwoMembersPinsMuPriorToMuMDPDE)
{
    const std::vector<rg::RHBMMemberDataset> member_datasets{
        { Eigen::MatrixXd::Identity(2, 2), MakeVector({ 1.0, 3.0 }) },
        { Eigen::MatrixXd::Identity(2, 2), MakeVector({ 5.0, 7.0 }) }
    };
    const Eigen::VectorXd mu_mdpde{ MakeVector({ 10.0, 20.0 }) };

    const auto result{
        rhbm_gem::rhbm_helper::EstimateWEB(
            member_datasets,
            { MakeDiagonal({ 1.0, 1.0 }), MakeDiagonal({ 1.0, 1.0 }) },
            mu_mdpde,
            { Eigen::MatrixXd::Identity(2, 2), Eigen::MatrixXd::Identity(2, 2) }
        )
    };

    EXPECT_EQ(rg::RHBMEstimationStatus::SUCCESS, result.status);
    EXPECT_TRUE(result.mu_prior.isApprox(mu_mdpde, 1e-12));
}

TEST(RHBMHelperTest, EstimateWEBRejectsMismatchedMemberCounts)
{
    const std::vector<rg::RHBMMemberDataset> member_datasets{
        { Eigen::MatrixXd::Identity(2, 2), MakeVector({ 1.0, 2.0 }) },
        { Eigen::MatrixXd::Identity(2, 2), MakeVector({ 3.0, 4.0 }) }
    };

    EXPECT_THROW(
        rhbm_gem::rhbm_helper::EstimateWEB(
            member_datasets,
            { MakeDiagonal({ 1.0, 1.0 }) },
            MakeVector({ 1.0, 2.0 }),
            { Eigen::MatrixXd::Identity(2, 2), Eigen::MatrixXd::Identity(2, 2) }
        ),
        std::invalid_argument
    );
}

TEST(RHBMHelperTest, CalculateMemberStatisticalDistanceRejectsInvalidCovarianceShapeWithHelperMessage)
{
    ExpectInvalidArgumentMessage(
        [&]()
        {
            rhbm_gem::rhbm_helper::CalculateMemberStatisticalDistance(
                MakeVector({ 1.0, 2.0 }),
                Eigen::MatrixXd::Identity(1, 1),
                Eigen::MatrixXd::Identity(2, 2)
            );
        },
        "Statistical distance input is invalid. Details: capital_lambda must have shape 2x2."
    );
}

TEST(RHBMHelperTest, CalculateOutlierMemberFlagMatchesThresholdRule)
{
    Eigen::ArrayXd distances(3);
    distances << 1.0, 12.0, 100.0;

    const auto flags{ rhbm_gem::rhbm_helper::CalculateOutlierMemberFlag(3, distances) };

    EXPECT_FALSE(flags(0));
    EXPECT_TRUE(flags(1));
    EXPECT_TRUE(flags(2));
}
