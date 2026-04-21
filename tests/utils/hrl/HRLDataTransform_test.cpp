#include <gtest/gtest.h>

#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include <rhbm_gem/utils/hrl/HRLDataTransform.hpp>

namespace
{
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

HRLDiagonalMatrix MakeDiagonal(std::initializer_list<double> values)
{
    HRLDiagonalMatrix result(static_cast<Eigen::Index>(values.size()));
    result.diagonal() = MakeVector(values);
    return result;
}

SeriesPoint MakeSeriesPoint(std::initializer_list<double> values)
{
    std::vector<double> row(values);
    const auto response{ row.back() };
    row.pop_back();
    return SeriesPoint{ std::move(row), response };
}

HRLMemberDataset MakeDataset(std::initializer_list<SeriesPoint> rows)
{
    return HRLDataTransform::BuildMemberDataset(SeriesPointList(rows));
}

HRLBetaEstimateResult MakeFitResult(
    std::initializer_list<double> beta_values,
    double sigma_square,
    std::initializer_list<double> weight_values,
    std::initializer_list<double> covariance_values,
    HRLEstimationStatus status = HRLEstimationStatus::SUCCESS)
{
    HRLBetaEstimateResult fit_result;
    fit_result.status = status;
    fit_result.beta_mdpde = MakeVector(beta_values);
    fit_result.sigma_square = sigma_square;
    fit_result.data_weight = MakeDiagonal(weight_values);
    fit_result.data_covariance = MakeDiagonal(covariance_values);
    return fit_result;
}
} // namespace

TEST(HRLDataTransformTest, BuildMemberDatasetSplitsPredictorsAndResponse)
{
    const SeriesPointList data_series{
        SeriesPoint({ 1.0, 2.0 }, 10.0, 0.5),
        SeriesPoint({ 3.0, 4.0 }, 20.0, 2.0)
    };

    const auto dataset{ HRLDataTransform::BuildMemberDataset(data_series) };

    Eigen::MatrixXd expected_X(2, 2);
    expected_X << 1.0, 2.0,
                  3.0, 4.0;
    EXPECT_TRUE(dataset.X.isApprox(expected_X, 1e-12));
    EXPECT_TRUE(dataset.y.isApprox(MakeVector({ 10.0, 20.0 }), 1e-12));
    EXPECT_TRUE(dataset.score.isApprox(MakeVector({ 0.5, 2.0 }), 1e-12));
}

TEST(HRLDataTransformTest, BuildMemberDatasetRejectsEmptyInput)
{
    const SeriesPointList data_series;
    EXPECT_THROW(HRLDataTransform::BuildMemberDataset(data_series), std::invalid_argument);
}

TEST(HRLDataTransformTest, BuildMemberDatasetRejectsNonFiniteValues)
{
    const SeriesPointList data_series{
        SeriesPoint({ 1.0, std::numeric_limits<double>::quiet_NaN() }, 10.0)
    };
    EXPECT_THROW(HRLDataTransform::BuildMemberDataset(data_series), std::invalid_argument);
    EXPECT_THROW(
        HRLDataTransform::BuildMemberDataset({
            SeriesPoint({ 1.0, 2.0 }, std::numeric_limits<double>::quiet_NaN())
        }),
        std::invalid_argument);
    EXPECT_THROW(
        HRLDataTransform::BuildMemberDataset({
            SeriesPoint({ 1.0, 2.0 }, 10.0, std::numeric_limits<double>::quiet_NaN())
        }),
        std::invalid_argument);
}

TEST(HRLDataTransformTest, BuildMemberDatasetRejectsInconsistentBasis)
{
    const SeriesPointList data_series{
        SeriesPoint({ 1.0, 2.0 }, 10.0),
        SeriesPoint({ 3.0, 4.0, 5.0 }, 20.0)
    };
    EXPECT_THROW(HRLDataTransform::BuildMemberDataset(data_series), std::invalid_argument);
}

TEST(HRLDataTransformTest, BuildMemberDatasetRejectsEmptyBasis)
{
    const SeriesPointList data_series{
        SeriesPoint(std::vector<double>{}, 10.0)
    };
    EXPECT_THROW(HRLDataTransform::BuildMemberDataset(data_series), std::invalid_argument);
}

TEST(HRLDataTransformTest, BuildBetaMatrixMapsVectorsToColumns)
{
    const std::vector<Eigen::VectorXd> beta_list{
        MakeVector({ 1.0, 2.0 }),
        MakeVector({ 3.0, 4.0 }),
        MakeVector({ 5.0, 6.0 })
    };

    const auto beta_matrix{ HRLDataTransform::BuildBetaMatrix(beta_list) };

    Eigen::MatrixXd expected(2, 3);
    expected << 1.0, 3.0, 5.0,
                2.0, 4.0, 6.0;
    EXPECT_TRUE(beta_matrix.isApprox(expected, 1e-12));
}

TEST(HRLDataTransformTest, BuildBetaMatrixRejectsEmptyInput)
{
    const std::vector<Eigen::VectorXd> beta_list;
    EXPECT_THROW(HRLDataTransform::BuildBetaMatrix(beta_list), std::invalid_argument);
}

TEST(HRLDataTransformTest, BuildBetaMatrixRejectsEmptyBasis)
{
    const std::vector<Eigen::VectorXd> beta_list{
        Eigen::VectorXd{}
    };
    EXPECT_THROW(HRLDataTransform::BuildBetaMatrix(beta_list), std::invalid_argument);
}

TEST(HRLDataTransformTest, BuildBetaMatrixRejectsNonFiniteValue)
{
    const std::vector<Eigen::VectorXd> beta_list{
        MakeVector({ 1.0, std::numeric_limits<double>::quiet_NaN() })
    };

    EXPECT_THROW(HRLDataTransform::BuildBetaMatrix(beta_list), std::invalid_argument);
}

TEST(HRLDataTransformTest, BuildBetaMatrixRejectsInconsistentBasisSize)
{
    const std::vector<Eigen::VectorXd> beta_list{
        MakeVector({ 1.0, 2.0 }),
        MakeVector({ 3.0, 4.0, 5.0 })
    };

    EXPECT_THROW(HRLDataTransform::BuildBetaMatrix(beta_list), std::invalid_argument);
}

TEST(HRLDataTransformTest, BuildGroupInputBuildsStructuredRequest)
{
    const auto input{
        HRLDataTransform::BuildGroupInput(
            {
                MakeDataset({ MakeSeriesPoint({ 1.0, 0.0, 1.0 }), MakeSeriesPoint({ 1.0, 1.0, 3.0 }) }),
                MakeDataset({ MakeSeriesPoint({ 1.0, 0.0, 2.0 }), MakeSeriesPoint({ 1.0, 1.0, 4.0 }) })
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

TEST(HRLDataTransformTest, BuildGroupInputRejectsMismatchedMemberCounts)
{
    EXPECT_THROW(
        HRLDataTransform::BuildGroupInput(
            { MakeDataset({ MakeSeriesPoint({ 1.0, 0.0, 1.0 }) }) },
            {
                MakeFitResult({ 1.0, 2.0 }, 0.25, { 1.0 }, { 1.0 }),
                MakeFitResult({ 3.0, 4.0 }, 0.5, { 1.0 }, { 1.0 })
            }
        ),
        std::invalid_argument
    );
}

TEST(HRLDataTransformTest, BuildGroupInputRejectsInconsistentWeightSize)
{
    EXPECT_THROW(
        HRLDataTransform::BuildGroupInput(
            {
                MakeDataset({ MakeSeriesPoint({ 1.0, 0.0, 1.0 }), MakeSeriesPoint({ 1.0, 1.0, 3.0 }) })
            },
            {
                MakeFitResult({ 1.0, 2.0 }, 0.25, { 1.0 }, { 0.25, 0.25 })
            }
        ),
        std::invalid_argument
    );
}

TEST(HRLDataTransformTest, BuildGroupInputRejectsInconsistentScoreSize)
{
    auto dataset{
        MakeDataset({
            MakeSeriesPoint({ 1.0, 0.0, 1.0 }),
            MakeSeriesPoint({ 1.0, 1.0, 3.0 })
        })
    };
    dataset.score = MakeVector({ 1.0 });

    EXPECT_THROW(
        HRLDataTransform::BuildGroupInput(
            { dataset },
            {
                MakeFitResult({ 1.0, 2.0 }, 0.25, { 1.0, 1.0 }, { 0.25, 0.25 })
            }
        ),
        std::invalid_argument
    );
}

TEST(HRLDataTransformTest, BuildGroupInputRejectsInconsistentMemberBetaBasisSize)
{
    EXPECT_THROW(
        HRLDataTransform::BuildGroupInput(
            {
                MakeDataset({
                    MakeSeriesPoint({ 1.0, 0.0, 1.0 }),
                    MakeSeriesPoint({ 1.0, 1.0, 3.0 })
                })
            },
            {
                MakeFitResult({ 1.0 }, 0.25, { 1.0, 1.0 }, { 0.25, 0.25 })
            }
        ),
        std::invalid_argument
    );
}
