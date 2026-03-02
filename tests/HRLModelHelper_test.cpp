#include <gtest/gtest.h>

#include <cmath>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <vector>

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wkeyword-macro"
#endif
#define private public
#include "HRLModelHelper.hpp"
#undef private
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

using DMatrixXd = Eigen::DiagonalMatrix<double, Eigen::Dynamic>;

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

DMatrixXd MakeDiagonal(std::initializer_list<double> values)
{
    DMatrixXd result(static_cast<Eigen::Index>(values.size()));
    result.diagonal() = MakeVector(values);
    return result;
}
} // namespace

TEST(HRLModelHelperTest, ConstructorRejectsNonPositiveSizes)
{
    EXPECT_THROW(HRLModelHelper(0, 1), std::invalid_argument);
    EXPECT_THROW(HRLModelHelper(1, 0), std::invalid_argument);
    EXPECT_THROW(HRLModelHelper(-1, 1), std::invalid_argument);
    EXPECT_THROW(HRLModelHelper(1, -1), std::invalid_argument);
}

TEST(HRLModelHelperTest, BuildBasisVectorAndResponseArraySplitsPredictorsAndResponse)
{
    const std::vector<Eigen::VectorXd> data_vector{
        MakeVector({ 1.0, 2.0, 10.0 }),
        MakeVector({ 3.0, 4.0, 20.0 })
    };

    const auto [X, y]{ HRLModelHelper::BuildBasisVectorAndResponseArray(data_vector, true) };

    Eigen::MatrixXd expected_X(2, 2);
    expected_X << 1.0, 2.0,
                  3.0, 4.0;
    const Eigen::VectorXd expected_y{ MakeVector({ 10.0, 20.0 }) };

    EXPECT_TRUE(X.isApprox(expected_X, 1e-12));
    EXPECT_TRUE(y.isApprox(expected_y, 1e-12));
}

TEST(HRLModelHelperTest, BuildBasisVectorAndResponseArrayRejectsNonFiniteValues)
{
    const std::vector<Eigen::VectorXd> data_vector{
        MakeVector({ 1.0, std::numeric_limits<double>::quiet_NaN(), 10.0 })
    };

    EXPECT_THROW(
        HRLModelHelper::BuildBasisVectorAndResponseArray(data_vector, true),
        std::invalid_argument);
}

TEST(HRLModelHelperTest, ConvertBetaListToMatrixMapsVectorsToColumns)
{
    const std::vector<Eigen::VectorXd> beta_list{
        MakeVector({ 1.0, 2.0 }),
        MakeVector({ 3.0, 4.0 }),
        MakeVector({ 5.0, 6.0 })
    };

    const auto beta_matrix{ HRLModelHelper::ConvertBetaListToMatrix(beta_list, true) };

    Eigen::MatrixXd expected(2, 3);
    expected << 1.0, 3.0, 5.0,
                2.0, 4.0, 6.0;
    EXPECT_TRUE(beta_matrix.isApprox(expected, 1e-12));
}

TEST(HRLModelHelperTest, SetMemberDataEntriesListRejectsUnexpectedMemberSize)
{
    HRLModelHelper helper(2, 2);
    const std::vector<std::vector<Eigen::VectorXd>> data_list{
        { MakeVector({ 1.0, 2.0, 10.0 }) }
    };

    EXPECT_THROW(helper.SetMemberDataEntriesList(data_list), std::invalid_argument);
}

TEST(HRLModelHelperTest, SetMemberDataEntriesListRejectsUnexpectedBasisSize)
{
    HRLModelHelper helper(2, 2);
    const std::vector<std::vector<Eigen::VectorXd>> data_list{
        { MakeVector({ 1.0, 2.0, 3.0, 10.0 }) },
        { MakeVector({ 1.0, 2.0, 10.0 }) }
    };

    EXPECT_THROW(helper.SetMemberDataEntriesList(data_list), std::invalid_argument);
}

TEST(HRLModelHelperTest, SetMemberBetaMDPDEListRejectsUnexpectedMemberSize)
{
    HRLModelHelper helper(2, 2);
    const std::vector<Eigen::VectorXd> beta_list{
        MakeVector({ 1.0, 2.0 })
    };
    const std::vector<double> sigma_square_list{ 0.5 };
    const std::vector<DMatrixXd> W_list{ MakeDiagonal({ 1.0, 1.0 }) };
    const std::vector<DMatrixXd> capital_sigma_list{ MakeDiagonal({ 1.0, 1.0 }) };

    EXPECT_THROW(
        helper.SetMemberBetaMDPDEList(beta_list, sigma_square_list, W_list, capital_sigma_list),
        std::invalid_argument);
}

TEST(HRLModelHelperTest, SetMemberBetaMDPDEListRejectsUnexpectedBasisSize)
{
    HRLModelHelper helper(2, 2);
    const std::vector<Eigen::VectorXd> beta_list{
        MakeVector({ 1.0, 2.0, 3.0 }),
        MakeVector({ 4.0, 5.0, 6.0 })
    };
    const std::vector<double> sigma_square_list{ 0.5, 0.25 };
    const std::vector<DMatrixXd> W_list{
        MakeDiagonal({ 1.0, 1.0 }),
        MakeDiagonal({ 1.0, 1.0 })
    };
    const std::vector<DMatrixXd> capital_sigma_list{
        MakeDiagonal({ 1.0, 1.0 }),
        MakeDiagonal({ 1.0, 1.0 })
    };

    EXPECT_THROW(
        helper.SetMemberBetaMDPDEList(beta_list, sigma_square_list, W_list, capital_sigma_list),
        std::invalid_argument);
}

TEST(HRLModelHelperTest, RunAlphaRTrainingWithExactLinearDataReturnsZeroError)
{
    const std::vector<Eigen::VectorXd> data_list{
        MakeVector({ 1.0, 0.0, 1.0 }),
        MakeVector({ 1.0, 1.0, 3.0 }),
        MakeVector({ 1.0, 2.0, 5.0 }),
        MakeVector({ 1.0, 3.0, 7.0 }),
        MakeVector({ 1.0, 4.0, 9.0 }),
        MakeVector({ 1.0, 5.0, 11.0 })
    };
    const std::vector<double> alpha_list{ 0.0, 0.5 };

    const auto error_sum_list{ HRLModelHelper::RunAlphaRTraining(data_list, 3, alpha_list) };

    ASSERT_EQ(2, error_sum_list.size());
    EXPECT_NEAR(0.0, error_sum_list(0), 1e-12);
    EXPECT_NEAR(0.0, error_sum_list(1), 1e-12);
}

TEST(HRLModelHelperTest, RunAlphaGTrainingWithIdenticalBetasReturnsZeroError)
{
    const std::vector<Eigen::VectorXd> data_list(6, MakeVector({ 1.5, -0.5 }));
    const std::vector<double> alpha_list{ 0.0, 0.5, 1.0 };

    const auto error_sum_list{ HRLModelHelper::RunAlphaGTraining(data_list, 4, alpha_list) };

    ASSERT_EQ(3, error_sum_list.size());
    EXPECT_NEAR(0.0, error_sum_list(0), 1e-12);
    EXPECT_NEAR(0.0, error_sum_list(1), 1e-12);
    EXPECT_NEAR(0.0, error_sum_list(2), 1e-12);
}

TEST(HRLModelHelperTest, AlgorithmBetaMDPDEReturnsFallbackForSingleDatum)
{
    const Eigen::MatrixXd X{ MakeVector({ 1.0, 2.0 }).transpose() };
    const Eigen::VectorXd y{ MakeVector({ 3.0 }) };
    Eigen::VectorXd beta_ols;
    Eigen::VectorXd beta;
    double sigma_square{ 0.0 };
    DMatrixXd W;
    DMatrixXd capital_sigma;

    const bool success{ HRLModelHelper::AlgorithmBetaMDPDE(
        0.2, X, y, beta_ols, beta, sigma_square, W, capital_sigma, true) };

    EXPECT_FALSE(success);
    EXPECT_TRUE(beta_ols.isApprox(Eigen::VectorXd::Zero(2), 1e-12));
    EXPECT_TRUE(beta.isApprox(Eigen::VectorXd::Zero(2), 1e-12));
    EXPECT_DOUBLE_EQ(std::numeric_limits<double>::max(), sigma_square);
    EXPECT_TRUE(W.diagonal().isApprox(MakeVector({ 1.0 }), 1e-12));
    EXPECT_TRUE(capital_sigma.diagonal().isApprox(MakeVector({ 1.0 }), 1e-12));
}

TEST(HRLModelHelperTest, AlgorithmBetaMDPDEWithAlphaZeroMatchesOLS)
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
    const Eigen::VectorXd residual{ y - (X * expected_beta) };
    const double expected_sigma_square{ residual.squaredNorm() / static_cast<double>(y.size()) };

    Eigen::VectorXd beta_ols;
    Eigen::VectorXd beta;
    double sigma_square{ 0.0 };
    DMatrixXd W;
    DMatrixXd capital_sigma;

    const bool success{ HRLModelHelper::AlgorithmBetaMDPDE(
        0.0, X, y, beta_ols, beta, sigma_square, W, capital_sigma, true) };

    EXPECT_TRUE(success);
    EXPECT_TRUE(beta_ols.isApprox(expected_beta, 1e-12));
    EXPECT_TRUE(beta.isApprox(expected_beta, 1e-12));
    EXPECT_NEAR(expected_sigma_square, sigma_square, 1e-12);
    EXPECT_TRUE(W.diagonal().isApprox(Eigen::VectorXd::Ones(4), 1e-12));
    EXPECT_TRUE(capital_sigma.diagonal().isApprox(
        Eigen::VectorXd::Constant(4, expected_sigma_square), 1e-12));
}

TEST(HRLModelHelperTest, AlgorithmBetaMDPDERejectsInvalidParameters)
{
    Eigen::MatrixXd X(2, 2);
    X << 1.0, 0.0,
         1.0, 1.0;
    const Eigen::VectorXd y{ MakeVector({ 1.0, 2.0 }) };
    Eigen::VectorXd beta_ols;
    Eigen::VectorXd beta;
    double sigma_square{ 0.0 };
    DMatrixXd W;
    DMatrixXd capital_sigma;

    EXPECT_THROW(
        HRLModelHelper::AlgorithmBetaMDPDE(
            -0.1, X, y, beta_ols, beta, sigma_square, W, capital_sigma, true),
        std::invalid_argument);
    EXPECT_THROW(
        HRLModelHelper::AlgorithmBetaMDPDE(
            0.1, X, y, beta_ols, beta, sigma_square, W, capital_sigma, true, 1, 0),
        std::invalid_argument);
    EXPECT_THROW(
        HRLModelHelper::AlgorithmBetaMDPDE(
            0.1,
            X,
            y,
            beta_ols,
            beta,
            sigma_square,
            W,
            capital_sigma,
            true,
            1,
            10,
            -1.0),
        std::invalid_argument);
}

TEST(HRLModelHelperTest, AlgorithmMuMDPDEReturnsFallbackForSingleMember)
{
    Eigen::MatrixXd beta_array(2, 1);
    beta_array << 1.0,
                  2.0;
    Eigen::VectorXd mu_mean;
    Eigen::VectorXd mu;
    Eigen::ArrayXd omega_array;
    double omega_sum{ 0.0 };
    Eigen::MatrixXd capital_lambda;
    std::vector<Eigen::MatrixXd> member_capital_lambda_list;

    const bool success{ HRLModelHelper::AlgorithmMuMDPDE(
        0.2,
        beta_array,
        mu_mean,
        mu,
        omega_array,
        omega_sum,
        capital_lambda,
        member_capital_lambda_list,
        true) };

    EXPECT_FALSE(success);
    EXPECT_TRUE(mu_mean.isApprox(MakeVector({ 1.0, 2.0 }), 1e-12));
    EXPECT_TRUE(mu.isApprox(MakeVector({ 1.0, 2.0 }), 1e-12));
    EXPECT_TRUE(omega_array.isApprox(Eigen::ArrayXd::Ones(1), 1e-12));
    EXPECT_DOUBLE_EQ(1.0, omega_sum);
    EXPECT_TRUE(capital_lambda.isApprox(Eigen::MatrixXd::Identity(2, 2), 1e-12));
    ASSERT_EQ(1u, member_capital_lambda_list.size());
}

TEST(HRLModelHelperTest, AlgorithmMuMDPDEWithAlphaZeroUsesMeanAndSampleCovariance)
{
    Eigen::MatrixXd beta_array(2, 3);
    beta_array << 1.0, 3.0, 5.0,
                  2.0, 4.0, 6.0;
    const Eigen::VectorXd expected_mean{ MakeVector({ 3.0, 4.0 }) };
    Eigen::MatrixXd expected_capital_lambda(2, 2);
    expected_capital_lambda << 8.0 / 3.0, 8.0 / 3.0,
                               8.0 / 3.0, 8.0 / 3.0;

    Eigen::VectorXd mu_mean;
    Eigen::VectorXd mu;
    Eigen::ArrayXd omega_array;
    double omega_sum{ 0.0 };
    Eigen::MatrixXd capital_lambda;
    std::vector<Eigen::MatrixXd> member_capital_lambda_list;

    const bool success{ HRLModelHelper::AlgorithmMuMDPDE(
        0.0,
        beta_array,
        mu_mean,
        mu,
        omega_array,
        omega_sum,
        capital_lambda,
        member_capital_lambda_list,
        true) };

    EXPECT_TRUE(success);
    EXPECT_TRUE(mu_mean.isApprox(expected_mean, 1e-12));
    EXPECT_TRUE(mu.isApprox(expected_mean, 1e-12));
    EXPECT_TRUE(omega_array.isApprox(Eigen::ArrayXd::Ones(3), 1e-12));
    EXPECT_DOUBLE_EQ(3.0, omega_sum);
    EXPECT_TRUE(capital_lambda.isApprox(expected_capital_lambda, 1e-12));
    ASSERT_EQ(3u, member_capital_lambda_list.size());
    EXPECT_TRUE(member_capital_lambda_list[0].isApprox(expected_capital_lambda, 1e-12));
    EXPECT_TRUE(member_capital_lambda_list[1].isApprox(expected_capital_lambda, 1e-12));
    EXPECT_TRUE(member_capital_lambda_list[2].isApprox(expected_capital_lambda, 1e-12));
}

TEST(HRLModelHelperTest, AlgorithmMuMDPDERejectsInvalidParameters)
{
    Eigen::MatrixXd beta_array(2, 2);
    beta_array << 1.0, 2.0,
                  3.0, 4.0;
    Eigen::VectorXd mu_mean;
    Eigen::VectorXd mu;
    Eigen::ArrayXd omega_array;
    double omega_sum{ 0.0 };
    Eigen::MatrixXd capital_lambda;
    std::vector<Eigen::MatrixXd> member_capital_lambda_list;

    EXPECT_THROW(
        HRLModelHelper::AlgorithmMuMDPDE(
            -0.1,
            beta_array,
            mu_mean,
            mu,
            omega_array,
            omega_sum,
            capital_lambda,
            member_capital_lambda_list,
            true),
        std::invalid_argument);
    EXPECT_THROW(
        HRLModelHelper::AlgorithmMuMDPDE(
            0.1,
            beta_array,
            mu_mean,
            mu,
            omega_array,
            omega_sum,
            capital_lambda,
            member_capital_lambda_list,
            true,
            1,
            0),
        std::invalid_argument);
    EXPECT_THROW(
        HRLModelHelper::AlgorithmMuMDPDE(
            0.1,
            beta_array,
            mu_mean,
            mu,
            omega_array,
            omega_sum,
            capital_lambda,
            member_capital_lambda_list,
            true,
            1,
            10,
            -1.0),
        std::invalid_argument);
}

TEST(HRLModelHelperTest, AlgorithmWEBReturnsFalseForSingleMember)
{
    const std::vector<Eigen::MatrixXd> X_list{ Eigen::MatrixXd::Identity(2, 2) };
    const std::vector<Eigen::VectorXd> y_list{ MakeVector({ 1.0, 2.0 }) };
    const std::vector<DMatrixXd> capital_sigma_list{ MakeDiagonal({ 1.0, 1.0 }) };
    const Eigen::VectorXd mu_mdpde{ MakeVector({ 3.0, 4.0 }) };
    const std::vector<Eigen::MatrixXd> member_capital_lambda_list{
        Eigen::MatrixXd::Identity(2, 2)
    };
    Eigen::VectorXd mu_prior;
    Eigen::MatrixXd beta_posterior_array;
    std::vector<Eigen::MatrixXd> capital_sigma_posterior_list;

    const bool success{ HRLModelHelper::AlgorithmWEB(
        X_list,
        y_list,
        capital_sigma_list,
        mu_mdpde,
        member_capital_lambda_list,
        mu_prior,
        beta_posterior_array,
        capital_sigma_posterior_list,
        true) };

    EXPECT_FALSE(success);
    EXPECT_TRUE(mu_prior.isApprox(Eigen::VectorXd::Zero(2), 1e-12));
    EXPECT_TRUE(beta_posterior_array.isApprox(Eigen::MatrixXd::Zero(2, 1), 1e-12));
    EXPECT_TRUE(capital_sigma_posterior_list.empty());
}

TEST(HRLModelHelperTest, AlgorithmWEBForTwoMembersPinsMuPriorToMuMDPDE)
{
    const std::vector<Eigen::MatrixXd> X_list{
        Eigen::MatrixXd::Identity(2, 2),
        Eigen::MatrixXd::Identity(2, 2)
    };
    const std::vector<Eigen::VectorXd> y_list{
        MakeVector({ 1.0, 3.0 }),
        MakeVector({ 5.0, 7.0 })
    };
    const std::vector<DMatrixXd> capital_sigma_list{
        MakeDiagonal({ 1.0, 1.0 }),
        MakeDiagonal({ 1.0, 1.0 })
    };
    const Eigen::VectorXd mu_mdpde{ MakeVector({ 10.0, 20.0 }) };
    const std::vector<Eigen::MatrixXd> member_capital_lambda_list{
        Eigen::MatrixXd::Identity(2, 2),
        Eigen::MatrixXd::Identity(2, 2)
    };
    Eigen::VectorXd mu_prior;
    Eigen::MatrixXd beta_posterior_array;
    std::vector<Eigen::MatrixXd> capital_sigma_posterior_list;

    const bool success{ HRLModelHelper::AlgorithmWEB(
        X_list,
        y_list,
        capital_sigma_list,
        mu_mdpde,
        member_capital_lambda_list,
        mu_prior,
        beta_posterior_array,
        capital_sigma_posterior_list,
        true) };

    Eigen::MatrixXd expected_beta_posterior(2, 2);
    expected_beta_posterior << 5.5, 7.5,
                               11.5, 13.5;

    EXPECT_TRUE(success);
    EXPECT_TRUE(mu_prior.isApprox(mu_mdpde, 1e-12));
    EXPECT_TRUE(beta_posterior_array.isApprox(expected_beta_posterior, 1e-12));
    ASSERT_EQ(2u, capital_sigma_posterior_list.size());
    EXPECT_TRUE(capital_sigma_posterior_list[0].isApprox(0.5 * Eigen::MatrixXd::Identity(2, 2), 1e-12));
    EXPECT_TRUE(capital_sigma_posterior_list[1].isApprox(0.5 * Eigen::MatrixXd::Identity(2, 2), 1e-12));
}

TEST(HRLModelHelperTest, RunGroupEstimationWithSingleMemberUsesFallbackOutputs)
{
    HRLModelHelper helper(2, 1);
    const std::vector<std::vector<Eigen::VectorXd>> data_list{
        {
            MakeVector({ 1.0, 0.0, 1.0 }),
            MakeVector({ 1.0, 1.0, 3.0 })
        }
    };
    const std::vector<Eigen::VectorXd> beta_list{ MakeVector({ 1.0, 2.0 }) };
    const std::vector<double> sigma_square_list{ 0.25 };
    const std::vector<DMatrixXd> W_list{ MakeDiagonal({ 1.0, 1.0 }) };
    const std::vector<DMatrixXd> capital_sigma_list{ MakeDiagonal({ 0.25, 0.25 }) };

    helper.SetMemberDataEntriesList(data_list);
    helper.SetMemberBetaMDPDEList(beta_list, sigma_square_list, W_list, capital_sigma_list);
    helper.RunGroupEstimation(0.0, true);

    EXPECT_TRUE(helper.GetMuVectorMean().isApprox(beta_list[0], 1e-12));
    EXPECT_TRUE(helper.GetMuVectorMDPDE().isApprox(beta_list[0], 1e-12));
    EXPECT_TRUE(helper.GetMuVectorPrior().isApprox(beta_list[0], 1e-12));
    EXPECT_TRUE(helper.GetCapitalLambdaMatrix().isApprox(Eigen::MatrixXd::Identity(2, 2), 1e-12));
    EXPECT_TRUE(helper.GetBetaPosterior(0).isApprox(beta_list[0], 1e-12));
    EXPECT_TRUE(helper.GetCapitalSigmaMatrixPosterior(0).isApprox(Eigen::MatrixXd::Zero(2, 2), 1e-12));
    EXPECT_DOUBLE_EQ(0.0, helper.GetStatisticalDistance(0));
    EXPECT_FALSE(helper.GetOutlierFlag(0));
}

TEST(HRLModelHelperTest, CalculateDataWeightAlphaZeroReturnsUnitWeights)
{
    Eigen::MatrixXd X(2, 2);
    X << 1.0, 0.0,
         1.0, 1.0;
    const Eigen::VectorXd y{ MakeVector({ 2.0, 3.0 }) };
    const Eigen::VectorXd beta{ MakeVector({ 2.0, 1.0 }) };

    const auto W{ HRLModelHelper::CalculateDataWeight(0.0, X, y, beta, 0.5) };

    EXPECT_TRUE(W.diagonal().isApprox(Eigen::VectorXd::Ones(2), 1e-12));
}

TEST(HRLModelHelperTest, CalculateDataWeightUsesMinimumFloorForInvalidSigmaSquare)
{
    const Eigen::MatrixXd X{ Eigen::MatrixXd::Identity(2, 2) };
    const Eigen::VectorXd y{ MakeVector({ 10.0, 0.0 }) };
    const Eigen::VectorXd beta{ Eigen::VectorXd::Zero(2) };

    const auto W{ HRLModelHelper::CalculateDataWeight(1.0, X, y, beta, 0.0, 0.2) };

    EXPECT_NEAR(0.2, W.diagonal()(0), 1e-12);
    EXPECT_NEAR(1.0, W.diagonal()(1), 1e-12);
}

TEST(HRLModelHelperTest, CalculateDataCovarianceFallsBackToIdentityForZeroInverseTrace)
{
    const auto capital_sigma{
        HRLModelHelper::CalculateDataCovariance(2.0, MakeDiagonal({ 0.0, 0.0, 0.0 }))
    };

    EXPECT_TRUE(capital_sigma.diagonal().isApprox(Eigen::VectorXd::Ones(3), 1e-12));
}

TEST(HRLModelHelperTest, CalculateMemberWeightClipsVerySmallValuesPerMember)
{
    Eigen::MatrixXd beta_array(2, 2);
    beta_array << 0.0, 10.0,
                  0.0, 10.0;
    const Eigen::VectorXd mu{ Eigen::VectorXd::Zero(2) };
    const Eigen::MatrixXd capital_lambda{ Eigen::MatrixXd::Identity(2, 2) };

    const auto omega_array{
        HRLModelHelper::CalculateMemberWeight(1.0, beta_array, mu, capital_lambda, 0.4)
    };

    EXPECT_NEAR(1.0, omega_array(0), 1e-12);
    EXPECT_NEAR(0.2, omega_array(1), 1e-12);
}

TEST(HRLModelHelperTest, DistancesFarFromThresholdAreFlaggedConsistently)
{
    Eigen::ArrayXd distances(3);
    distances << 1.0, 12.0, 100.0;

    const auto flags{ HRLModelHelper::CalculateOutlierMemberFlag(3, distances) };

    ASSERT_EQ(3, flags.size());
    EXPECT_FALSE(flags(0));
    EXPECT_TRUE(flags(1));
    EXPECT_TRUE(flags(2));
}

TEST(HRLModelHelperTest, BasisTwoBoundaryReflectsBoostOrLegacyThreshold)
{
    Eigen::ArrayXd distances(1);
    distances << 9.2102;

    const auto flags{ HRLModelHelper::CalculateOutlierMemberFlag(2, distances) };

#ifdef HAVE_BOOST
    EXPECT_FALSE(flags(0));
#else
    EXPECT_TRUE(flags(0));
#endif
}

TEST(HRLModelHelperTest, GettersRejectOutOfRangeMemberIds)
{
    HRLModelHelper helper(2, 2);

    EXPECT_THROW(helper.GetOutlierFlag(-1), std::out_of_range);
    EXPECT_THROW(helper.GetStatisticalDistance(2), std::out_of_range);
    EXPECT_THROW(helper.GetBetaPosterior(2), std::out_of_range);
    EXPECT_THROW(helper.GetCapitalSigmaMatrixPosterior(2), std::out_of_range);
}
