#include <gtest/gtest.h>

#include <initializer_list>
#include <limits>
#include <stdexcept>

#include <rhbm_gem/utils/hrl/HRLModelAlgorithms.hpp>

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

HRLMemberDataset MakeDataset(
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & y)
{
    return HRLMemberDataset{ X, y };
}
} // namespace

TEST(HRLModelAlgorithmsTest, EstimateBetaMDPDEAlphaZeroMatchesOLS)
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

    const auto result{ HRLModelAlgorithms::EstimateBetaMDPDE(0.0, MakeDataset(X, y)) };

    EXPECT_EQ(HRLEstimationStatus::SUCCESS, result.status);
    EXPECT_TRUE(result.beta_ols.isApprox(expected_beta, 1e-12));
    EXPECT_TRUE(result.beta_mdpde.isApprox(expected_beta, 1e-12));
}

TEST(HRLModelAlgorithmsTest, EstimateBetaMDPDESingleDatumReturnsInsufficientData)
{
    const Eigen::MatrixXd X{ MakeVector({ 1.0, 2.0 }).transpose() };
    const Eigen::VectorXd y{ MakeVector({ 3.0 }) };

    const auto result{ HRLModelAlgorithms::EstimateBetaMDPDE(0.2, MakeDataset(X, y)) };

    EXPECT_EQ(HRLEstimationStatus::INSUFFICIENT_DATA, result.status);
    EXPECT_TRUE(result.beta_mdpde.isApprox(Eigen::VectorXd::Zero(2), 1e-12));
    EXPECT_DOUBLE_EQ(std::numeric_limits<double>::max(), result.sigma_square);
}

TEST(HRLModelAlgorithmsTest, EstimateBetaMDPDERejectsInvalidParameters)
{
    const Eigen::MatrixXd X{ Eigen::MatrixXd::Identity(2, 2) };
    const Eigen::VectorXd y{ MakeVector({ 1.0, 2.0 }) };
    HRLExecutionOptions options;

    EXPECT_THROW(
        HRLModelAlgorithms::EstimateBetaMDPDE(-0.1, MakeDataset(X, y), options),
        std::invalid_argument
    );

    options.max_iterations = 0;
    EXPECT_THROW(
        HRLModelAlgorithms::EstimateBetaMDPDE(0.0, MakeDataset(X, y), options),
        std::invalid_argument
    );

    options = HRLExecutionOptions{};
    options.tolerance = -1.0;
    EXPECT_THROW(
        HRLModelAlgorithms::EstimateBetaMDPDE(0.0, MakeDataset(X, y), options),
        std::invalid_argument
    );
}

TEST(HRLModelAlgorithmsTest, EstimateBetaMDPDERejectsInvalidDatasetShapeAndValues)
{
    Eigen::MatrixXd X(2, 2);
    X << 1.0, 0.0,
         1.0, 1.0;

    EXPECT_THROW(
        HRLModelAlgorithms::EstimateBetaMDPDE(0.0, MakeDataset(X, MakeVector({ 1.0 }))),
        std::invalid_argument
    );

    Eigen::MatrixXd X_with_nan{ X };
    X_with_nan(0, 0) = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(
        HRLModelAlgorithms::EstimateBetaMDPDE(0.0, MakeDataset(X_with_nan, MakeVector({ 1.0, 2.0 }))),
        std::invalid_argument
    );
}

TEST(HRLModelAlgorithmsTest, EstimateMuMDPDESingleMemberReturnsSingleMember)
{
    Eigen::MatrixXd beta_array(2, 1);
    beta_array << 1.0,
                  2.0;

    const auto result{ HRLModelAlgorithms::EstimateMuMDPDE(0.2, beta_array) };

    EXPECT_EQ(HRLEstimationStatus::SINGLE_MEMBER, result.status);
    EXPECT_TRUE(result.mu_mean.isApprox(MakeVector({ 1.0, 2.0 }), 1e-12));
    EXPECT_TRUE(result.mu_mdpde.isApprox(MakeVector({ 1.0, 2.0 }), 1e-12));
}

TEST(HRLModelAlgorithmsTest, EstimateMuMDPDERejectsInvalidParameters)
{
    Eigen::MatrixXd beta_array(2, 2);
    beta_array << 1.0, 2.0,
                  3.0, 4.0;
    HRLExecutionOptions options;

    EXPECT_THROW(
        HRLModelAlgorithms::EstimateMuMDPDE(std::numeric_limits<double>::quiet_NaN(), beta_array, options),
        std::invalid_argument
    );

    options.max_iterations = 0;
    EXPECT_THROW(HRLModelAlgorithms::EstimateMuMDPDE(0.0, beta_array, options), std::invalid_argument);

    options = HRLExecutionOptions{};
    options.member_weight_min = 0.0;
    EXPECT_THROW(HRLModelAlgorithms::EstimateMuMDPDE(0.0, beta_array, options), std::invalid_argument);
}

TEST(HRLModelAlgorithmsTest, EstimateWEBReturnsSingleMemberStatus)
{
    const std::vector<HRLMemberDataset> member_datasets{
        { Eigen::MatrixXd::Identity(2, 2), MakeVector({ 1.0, 2.0 }) }
    };

    const auto result{
        HRLModelAlgorithms::EstimateWEB(
            member_datasets,
            { MakeDiagonal({ 1.0, 1.0 }) },
            MakeVector({ 3.0, 4.0 }),
            { Eigen::MatrixXd::Identity(2, 2) }
        )
    };

    EXPECT_EQ(HRLEstimationStatus::SINGLE_MEMBER, result.status);
}

TEST(HRLModelAlgorithmsTest, EstimateWEBForTwoMembersPinsMuPriorToMuMDPDE)
{
    const std::vector<HRLMemberDataset> member_datasets{
        { Eigen::MatrixXd::Identity(2, 2), MakeVector({ 1.0, 3.0 }) },
        { Eigen::MatrixXd::Identity(2, 2), MakeVector({ 5.0, 7.0 }) }
    };
    const Eigen::VectorXd mu_mdpde{ MakeVector({ 10.0, 20.0 }) };

    const auto result{
        HRLModelAlgorithms::EstimateWEB(
            member_datasets,
            { MakeDiagonal({ 1.0, 1.0 }), MakeDiagonal({ 1.0, 1.0 }) },
            mu_mdpde,
            { Eigen::MatrixXd::Identity(2, 2), Eigen::MatrixXd::Identity(2, 2) }
        )
    };

    EXPECT_EQ(HRLEstimationStatus::SUCCESS, result.status);
    EXPECT_TRUE(result.mu_prior.isApprox(mu_mdpde, 1e-12));
}

TEST(HRLModelAlgorithmsTest, EstimateWEBRejectsMismatchedMemberCounts)
{
    const std::vector<HRLMemberDataset> member_datasets{
        { Eigen::MatrixXd::Identity(2, 2), MakeVector({ 1.0, 2.0 }) },
        { Eigen::MatrixXd::Identity(2, 2), MakeVector({ 3.0, 4.0 }) }
    };

    EXPECT_THROW(
        HRLModelAlgorithms::EstimateWEB(
            member_datasets,
            { MakeDiagonal({ 1.0, 1.0 }) },
            MakeVector({ 1.0, 2.0 }),
            { Eigen::MatrixXd::Identity(2, 2), Eigen::MatrixXd::Identity(2, 2) }
        ),
        std::invalid_argument
    );
}

TEST(HRLModelAlgorithmsTest, CalculateOutlierMemberFlagMatchesThresholdRule)
{
    Eigen::ArrayXd distances(3);
    distances << 1.0, 12.0, 100.0;

    const auto flags{ HRLModelAlgorithms::CalculateOutlierMemberFlag(3, distances) };

    EXPECT_FALSE(flags(0));
    EXPECT_TRUE(flags(1));
    EXPECT_TRUE(flags(2));
}
