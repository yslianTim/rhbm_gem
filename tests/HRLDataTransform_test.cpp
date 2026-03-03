#include <gtest/gtest.h>

#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <vector>

#include "HRLDataTransform.hpp"

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
} // namespace

TEST(HRLDataTransformTest, BuildMemberDatasetSplitsPredictorsAndResponse)
{
    const std::vector<Eigen::VectorXd> data_vector{
        MakeVector({ 1.0, 2.0, 10.0 }),
        MakeVector({ 3.0, 4.0, 20.0 })
    };

    const auto dataset{ HRLDataTransform::BuildMemberDataset(data_vector) };

    Eigen::MatrixXd expected_X(2, 2);
    expected_X << 1.0, 2.0,
                  3.0, 4.0;
    EXPECT_TRUE(dataset.X.isApprox(expected_X, 1e-12));
    EXPECT_TRUE(dataset.y.isApprox(MakeVector({ 10.0, 20.0 }), 1e-12));
}

TEST(HRLDataTransformTest, BuildMemberDatasetRejectsEmptyInput)
{
    const std::vector<Eigen::VectorXd> data_vector;
    EXPECT_THROW(HRLDataTransform::BuildMemberDataset(data_vector), std::invalid_argument);
}

TEST(HRLDataTransformTest, BuildMemberDatasetRejectsNonFiniteValues)
{
    const std::vector<Eigen::VectorXd> data_vector{
        MakeVector({ 1.0, std::numeric_limits<double>::quiet_NaN(), 10.0 })
    };
    EXPECT_THROW(HRLDataTransform::BuildMemberDataset(data_vector), std::invalid_argument);
}

TEST(HRLDataTransformTest, BuildMemberDatasetRejectsInconsistentBasis)
{
    const std::vector<Eigen::VectorXd> data_vector{
        MakeVector({ 1.0, 2.0, 10.0 }),
        MakeVector({ 3.0, 4.0, 5.0, 20.0 })
    };
    EXPECT_THROW(HRLDataTransform::BuildMemberDataset(data_vector), std::invalid_argument);
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

TEST(HRLDataTransformTest, BuildGroupInputBuildsStructuredRequest)
{
    const auto input{
        HRLDataTransform::BuildGroupInput(
            2,
            {
                { MakeVector({ 1.0, 0.0, 1.0 }), MakeVector({ 1.0, 1.0, 3.0 }) },
                { MakeVector({ 1.0, 0.0, 2.0 }), MakeVector({ 1.0, 1.0, 4.0 }) }
            },
            { MakeVector({ 1.0, 2.0 }), MakeVector({ 2.0, 2.0 }) },
            { 0.25, 0.5 },
            { MakeDiagonal({ 1.0, 1.0 }), MakeDiagonal({ 1.0, 1.0 }) },
            { MakeDiagonal({ 0.25, 0.25 }), MakeDiagonal({ 0.5, 0.5 }) }
        )
    };

    ASSERT_EQ(2u, input.member_datasets.size());
    ASSERT_EQ(2u, input.member_estimates.size());
    EXPECT_EQ(2, input.basis_size);
    EXPECT_TRUE(input.member_estimates[0].beta_mdpde.isApprox(MakeVector({ 1.0, 2.0 }), 1e-12));
}

TEST(HRLDataTransformTest, BuildGroupInputRejectsMismatchedMemberCounts)
{
    EXPECT_THROW(
        HRLDataTransform::BuildGroupInput(
            2,
            { { MakeVector({ 1.0, 0.0, 1.0 }) } },
            { MakeVector({ 1.0, 2.0 }), MakeVector({ 3.0, 4.0 }) },
            { 0.25 },
            { MakeDiagonal({ 1.0 }) },
            { MakeDiagonal({ 1.0 }) }
        ),
        std::invalid_argument
    );
}

TEST(HRLDataTransformTest, BuildGroupInputRejectsInconsistentWeightSize)
{
    EXPECT_THROW(
        HRLDataTransform::BuildGroupInput(
            2,
            {
                { MakeVector({ 1.0, 0.0, 1.0 }), MakeVector({ 1.0, 1.0, 3.0 }) }
            },
            { MakeVector({ 1.0, 2.0 }) },
            { 0.25 },
            { MakeDiagonal({ 1.0 }) },
            { MakeDiagonal({ 0.25, 0.25 }) }
        ),
        std::invalid_argument
    );
}
