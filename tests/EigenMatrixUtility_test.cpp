#include <gtest/gtest.h>

#include "EigenMatrixUtility.hpp"

TEST(EigenMatrixUtilityTest, GetMedianOddNumberOfElements)
{
    Eigen::Matrix<double, 3, 3> matrix;
    matrix << 1, 2, 3,
              4, 5, 6,
              7, 8, 9;
    const double median{ EigenMatrixUtility::GetMedian(matrix) };
    EXPECT_DOUBLE_EQ(median, 5.0);
}

TEST(EigenMatrixUtilityTest, GetMedianEvenNumberOfElementsConst)
{
    const Eigen::Matrix<double, 2, 2> matrix{
        (Eigen::Matrix<double, 2, 2>() << 1, 2, 3, 4).finished()
    };
    const double median{ EigenMatrixUtility::GetMedian(matrix) };
    EXPECT_DOUBLE_EQ(median, 2.5);
}

TEST(EigenMatrixUtilityTest, GetMedianEmptyMatrixWarnsAndReturnsZero)
{
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix{ 0, 0 };
    testing::internal::CaptureStderr();
    const double median{ EigenMatrixUtility::GetMedian(matrix) };
    const std::string output{ testing::internal::GetCapturedStderr() };
    EXPECT_DOUBLE_EQ(median, 0.0);
    EXPECT_NE(
        output.find("[Warning]  [EigenMatrixUtility::GetMedian] The input data size is zero, return 0..."),
        std::string::npos);
}

TEST(EigenMatrixUtilityTest, GetMedianSortsNonConstMatrix)
{
    Eigen::Matrix<double, 1, 5> matrix;
    matrix << 5.0, 3.0, 1.0, 4.0, 2.0;
    const double median{ EigenMatrixUtility::GetMedian(matrix) };
    EXPECT_DOUBLE_EQ(median, 3.0);
    Eigen::Matrix<double, 1, 5> expected;
    expected << 1.0, 2.0, 3.0, 4.0, 5.0;
    EXPECT_TRUE(matrix.isApprox(expected));
}

TEST(EigenMatrixUtilityTest, GetMedianDoesNotModifyConstMatrix)
{
    Eigen::Matrix<double, 1, 5> original;
    original << 5.0, 3.0, 1.0, 4.0, 2.0;
    const Eigen::Matrix<double, 1, 5> matrix{ original };
    const double median{ EigenMatrixUtility::GetMedian(matrix) };
    EXPECT_DOUBLE_EQ(median, 3.0);
    EXPECT_TRUE(matrix.isApprox(original));
}

TEST(EigenMatrixUtilityTest, GetStandardDeviationMultiElement)
{
    Eigen::ArrayXd data(4);
    data << 1.0, 2.0, 3.0, 4.0;
    const double expected{ std::sqrt(5.0 / 3.0) };

    EXPECT_NEAR(expected, EigenMatrixUtility::GetStandardDeviation(data), 1e-9);

    const Eigen::ArrayXd const_data{ data };
    EXPECT_NEAR(expected, EigenMatrixUtility::GetStandardDeviation(const_data), 1e-9);
}

TEST(EigenMatrixUtilityTest, GetStandardDeviationEmptyMatrixWarnsAndReturnsZero)
{
    Eigen::ArrayXd data(0);
    testing::internal::CaptureStderr();
    const double value{ EigenMatrixUtility::GetStandardDeviation(data) };
    const std::string output{ testing::internal::GetCapturedStderr() };
    EXPECT_DOUBLE_EQ(0.0, value);
    EXPECT_NE(output.find("[Warning]  [EigenMatrixUtility::GetStandardDeviation] "
                          "The input data size is zero, return 0..."),
              std::string::npos);
    const Eigen::ArrayXd const_data(0);
    testing::internal::CaptureStderr();
    const double const_value{ EigenMatrixUtility::GetStandardDeviation(const_data) };
    const std::string const_output{ testing::internal::GetCapturedStderr() };
    EXPECT_DOUBLE_EQ(0.0, const_value);
    EXPECT_NE(const_output.find("[Warning]  [EigenMatrixUtility::GetStandardDeviation] "
                                "The input data size is zero, return 0..."),
              std::string::npos);
}

TEST(EigenMatrixUtilityTest, GetStandardDeviationSingleElement)
{
    Eigen::ArrayXd data(1);
    data << 42.0;
    EXPECT_DOUBLE_EQ(0.0, EigenMatrixUtility::GetStandardDeviation(data));

    const Eigen::ArrayXd const_data{ data };
    EXPECT_DOUBLE_EQ(0.0, EigenMatrixUtility::GetStandardDeviation(const_data));
}

TEST(EigenMatrixUtilityTest, InvertibleMatrixReturnsCorrectInverse)
{
    Eigen::Matrix2d matrix;
    matrix << 1.0, 2.0, 3.0, 4.0;
    const auto result{ EigenMatrixUtility::GetInverseMatrix(matrix) };
    const auto expected{ matrix.inverse() };
    EXPECT_TRUE(result.isApprox(expected));
}

TEST(EigenMatrixUtilityTest, SingularMatrixReturnsPseudoInverse)
{
    Eigen::Matrix2d matrix;
    matrix << 1.0, 2.0, 2.0, 4.0;
    const auto result{ EigenMatrixUtility::GetInverseMatrix(matrix) };
    const auto expected{ matrix.completeOrthogonalDecomposition().pseudoInverse() };
    EXPECT_TRUE(result.isApprox(expected));
}

TEST(EigenMatrixUtilityTest, OneByOneNonZeroMatrixReturnsInverse)
{
    Eigen::Matrix<double, 1, 1> matrix;
    matrix << 2.0;
    const auto result{ EigenMatrixUtility::GetInverseMatrix(matrix) };
    EXPECT_DOUBLE_EQ(result(0, 0), 0.5);
}

TEST(EigenMatrixUtilityTest, OneByOneZeroMatrixReturnsZero)
{
    Eigen::Matrix<double, 1, 1> matrix;
    matrix << 0.0;
    const auto result{ EigenMatrixUtility::GetInverseMatrix(matrix) };
    EXPECT_DOUBLE_EQ(result(0, 0), 0.0);
}

TEST(EigenMatrixUtilityTest, IdentityMatrixRemainsUnchangedAfterInverse)
{
    Eigen::Matrix3d matrix{ Eigen::Matrix3d::Identity() };
    const auto result{ EigenMatrixUtility::GetInverseMatrix(matrix) };
    EXPECT_TRUE(result.isApprox(Eigen::Matrix3d::Identity()));
}

TEST(EigenMatrixUtilityTest, InverseDiagonalHandlesZeros)
{
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> diag(4);
    diag.diagonal() << 2.0, 0.0, -4.0, 5.0;
    auto inv{ EigenMatrixUtility::GetInverseDiagonalMatrix(diag) };
    EXPECT_DOUBLE_EQ(inv.diagonal()[0], 0.5);
    EXPECT_DOUBLE_EQ(inv.diagonal()[1], 0.0);
    EXPECT_DOUBLE_EQ(inv.diagonal()[2], -0.25);
    EXPECT_DOUBLE_EQ(inv.diagonal()[3], 0.2);
}
