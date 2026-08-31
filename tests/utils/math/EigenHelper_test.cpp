#include <gtest/gtest.h>

#include <array>
#include <type_traits>

#include <rhbm_gem/utils/math/EigenHelper.hpp>

namespace
{
using ScopedEigenThreadCount = rhbm_gem::eigen_helper::ScopedEigenThreadCount;

static_assert(!std::is_copy_constructible_v<ScopedEigenThreadCount>);
static_assert(!std::is_copy_assignable_v<ScopedEigenThreadCount>);

class EigenThreadCountTest : public testing::Test
{
protected:
    void SetUp() override
    {
        m_original_thread_count = Eigen::nbThreads();
    }

    void TearDown() override
    {
        Eigen::setNbThreads(m_original_thread_count);
    }

    int AlternateThreadCount() const
    {
        return (m_original_thread_count == 1) ? 2 : 1;
    }

    bool SupportsThreadCount(int thread_count)
    {
        Eigen::setNbThreads(thread_count);
        const bool supported{ Eigen::nbThreads() == thread_count };
        Eigen::setNbThreads(m_original_thread_count);
        return supported;
    }

    int m_original_thread_count{ 1 };
};
} // namespace

TEST_F(EigenThreadCountTest, PositiveRequestAppliesWithinScopeAndRestoresPreviousCount)
{
    const int requested_thread_count{ AlternateThreadCount() };
    if (!SupportsThreadCount(requested_thread_count))
    {
        GTEST_SKIP() << "Eigen thread-count changes are unavailable in this build";
    }

    {
        const ScopedEigenThreadCount guard{ requested_thread_count };
        EXPECT_EQ(requested_thread_count, Eigen::nbThreads());
    }

    EXPECT_EQ(m_original_thread_count, Eigen::nbThreads());
}

TEST_F(EigenThreadCountTest, NonPositiveRequestLeavesCurrentCountUnchanged)
{
    {
        const ScopedEigenThreadCount zero_guard{ 0 };
        EXPECT_EQ(m_original_thread_count, Eigen::nbThreads());
    }
    {
        const ScopedEigenThreadCount negative_guard{ -1 };
        EXPECT_EQ(m_original_thread_count, Eigen::nbThreads());
    }

    EXPECT_EQ(m_original_thread_count, Eigen::nbThreads());
}

TEST_F(EigenThreadCountTest, NestedGuardsRestoreEachPreviousCount)
{
    const int outer_thread_count{ AlternateThreadCount() };
    if (!SupportsThreadCount(outer_thread_count))
    {
        GTEST_SKIP() << "Eigen thread-count changes are unavailable in this build";
    }

    {
        const ScopedEigenThreadCount outer_guard{ outer_thread_count };
        EXPECT_EQ(outer_thread_count, Eigen::nbThreads());
        {
            const ScopedEigenThreadCount inner_guard{ m_original_thread_count };
            EXPECT_EQ(m_original_thread_count, Eigen::nbThreads());
        }
        EXPECT_EQ(outer_thread_count, Eigen::nbThreads());
    }

    EXPECT_EQ(m_original_thread_count, Eigen::nbThreads());
}

TEST(EigenHelperTest, ToEigenVectorConvertsFloatArray3D)
{
    const std::array<float, 3> value{ 1.0f, -2.5f, 3.25f };
    const auto vector{ rhbm_gem::eigen_helper::ToEigenVector(value) };

    ASSERT_EQ(vector.size(), 3);
    EXPECT_DOUBLE_EQ(vector(0), 1.0);
    EXPECT_DOUBLE_EQ(vector(1), -2.5);
    EXPECT_DOUBLE_EQ(vector(2), 3.25);
}

TEST(EigenHelperTest, ToEigenVectorConvertsDoubleArray4D)
{
    const std::array<double, 4> value{ 1.0, 2.0, -3.0, 4.5 };
    const auto vector{ rhbm_gem::eigen_helper::ToEigenVector(value) };

    ASSERT_EQ(vector.size(), 4);
    EXPECT_DOUBLE_EQ(vector(0), 1.0);
    EXPECT_DOUBLE_EQ(vector(1), 2.0);
    EXPECT_DOUBLE_EQ(vector(2), -3.0);
    EXPECT_DOUBLE_EQ(vector(3), 4.5);
}

TEST(EigenHelperTest, GetMedianOddNumberOfElements)
{
    Eigen::Matrix<double, 3, 3> matrix;
    matrix << 1, 2, 3,
              4, 5, 6,
              7, 8, 9;
    const double median{ rhbm_gem::eigen_helper::GetMedian(matrix) };
    EXPECT_DOUBLE_EQ(median, 5.0);
}

TEST(EigenHelperTest, GetMedianEvenNumberOfElementsConst)
{
    const Eigen::Matrix<double, 2, 2> matrix{
        (Eigen::Matrix<double, 2, 2>() << 1, 2, 3, 4).finished()
    };
    const double median{ rhbm_gem::eigen_helper::GetMedian(matrix) };
    EXPECT_DOUBLE_EQ(median, 2.5);
}

TEST(EigenHelperTest, GetMedianEmptyMatrixReturnsZero)
{
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix{ 0, 0 };
    const double median{ rhbm_gem::eigen_helper::GetMedian(matrix) };
    EXPECT_DOUBLE_EQ(median, 0.0);
}

TEST(EigenHelperTest, GetMedianSortsNonConstMatrix)
{
    Eigen::Matrix<double, 1, 5> matrix;
    matrix << 5.0, 3.0, 1.0, 4.0, 2.0;
    const double median{ rhbm_gem::eigen_helper::GetMedian(matrix) };
    EXPECT_DOUBLE_EQ(median, 3.0);
    Eigen::Matrix<double, 1, 5> expected;
    expected << 1.0, 2.0, 3.0, 4.0, 5.0;
    EXPECT_TRUE(matrix.isApprox(expected));
}

TEST(EigenHelperTest, GetMedianDoesNotModifyConstMatrix)
{
    Eigen::Matrix<double, 1, 5> original;
    original << 5.0, 3.0, 1.0, 4.0, 2.0;
    const Eigen::Matrix<double, 1, 5> matrix{ original };
    const double median{ rhbm_gem::eigen_helper::GetMedian(matrix) };
    EXPECT_DOUBLE_EQ(median, 3.0);
    EXPECT_TRUE(matrix.isApprox(original));
}

TEST(EigenHelperTest, GetStandardDeviationMultiElement)
{
    Eigen::ArrayXd data(4);
    data << 1.0, 2.0, 3.0, 4.0;
    const double expected{ std::sqrt(5.0 / 3.0) };

    EXPECT_NEAR(expected, rhbm_gem::eigen_helper::GetStandardDeviation(data), 1e-9);

    const Eigen::ArrayXd const_data{ data };
    EXPECT_NEAR(expected, rhbm_gem::eigen_helper::GetStandardDeviation(const_data), 1e-9);
}

TEST(EigenHelperTest, GetStandardDeviationEmptyMatrixReturnsZero)
{
    Eigen::ArrayXd data(0);
    const double value{ rhbm_gem::eigen_helper::GetStandardDeviation(data) };
    EXPECT_DOUBLE_EQ(0.0, value);
    const Eigen::ArrayXd const_data(0);
    const double const_value{ rhbm_gem::eigen_helper::GetStandardDeviation(const_data) };
    EXPECT_DOUBLE_EQ(0.0, const_value);
}

TEST(EigenHelperTest, GetStandardDeviationSingleElement)
{
    Eigen::ArrayXd data(1);
    data << 42.0;
    EXPECT_DOUBLE_EQ(0.0, rhbm_gem::eigen_helper::GetStandardDeviation(data));

    const Eigen::ArrayXd const_data{ data };
    EXPECT_DOUBLE_EQ(0.0, rhbm_gem::eigen_helper::GetStandardDeviation(const_data));
}

TEST(EigenHelperTest, GetStandardDeviationIdenticalElementsReturnsZero)
{
    Eigen::ArrayXd data(5);
    data << 7.0, 7.0, 7.0, 7.0, 7.0;
    const double value{ rhbm_gem::eigen_helper::GetStandardDeviation(data) };
    EXPECT_DOUBLE_EQ(0.0, value);
    const Eigen::ArrayXd const_data{ data };
    const double const_value{ rhbm_gem::eigen_helper::GetStandardDeviation(const_data) };
    EXPECT_DOUBLE_EQ(0.0, const_value);
}

TEST(EigenHelperTest, InvertibleMatrixReturnsCorrectInverse)
{
    Eigen::Matrix2d matrix;
    matrix << 1.0, 2.0, 3.0, 4.0;
    const auto result{ rhbm_gem::eigen_helper::GetInverseMatrix(matrix) };
    const auto expected{ matrix.inverse() };
    EXPECT_TRUE(result.isApprox(expected));
}

TEST(EigenHelperTest, SingularMatrixReturnsPseudoInverse)
{
    Eigen::Matrix2d matrix;
    matrix << 1.0, 2.0, 2.0, 4.0;
    const auto result{ rhbm_gem::eigen_helper::GetInverseMatrix(matrix) };
    const auto expected{ matrix.completeOrthogonalDecomposition().pseudoInverse() };
    EXPECT_TRUE(result.isApprox(expected));
}

TEST(EigenHelperTest, DynamicInvertibleMatrixReturnsCorrectInverse)
{
    Eigen::MatrixXd matrix(3, 3);
    matrix << 2.0, -1.0, 0.0,
              -1.0, 2.0, -1.0,
              0.0, -1.0, 2.0;
    const auto result{ rhbm_gem::eigen_helper::GetInverseMatrix(matrix) };
    const auto expected{ matrix.inverse() };
    EXPECT_TRUE(result.isApprox(expected));
}

TEST(EigenHelperTest, DynamicSingularMatrixReturnsPseudoInverse)
{
    Eigen::MatrixXd matrix(2, 3);
    matrix << 1.0, 2.0, 3.0,
              2.0, 4.0, 6.0;
    const auto result{ rhbm_gem::eigen_helper::GetInverseMatrix(matrix) };
    const auto expected{ matrix.completeOrthogonalDecomposition().pseudoInverse() };
    EXPECT_TRUE(result.isApprox(expected));
    EXPECT_EQ(result.rows(), matrix.cols());
    EXPECT_EQ(result.cols(), matrix.rows());
}

TEST(EigenHelperTest, OneByOneNonZeroMatrixReturnsInverse)
{
    Eigen::Matrix<double, 1, 1> matrix;
    matrix << 2.0;
    const auto result{ rhbm_gem::eigen_helper::GetInverseMatrix(matrix) };
    EXPECT_DOUBLE_EQ(result(0, 0), 0.5);
}

TEST(EigenHelperTest, OneByOneZeroMatrixReturnsZero)
{
    Eigen::Matrix<double, 1, 1> matrix;
    matrix << 0.0;
    const auto result{ rhbm_gem::eigen_helper::GetInverseMatrix(matrix) };
    EXPECT_DOUBLE_EQ(result(0, 0), 0.0);
}

TEST(EigenHelperTest, IdentityMatrixRemainsUnchangedAfterInverse)
{
    Eigen::Matrix3d matrix{ Eigen::Matrix3d::Identity() };
    const auto result{ rhbm_gem::eigen_helper::GetInverseMatrix(matrix) };
    EXPECT_TRUE(result.isApprox(Eigen::Matrix3d::Identity()));
}

TEST(EigenHelperTest, InverseDiagonalHandlesZeros)
{
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> diag(4);
    diag.diagonal() << 2.0, 0.0, -4.0, 5.0;
    auto inv{ rhbm_gem::eigen_helper::GetInverseDiagonalMatrix(diag) };
    EXPECT_DOUBLE_EQ(inv.diagonal()[0], 0.5);
    EXPECT_DOUBLE_EQ(inv.diagonal()[1], 0.0);
    EXPECT_DOUBLE_EQ(inv.diagonal()[2], -0.25);
    EXPECT_DOUBLE_EQ(inv.diagonal()[3], 0.2);
}

TEST(EigenHelperTest, ZeroDiagonalMatrixRemainsZeroAfterInverse)
{
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> diag(3);
    diag.diagonal().setZero();
    auto inv{ rhbm_gem::eigen_helper::GetInverseDiagonalMatrix(diag) };
    EXPECT_DOUBLE_EQ(inv.diagonal()[0], 0.0);
    EXPECT_DOUBLE_EQ(inv.diagonal()[1], 0.0);
    EXPECT_DOUBLE_EQ(inv.diagonal()[2], 0.0);
}
