#include <gtest/gtest.h>

#include <stdexcept>
#include <string>

#include <rhbm_gem/utils/math/EigenValidation.hpp>

namespace rg = rhbm_gem;

namespace {

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

} // namespace

TEST(EigenValidationTest, RequireSameSizeAcceptsSameSizedVectors)
{
    const Eigen::VectorXd lhs{ Eigen::VectorXd::Zero(3) };
    const Eigen::VectorXd rhs{ Eigen::VectorXd::Zero(3) };

    EXPECT_NO_THROW(rg::EigenValidation::RequireSameSize(lhs, rhs, "vector pair"));
}

TEST(EigenValidationTest, RequireSameSizeRejectsMismatchedVectors)
{
    const Eigen::VectorXd lhs{ Eigen::VectorXd::Zero(3) };
    const Eigen::VectorXd rhs{ Eigen::VectorXd::Zero(2) };

    ExpectInvalidArgumentMessage(
        [&]()
        {
            rg::EigenValidation::RequireSameSize(lhs, rhs, "vector pair");
        },
        "vector pair must have matching sizes."
    );
}

TEST(EigenValidationTest, RequireVectorSizeWithContextPrefixesOriginalMessage)
{
    const Eigen::VectorXd vector{ Eigen::VectorXd::Zero(2) };

    ExpectInvalidArgumentMessage(
        [&]()
        {
            rg::EigenValidation::RequireVectorSize(
                vector,
                3,
                "vector",
                "Vector shape is invalid.");
        },
        "Vector shape is invalid. Details: vector must have size 3."
    );
}

TEST(EigenValidationTest, RequireRowsAcceptsExpectedRowCount)
{
    const Eigen::MatrixXd matrix{ Eigen::MatrixXd::Zero(3, 2) };

    EXPECT_NO_THROW(rg::EigenValidation::RequireRows(matrix, 3, "matrix"));
}

TEST(EigenValidationTest, RequireRowsWithContextPrefixesOriginalMessage)
{
    const Eigen::MatrixXd matrix{ Eigen::MatrixXd::Zero(2, 3) };

    ExpectInvalidArgumentMessage(
        [&]()
        {
            rg::EigenValidation::RequireRows(
                matrix,
                4,
                "matrix",
                "Row validation failed.");
        },
        "Row validation failed. Details: matrix must have 4 rows."
    );
}

TEST(EigenValidationTest, RequireColsAcceptsExpectedColumnCount)
{
    const Eigen::MatrixXd matrix{ Eigen::MatrixXd::Zero(3, 2) };

    EXPECT_NO_THROW(rg::EigenValidation::RequireCols(matrix, 2, "matrix"));
}

TEST(EigenValidationTest, RequireColsWithContextPrefixesOriginalMessage)
{
    const Eigen::MatrixXd matrix{ Eigen::MatrixXd::Zero(2, 3) };

    ExpectInvalidArgumentMessage(
        [&]()
        {
            rg::EigenValidation::RequireCols(
                matrix,
                4,
                "matrix",
                "Column validation failed.");
        },
        "Column validation failed. Details: matrix must have 4 columns."
    );
}

TEST(EigenValidationTest, RequireSameSizeWithContextPrefixesOriginalMessage)
{
    const Eigen::VectorXd lhs{ Eigen::VectorXd::Zero(3) };
    const Eigen::VectorXd rhs{ Eigen::VectorXd::Zero(2) };

    ExpectInvalidArgumentMessage(
        [&]()
        {
            rg::EigenValidation::RequireSameSize(
                lhs,
                rhs,
                "vector pair",
                "Pair size check failed.");
        },
        "Pair size check failed. Details: vector pair must have matching sizes."
    );
}
