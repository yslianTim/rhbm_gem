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
