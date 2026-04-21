#include <gtest/gtest.h>

#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <rhbm_gem/utils/math/NumericValidation.hpp>

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

TEST(NumericValidationTest, RequirePositiveAcceptsPositiveSignedAndUnsignedValues)
{
    EXPECT_EQ(3, rg::NumericValidation::RequirePositive(3, "count"));
    EXPECT_EQ(4u, rg::NumericValidation::RequirePositive(4u, "sample_count"));
}

TEST(NumericValidationTest, RequirePositiveRejectsZeroAndNegativeValues)
{
    ExpectInvalidArgumentMessage(
        []()
        {
            (void)rg::NumericValidation::RequirePositive(0, "count");
        },
        "count must be positive.");
    ExpectInvalidArgumentMessage(
        []()
        {
            (void)rg::NumericValidation::RequirePositive(-1, "count");
        },
        "count must be positive.");
}

TEST(NumericValidationTest, RequireFiniteRejectsNonFiniteValue)
{
    ExpectInvalidArgumentMessage(
        []()
        {
            (void)rg::NumericValidation::RequireFinite(
                std::numeric_limits<double>::infinity(),
                "alpha");
        },
        "alpha must be finite.");
}

TEST(NumericValidationTest, RequireFiniteWithContextPrefixesOriginalMessage)
{
    ExpectInvalidArgumentMessage(
        []()
        {
            (void)rg::NumericValidation::RequireFinite(
                std::numeric_limits<double>::infinity(),
                "alpha",
                "Invalid alpha input.");
        },
        "Invalid alpha input. Details: alpha must be finite."
    );
}

TEST(NumericValidationTest, RequireAllFiniteAcceptsFiniteVector)
{
    const std::vector<double> values{ 0.0, 1.0, 2.0 };

    EXPECT_EQ(&values, &rg::NumericValidation::RequireAllFinite(values, "values"));
}

TEST(NumericValidationTest, RequireAllFiniteRejectsNonFiniteVectorValue)
{
    const std::vector<double> values{ 0.0, std::numeric_limits<double>::infinity() };

    ExpectInvalidArgumentMessage(
        [&]()
        {
            (void)rg::NumericValidation::RequireAllFinite(values, "values");
        },
        "values must contain only finite values."
    );
}

TEST(NumericValidationTest, RequireAllFiniteWithContextPrefixesOriginalMessage)
{
    const std::vector<double> values{ 0.0, std::numeric_limits<double>::infinity() };

    ExpectInvalidArgumentMessage(
        [&]()
        {
            (void)rg::NumericValidation::RequireAllFinite(
                values,
                "values",
                "Invalid value list.");
        },
        "Invalid value list. Details: values must contain only finite values."
    );
}

TEST(NumericValidationTest, RequireFinitePositiveRejectsNonPositiveAndNonFiniteValues)
{
    ExpectInvalidArgumentMessage(
        []()
        {
            (void)rg::NumericValidation::RequireFinitePositive(0.0, "width");
        },
        "width must be finite and positive.");
    ExpectInvalidArgumentMessage(
        []()
        {
            (void)rg::NumericValidation::RequireFinitePositive(
                std::numeric_limits<double>::quiet_NaN(),
                "width");
        },
        "width must be finite and positive.");
}

TEST(NumericValidationTest, RequireFiniteNonNegativeAcceptsZeroAndRejectsNegative)
{
    EXPECT_DOUBLE_EQ(0.0, rg::NumericValidation::RequireFiniteNonNegative(0.0, "alpha"));
    ExpectInvalidArgumentMessage(
        []()
        {
            (void)rg::NumericValidation::RequireFiniteNonNegative(-0.1, "alpha");
        },
        "alpha must be finite and non-negative.");
}

TEST(NumericValidationTest, RequireFiniteNonNegativeRangeAcceptsBoundaryValues)
{
    const auto range{
        rg::NumericValidation::RequireFiniteNonNegativeRange(0.0, 0.0, "distance range")
    };

    EXPECT_DOUBLE_EQ(0.0, range.first);
    EXPECT_DOUBLE_EQ(0.0, range.second);
}

TEST(NumericValidationTest, RequireFiniteNonNegativeRangeRejectsInvalidInputs)
{
    ExpectInvalidArgumentMessage(
        []()
        {
            (void)rg::NumericValidation::RequireFiniteNonNegativeRange(
                2.0,
                1.0,
                "distance range");
        },
        "distance range must be finite and satisfy 0 <= min <= max.");
    ExpectInvalidArgumentMessage(
        []()
        {
            (void)rg::NumericValidation::RequireFiniteNonNegativeRange(
                0.0,
                std::numeric_limits<double>::infinity(),
                "distance range");
        },
        "distance range must be finite and satisfy 0 <= min <= max.");
}

TEST(NumericValidationTest, RequireFiniteInclusiveRangeAcceptsBoundsAndRejectsOutsideValues)
{
    EXPECT_DOUBLE_EQ(
        180.0,
        rg::NumericValidation::RequireFiniteInclusiveRange(180.0, 0.0, 180.0, "angle"));

    ExpectInvalidArgumentMessage(
        []()
        {
            (void)rg::NumericValidation::RequireFiniteInclusiveRange(181.0, 0.0, 180.0, "angle");
        },
        "angle must be finite and within [0, 180].");
    ExpectInvalidArgumentMessage(
        []()
        {
            (void)rg::NumericValidation::RequireFiniteInclusiveRange(
                std::numeric_limits<double>::quiet_NaN(),
                0.0,
                180.0,
                "angle");
        },
        "angle must be finite and within [0, 180].");
}
