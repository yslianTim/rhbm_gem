#include <gtest/gtest.h>

#include <initializer_list>
#include <stdexcept>

#include "HRLAlphaTrainer.hpp"

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
} // namespace

TEST(HRLAlphaTrainerTest, EvaluateAlphaRWithExactLinearDataReturnsZeroError)
{
    const std::vector<Eigen::VectorXd> data_list{
        MakeVector({ 1.0, 0.0, 1.0 }),
        MakeVector({ 1.0, 1.0, 3.0 }),
        MakeVector({ 1.0, 2.0, 5.0 }),
        MakeVector({ 1.0, 3.0, 7.0 }),
        MakeVector({ 1.0, 4.0, 9.0 }),
        MakeVector({ 1.0, 5.0, 11.0 })
    };

    const auto error_sum_list{
        HRLAlphaTrainer::EvaluateAlphaR(data_list, 3, { 0.0, 0.5 })
    };

    EXPECT_NEAR(0.0, error_sum_list(0), 1e-12);
    EXPECT_NEAR(0.0, error_sum_list(1), 1e-12);
}

TEST(HRLAlphaTrainerTest, EvaluateAlphaGWithIdenticalBetasReturnsZeroError)
{
    const std::vector<Eigen::VectorXd> beta_list(6, MakeVector({ 1.5, -0.5 }));

    const auto error_sum_list{
        HRLAlphaTrainer::EvaluateAlphaG(beta_list, 4, { 0.0, 0.5, 1.0 })
    };

    EXPECT_NEAR(0.0, error_sum_list(0), 1e-12);
    EXPECT_NEAR(0.0, error_sum_list(1), 1e-12);
    EXPECT_NEAR(0.0, error_sum_list(2), 1e-12);
}

TEST(HRLAlphaTrainerTest, EvaluateAlphaRRejectsInvalidSubsetSize)
{
    const std::vector<Eigen::VectorXd> data_list{
        MakeVector({ 1.0, 0.0, 1.0 }),
        MakeVector({ 1.0, 1.0, 3.0 })
    };

    EXPECT_THROW(HRLAlphaTrainer::EvaluateAlphaR(data_list, 0, { 0.0 }), std::invalid_argument);
    EXPECT_THROW(HRLAlphaTrainer::EvaluateAlphaR(data_list, 3, { 0.0 }), std::invalid_argument);
}

TEST(HRLAlphaTrainerTest, EvaluateAlphaGWithSeedIsDeterministic)
{
    const std::vector<Eigen::VectorXd> beta_list{
        MakeVector({ 1.0, 2.0 }),
        MakeVector({ 1.5, 2.5 }),
        MakeVector({ 2.0, 3.0 }),
        MakeVector({ 2.5, 3.5 }),
        MakeVector({ 3.0, 4.0 }),
        MakeVector({ 3.5, 4.5 })
    };
    HRLExecutionOptions options;
    options.random_seed = 42U;

    const auto first{ HRLAlphaTrainer::EvaluateAlphaG(beta_list, 3, { 0.0, 0.5 }, options) };
    const auto second{ HRLAlphaTrainer::EvaluateAlphaG(beta_list, 3, { 0.0, 0.5 }, options) };

    EXPECT_TRUE(first.isApprox(second, 1e-12));
}
