#include <gtest/gtest.h>
#include "core/detail/joint_component/HaloEligibility.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using rhbm_gem::core::JointProblem;
using rhbm_gem::core::JointProblemInput;
using rhbm_gem::core::JointSupport;
using rhbm_gem::core::joint_component::EvaluateFullABCHaloEligibility;

std::vector<JointSupport> Support(std::initializer_list<double> squares)
{
    std::vector<JointSupport> out;
    std::size_t row{};
    for(const double square:squares) out.push_back({row++,square});
    return out;
}
}

TEST(JointHaloEligibilityTest, RejectsZeroSupport)
{
    std::vector<JointSupport> support;
    const auto result=EvaluateFullABCHaloEligibility(support);
    EXPECT_FALSE(result.eligible);
    EXPECT_EQ(result.support_rows,0u);
    EXPECT_EQ(result.distinct_squared_distances,0u);
    EXPECT_EQ(result.reason,"insufficient-observations");
}

TEST(JointHaloEligibilityTest, RejectsSingletonAndTwoRowSupport)
{
    auto singleton=Support({0.0});
    const auto one=EvaluateFullABCHaloEligibility(singleton);
    EXPECT_FALSE(one.eligible);
    EXPECT_EQ(one.support_rows,1u);
    EXPECT_EQ(one.reason,"insufficient-observations");

    auto two_rows=Support({0.0,1.0});
    const auto two=EvaluateFullABCHaloEligibility(two_rows);
    EXPECT_FALSE(two.eligible);
    EXPECT_EQ(two.support_rows,2u);
    EXPECT_EQ(two.reason,"insufficient-observations");
}

TEST(JointHaloEligibilityTest, RejectsThreeRowsWithOneOrTwoExactDistances)
{
    auto equal=Support({1.0,1.0,1.0});
    const auto one=EvaluateFullABCHaloEligibility(equal);
    EXPECT_FALSE(one.eligible);
    EXPECT_EQ(one.distinct_squared_distances,1u);
    EXPECT_EQ(one.reason,"insufficient-radial-information");

    auto two=Support({1.0,1.0,2.0});
    const auto two_distances=EvaluateFullABCHaloEligibility(two);
    EXPECT_FALSE(two_distances.eligible);
    EXPECT_EQ(two_distances.distinct_squared_distances,2u);
    EXPECT_EQ(two_distances.reason,"insufficient-radial-information");
}

TEST(JointHaloEligibilityTest, AcceptsThreeRowsWithThreeExactDistances)
{
    auto support=Support({0.0,1.0,4.0});
    const auto result=EvaluateFullABCHaloEligibility(support);
    EXPECT_TRUE(result.eligible);
    EXPECT_EQ(result.support_rows,3u);
    EXPECT_EQ(result.distinct_squared_distances,3u);
    EXPECT_EQ(result.reason,"eligible");
}

TEST(JointHaloEligibilityTest, ExactNearEqualDistancesAreDeterministicAndHaveNoTolerance)
{
    const double adjacent=std::nextafter(1.0,2.0);
    auto support=Support({1.0,adjacent,4.0,1.0,adjacent,4.0,1.0,adjacent,4.0,1.0,adjacent,4.0});
    auto reversed=support;
    std::reverse(reversed.begin(),reversed.end());
    const auto first=EvaluateFullABCHaloEligibility(support,8);
    const auto second=EvaluateFullABCHaloEligibility(reversed,8);
    EXPECT_TRUE(first.eligible);
    EXPECT_EQ(first.support_rows,12u);
    EXPECT_EQ(first.distinct_squared_distances,3u);
    EXPECT_EQ(first.support_rows,second.support_rows);
    EXPECT_EQ(first.distinct_squared_distances,second.distinct_squared_distances);
    EXPECT_EQ(first.reason,second.reason);
}

TEST(JointHaloEligibilityTest, SelectedTargetWithInsufficientSupportRemainsInExplicitProblem)
{
    JointProblemInput input;
    input.atom_ids={"target"};
    input.row_ids={"row"};
    input.observations={1.0};
    input.support={{{0,0.0}}};
    input.selection_domain.emplace();
    input.selection_domain->target_indices={0};
    const JointProblem problem(std::move(input));
    EXPECT_EQ(problem.Input().atom_ids,(std::vector<std::string>{"target"}));
    EXPECT_EQ(problem.ParameterLayout().full_atoms,(std::vector<std::size_t>{0}));
}

TEST(JointHaloEligibilityTest, ExplicitProblemWithoutSelectionMetadataKeepsFullABCContributors)
{
    JointProblemInput input;
    input.atom_ids={"a","b","c"};
    input.row_ids={"row"};
    input.observations={1.0};
    input.support={{{0,0.0}},{{0,1.0}},{{0,4.0}}};
    const JointProblem problem(std::move(input));
    EXPECT_FALSE(problem.Input().selection_domain.has_value());
    EXPECT_EQ(problem.ParameterLayout().full_atoms,(std::vector<std::size_t>{0,1,2}));
    EXPECT_TRUE(problem.ParameterLayout().groups.empty());
}

TEST(JointHaloEligibilityTest, RejectsConfiguredObservationMinimumBelowThree)
{
    auto support=Support({0.0,1.0,4.0});
    EXPECT_THROW(EvaluateFullABCHaloEligibility(support,2),std::invalid_argument);
}
