#include <gtest/gtest.h>
#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/joint_component/Numerics.hpp"
#include <cmath>
#include <limits>
#include <string>

namespace {
rhbm_gem::JointInitializationAtom FirstStageHalo()
{
    rhbm_gem::JointInitializationAtom atom;
    atom.id = "42";
    atom.mdpde = {1.25, 0.8, -0.3};
    atom.sample_count = 17;
    atom.native_status = 0;
    atom.first_stage_fit_qualified = true;
    atom.reason = "valid-width";
    return atom;
}
}

TEST(JointFixedHaloInitializationTest, CapturesFirstStageAmplitudeAndWidthWithZeroCharge)
{
    std::string reason;
    const auto fixed = rhbm_gem::core::detail::MakeFixedHaloParameter(
        FirstStageHalo(), 3, "42", reason);
    ASSERT_TRUE(fixed) << reason;
    EXPECT_EQ(fixed->contributor_index, 3u);
    EXPECT_EQ(fixed->atom_id, "42");
    EXPECT_DOUBLE_EQ(fixed->amplitude, 1.25);
    EXPECT_DOUBLE_EQ(fixed->width, 0.8);
    EXPECT_DOUBLE_EQ(fixed->charge, 0.0);
    EXPECT_EQ(fixed->first_stage_sample_count, 17u);
    EXPECT_EQ(fixed->first_stage_status, 0);

    for (const double distance : {0.0, 0.4, 1.5, 2.5})
    {
        const auto basis = rhbm_gem::core::joint_component::EvaluateKernel(
            distance * distance, fixed->width, 2.5);
        const rhbm_gem::GaussianModel3D first_stage_model{fixed->amplitude, fixed->width, 0.0};
        EXPECT_NEAR(first_stage_model.SignalAtDistance(distance),
            fixed->amplitude * basis.gaussian,
            1e-14 * (1.0 + std::abs(fixed->amplitude * basis.gaussian)));
    }
}

TEST(JointFixedHaloInitializationTest, ReportsMissingAndMismatchedIdentity)
{
    std::string reason;
    auto atom = FirstStageHalo();
    EXPECT_FALSE(rhbm_gem::core::detail::MakeFixedHaloParameter(atom, 0, "", reason));
    EXPECT_EQ(reason, "missing-atom-identity");
    EXPECT_FALSE(rhbm_gem::core::detail::MakeFixedHaloParameter(atom, 0, "43", reason));
    EXPECT_EQ(reason, "contributor-identity-mismatch");
}

TEST(JointFixedHaloInitializationTest, RejectsFailedFitAndMedianWidthFallback)
{
    std::string reason;
    auto atom = FirstStageHalo();
    atom.first_stage_fit_qualified = false;
    EXPECT_FALSE(rhbm_gem::core::detail::MakeFixedHaloParameter(atom, 0, "42", reason));
    EXPECT_EQ(reason, "first-stage-fitting-failure");

    atom = FirstStageHalo();
    atom.seed_source = "median-fallback";
    EXPECT_FALSE(rhbm_gem::core::detail::MakeFixedHaloParameter(atom, 0, "42", reason));
    EXPECT_EQ(reason, "median-width-fallback");
}

TEST(JointFixedHaloInitializationTest, RejectsNonFiniteOrInvalidFirstStageAmplitudeAndWidth)
{
    std::string reason;
    auto atom = FirstStageHalo();
    atom.mdpde[0] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_FALSE(rhbm_gem::core::detail::MakeFixedHaloParameter(atom, 0, "42", reason));
    EXPECT_EQ(reason, "non-finite-parameters");

    atom = FirstStageHalo();
    atom.mdpde[0] = -0.1;
    EXPECT_FALSE(rhbm_gem::core::detail::MakeFixedHaloParameter(atom, 0, "42", reason));
    EXPECT_EQ(reason, "invalid-amplitude");

    atom = FirstStageHalo();
    atom.mdpde[1] = 0.0;
    EXPECT_FALSE(rhbm_gem::core::detail::MakeFixedHaloParameter(atom, 0, "42", reason));
    EXPECT_EQ(reason, "invalid-width");
}

TEST(JointFixedHaloInitializationTest, RejectsTargetRecords)
{
    std::string reason;
    auto atom = FirstStageHalo();
    atom.target = true;
    EXPECT_FALSE(rhbm_gem::core::detail::MakeFixedHaloParameter(atom, 0, "42", reason));
    EXPECT_EQ(reason, "target-is-not-fixed-halo");
}
