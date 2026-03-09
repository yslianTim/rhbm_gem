#include <gtest/gtest.h>
#include <cmath>
#include <array>
#include <tuple>
#include <stdexcept>

#include <rhbm_gem/utils/math/ElectricPotential.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>

namespace
{
constexpr double PI{ 3.14159265358979323846264338327950288 };
constexpr double F_0{ 47.87764193 };
constexpr double F_1{ 14.39964393 };

double ComputeSingleGaus(Element element, double distance, double blurring_width)
{
    const double atomic_number{ static_cast<double>(ChemicalDataHelper::GetAtomicNumber(element)) };
    const double inv_atomic_number{ 1.0 / atomic_number };
    const double sigma_total_square{ inv_atomic_number * inv_atomic_number +
                                     blurring_width * blurring_width };
    const double exp_index{ -(distance * distance) / (2.0 * sigma_total_square) };
    return std::pow(2.0 * PI * sigma_total_square, -1.5) * std::exp(exp_index);
}

double ComputeFiveGausChargeIntrinsicTerm(
    ElectricPotential & pot, Element element, double distance,
    int delta_z, double blurring_width)
{
    const auto & a{ pot.GetModelParameterAList(element, delta_z) };
    const auto & b{ pot.GetModelParameterBList(element, delta_z) };
    double potential_total{ 0.0 };
    const double scalar{ 8.0 * PI * PI };
    const double r_square{ distance * distance };
    const double blurring_width_square{ blurring_width * blurring_width };
    for (size_t i = 0; i < 5; ++i)
    {
        const double factor{ b[i] / scalar + blurring_width_square };
        potential_total +=
            a[i] * std::pow(2.0 * PI * factor, -1.5) *
            std::exp(-r_square / factor / 2.0);
    }
    potential_total *= F_0;
    return potential_total;
}

double ComputeFiveGausCharge(ElectricPotential & pot, Element element,
                             double distance, double charge,
                             double blurring_width)
{
    const double potential_0{
        ComputeFiveGausChargeIntrinsicTerm(pot, element, distance, 0, blurring_width)
    };
    double potential_p{
        ComputeFiveGausChargeIntrinsicTerm(pot, element, distance, 1, blurring_width)
    };
    double potential_n{
        ComputeFiveGausChargeIntrinsicTerm(pot, element, distance, -1, blurring_width)
    };
    const double potential_delta{ F_1 / distance };
    potential_p += potential_delta;
    potential_n -= potential_delta;
    if (charge >= 0.0)
    {
        return potential_0 + charge * (potential_p - potential_0);
    }
    return potential_0 + charge * (potential_0 - potential_n);
}

double ComputeSingleGausUser(double distance, double amplitude, double width)
{
    const double sigma_square{ width * width };
    const double exp_index{ -(distance * distance) / (2.0 * sigma_square) };
    return amplitude * std::pow(2.0 * PI * sigma_square, -1.5) *
           std::exp(exp_index);
}

double ComputeFiveGausChargeDeltaTerm(double distance, double charge, double blurring_width)
{
    double bw{ blurring_width };
    if (bw == 0.0) return F_1 * charge / distance;
    if (bw <= 1.0e-5) bw = 1.0e-5;
    if (distance > 3.0) return 0.0;
    if (distance < 1.0e-5) return F_1 * charge * std::sqrt(2.0 / PI) / bw;
    return F_1 * charge / distance * std::erf(distance / bw / std::sqrt(2.0));
}
} // namespace

TEST(ElectricPotentialTest, SetModelChoiceSelectsModel)
{
    ElectricPotential potential;

    // Model 0: SINGLE_GAUS
    potential.SetBlurringWidth(0.0);
    potential.SetModelChoice(0);
    const double expected_single_gaus{ ComputeSingleGaus(Element::CARBON, 1.0, 0.0) };
    const double actual_single_gaus{ potential.GetPotentialValue(Element::CARBON, 1.0, 0.0) };
    EXPECT_DOUBLE_EQ(expected_single_gaus, actual_single_gaus);

    // Model 1: FIVE_GAUS_CHARGE
    potential.SetBlurringWidth(0.0);
    potential.SetModelChoice(1);
    const double expected_five_gaus{
        ComputeFiveGausCharge(potential, Element::CARBON, 1.0, 1.0, 0.0)
    };
    const double actual_five_gaus{
        potential.GetPotentialValue(Element::CARBON, 1.0, 1.0)
    };
    EXPECT_DOUBLE_EQ(expected_five_gaus, actual_five_gaus);

    // Model 2: SINGLE_GAUS_USER
    potential.SetModelChoice(2);
    const double amplitude{ 1.0 };
    const double width{ 1.0 };
    const double expected_single_gaus_user{
        ComputeSingleGausUser(1.0, amplitude, width)
    };
    const double actual_single_gaus_user{
        potential.GetPotentialValue(Element::CARBON, 1.0, 0.0, amplitude, width)
    };
    EXPECT_DOUBLE_EQ(expected_single_gaus_user, actual_single_gaus_user);
}

TEST(ElectricPotentialTest, DefaultConstructorSetsExpectedState)
{
    ElectricPotential potential; // uses default SINGLE_GAUS model and blurring width 0.5

    constexpr double distance{ 1.0 };
    const Element element{ Element::CARBON };

    // Value calculated by the model under test using defaults
    double computed{ potential.GetPotentialValue(element, distance, 0.0) };

    // Expected value computed independently using blurring width 0.5
    int atomic_number{ ChemicalDataHelper::GetAtomicNumber(element) };
    double inv_atomic_number{ 1.0 / static_cast<double>(atomic_number) };
    double sigma_total_square{ inv_atomic_number * inv_atomic_number + 0.5 * 0.5 };
    double distance_square{ distance * distance };
    double exp_index{ -distance_square / (2.0 * sigma_total_square) };
    double expected{
        std::pow(2.0 * M_PI * sigma_total_square, -1.5) * std::exp(exp_index)
    };

    EXPECT_NEAR(expected, computed, 1e-12);
}

TEST(ElectricPotentialTest, SingleGaussianCarbon)
{
    ElectricPotential potential;
    potential.SetModelChoice(0);        // ModelChoice::SINGLE_GAUS
    potential.SetBlurringWidth(0.1);    // Known blurring width

    constexpr double distance{ 1.0 };
    const Element element{ Element::CARBON };

    // Value calculated by the model under test
    double computed{ potential.GetPotentialValue(element, distance, 0.0) };

    // Expected value using the formula from CalculateSingleGausModel
    int atomic_number{ ChemicalDataHelper::GetAtomicNumber(element) };
    double inv_atomic_number{ 1.0 / static_cast<double>(atomic_number) };
    double sigma_total_square{ inv_atomic_number * inv_atomic_number + 0.1 * 0.1 };
    double distance_square{ distance * distance };
    double exp_index{ -distance_square / (2.0 * sigma_total_square) };
    double expected{
        std::pow(2.0 * M_PI * sigma_total_square, -1.5) * std::exp(exp_index)
    };

    EXPECT_NEAR(expected, computed, 1e-12);
}

TEST(ElectricPotentialTest, ReturnsKnownParameterListsForCarbon)
{
    ElectricPotential ep;
    const std::array<double, 5> expected_a0{ 0.0489, 0.2091, 0.7537, 1.1420, 0.3555 };
    EXPECT_EQ(expected_a0, ep.GetModelParameterAList(Element::CARBON, 0));

    const std::array<double, 5> expected_b0{ 0.1140, 1.0825, 5.4281, 17.8811, 51.1341 };
    EXPECT_EQ(expected_b0, ep.GetModelParameterBList(Element::CARBON, 0));

    const std::array<double, 5> expected_a1{ 2.079e-2, 9.266e-2, 2.949e-1, 6.812e-1, 3.304e-1 };
    EXPECT_EQ(expected_a1, ep.GetModelParameterAList(Element::CARBON, 1));

    const std::array<double, 5> expected_b1{ 5.950e-2, 5.359e-1, 2.760e+0, 9.283e+0, 2.442e+1 };
    EXPECT_EQ(expected_b1, ep.GetModelParameterBList(Element::CARBON, 1));

    const std::array<double, 5> expected_a_neg{ 2.248e-1, 8.254e-1, 1.796e+0, 1.690e+0, 6.994e-1 };
    EXPECT_EQ(expected_a_neg, ep.GetModelParameterAList(Element::CARBON, -1));

    const std::array<double, 5> expected_b_neg{ 5.518e-1, 4.308e+0, 1.600e+1, 5.196e+1, 1.708e+2 };
    EXPECT_EQ(expected_b_neg, ep.GetModelParameterBList(Element::CARBON, -1));
}

TEST(ElectricPotentialTest, UnsupportedElementThrows)
{
    ElectricPotential ep;
    EXPECT_THROW(ep.GetModelParameterAList(Element::ALUMINUM, 0), std::out_of_range);
    EXPECT_THROW(ep.GetModelParameterBList(Element::ALUMINUM, 0), std::out_of_range);
}

TEST(ElectricPotentialTest, UnsupportedElementPositiveNegativeThrows)
{
    ElectricPotential ep;
    EXPECT_THROW(ep.GetModelParameterAList(Element::ALUMINUM, 1), std::out_of_range);
    EXPECT_THROW(ep.GetModelParameterBList(Element::ALUMINUM, 1), std::out_of_range);
    EXPECT_THROW(ep.GetModelParameterAList(Element::ALUMINUM, -1), std::out_of_range);
    EXPECT_THROW(ep.GetModelParameterBList(Element::ALUMINUM, -1), std::out_of_range);
}

TEST(ElectricPotentialTest, InvalidDeltaZThrows)
{
    ElectricPotential ep;
    EXPECT_THROW(ep.GetModelParameterAList(Element::CARBON, 2), std::out_of_range);
    EXPECT_THROW(ep.GetModelParameterBList(Element::CARBON, 2), std::out_of_range);
}

TEST(ElectricPotentialTest, SetModelChoiceThrowsOnInvalidValue)
{
    ElectricPotential potential;
    EXPECT_THROW(potential.SetModelChoice(3), std::out_of_range);
}

TEST(ElectricPotentialTest, FiveGausChargeInterpolatesWithCharge)
{
    ElectricPotential potential;
    potential.SetModelChoice(1); // ModelChoice::FIVE_GAUS_CHARGE
    constexpr double blurring_width{ 0.2 };
    potential.SetBlurringWidth(blurring_width);
    const Element element{ Element::CARBON };
    constexpr double distance{ 1.0 };
    const double actual_neutral{ potential.GetPotentialValue(element, distance, 0.0) };
    const double actual_half{ potential.GetPotentialValue(element, distance, 0.5) };
    const double expected_neutral{
        ComputeFiveGausChargeIntrinsicTerm(potential, element, distance, 0, blurring_width)
    };
    const double expected_positive{
        ComputeFiveGausChargeIntrinsicTerm(potential, element, distance, 1, blurring_width) +
        ComputeFiveGausChargeDeltaTerm(distance, 1.0, blurring_width)
    };
    const double expected_half{ expected_neutral + 0.5 * (expected_positive - expected_neutral) };
    EXPECT_NEAR(expected_neutral, actual_neutral, 1e-12);
    EXPECT_NEAR(expected_half, actual_half, 1e-12);
}

TEST(ElectricPotentialTest, NegativeChargeFiveGaus)
{
    ElectricPotential potential;
    potential.SetBlurringWidth(0.0);
    potential.SetModelChoice(1); // ModelChoice::FIVE_GAUS_CHARGE
    const double expected{
        ComputeFiveGausCharge(potential, Element::CARBON, 1.0, -1.0, 0.0)
    };
    const double actual{
        potential.GetPotentialValue(Element::CARBON, 1.0, -1.0)
    };
    EXPECT_DOUBLE_EQ(expected, actual);
}

TEST(ElectricPotentialTest, NeutralChargeFiveGausBlurringWidth)
{
    ElectricPotential potential;
    potential.SetModelChoice(1); // ModelChoice::FIVE_GAUS_CHARGE

    const Element element{ Element::CARBON };
    constexpr double distance{ 1.0 };
    potential.SetBlurringWidth(0.0);
    const double neutral_no_blur{ potential.GetPotentialValue(element, distance, 0.0) };
    potential.SetBlurringWidth(0.5);
    const double neutral_with_blur{ potential.GetPotentialValue(element, distance, 0.0) };
    const double expected_no_blur{
        ComputeFiveGausChargeIntrinsicTerm(potential, element, distance, 0, 0.0)
    };
    const double expected_with_blur{
        ComputeFiveGausChargeIntrinsicTerm(potential, element, distance, 0, 0.5)
    };

    const double actual_diff{ neutral_with_blur - neutral_no_blur };
    const double expected_diff{ expected_with_blur - expected_no_blur };
    EXPECT_NEAR(expected_no_blur, neutral_no_blur, 1e-12);
    EXPECT_NEAR(expected_with_blur, neutral_with_blur, 1e-12);
    EXPECT_NEAR(expected_diff, actual_diff, 1e-12);
    EXPECT_NE(0.0, actual_diff);
}

using AmplitudeWidthPair = std::pair<double, double>;

class SingleGausUserModelTest : public ::testing::Test,
                                public ::testing::WithParamInterface<AmplitudeWidthPair>
{
};

TEST_P(SingleGausUserModelTest, ComputesExpectedPotential)
{
    const auto [amplitude, width]{ GetParam() };
    ElectricPotential potential;
    potential.SetModelChoice(2); // ModelChoice::SINGLE_GAUS_USER
    constexpr double distance{ 1.0 };
    const double expected{ ComputeSingleGausUser(distance, amplitude, width) };
    const double actual{
        potential.GetPotentialValue(Element::CARBON, distance, 0.0, amplitude, width)
    };
    EXPECT_NEAR(expected, actual, 1e-12);
}

INSTANTIATE_TEST_SUITE_P(SingleGausUserCases, SingleGausUserModelTest,
    ::testing::Values(
        AmplitudeWidthPair{ 1.0, 1.0 },
        AmplitudeWidthPair{ 0.5, 0.2 },
        AmplitudeWidthPair{ 2.0, 0.8 },
        AmplitudeWidthPair{ 1.5, 1.5 }
    )
);

using FiveGausChargeDeltaTermCase = std::tuple<double, double, double>;

class FiveGausChargeDeltaTermTest : public ::testing::Test,
                                   public ::testing::WithParamInterface<FiveGausChargeDeltaTermCase>
{
};

TEST_P(FiveGausChargeDeltaTermTest, ComputesExpectedDeltaTerm)
{
    const auto [blurring_width, distance, charge]{ GetParam() };
    ElectricPotential potential;
    potential.SetModelChoice(1); // FIVE_GAUS_CHARGE
    potential.SetBlurringWidth(blurring_width);

    const Element element{ Element::CARBON };

    const double expected{
        ComputeFiveGausChargeDeltaTerm(distance, charge, blurring_width)
    };

    const double with_charge{
        potential.GetPotentialValue(element, distance, charge)
    };
    const double without_charge{
        potential.GetPotentialValue(element, distance, 0.0)
    };

    const double intrinsic_p{
        ComputeFiveGausChargeIntrinsicTerm(potential, element, distance, 1, blurring_width)
    };
    const double intrinsic_0{
        ComputeFiveGausChargeIntrinsicTerm(potential, element, distance, 0, blurring_width)
    };
    const double intrinsic_n{
        ComputeFiveGausChargeIntrinsicTerm(potential, element, distance, -1, blurring_width)
    };

    const double intrinsic_delta{
        charge >= 0.0 ? (intrinsic_p - intrinsic_0) : -(intrinsic_0 - intrinsic_n)
    };

    const double delta_actual{
        (with_charge - without_charge) - intrinsic_delta
    };
    EXPECT_NEAR(expected, delta_actual, 1e-12);
}

INSTANTIATE_TEST_SUITE_P(FiveGausChargeDeltaTermCases, FiveGausChargeDeltaTermTest,
    ::testing::Values(
        FiveGausChargeDeltaTermCase{ 0.0,     1.0,     1.0 },    // blurring_width == 0
        FiveGausChargeDeltaTermCase{ 1.0e-6,  1.0,     1.0 },    // blurring_width <= 1e-5
        FiveGausChargeDeltaTermCase{ 0.5,     4.0,     1.0 },    // distance > 3.0
        FiveGausChargeDeltaTermCase{ 0.5,     1.0e-6,  1.0 },    // distance < 1e-5
        FiveGausChargeDeltaTermCase{ 0.5,     1.0,     1.0 },    // general case
        FiveGausChargeDeltaTermCase{ 0.5,     1.0,    -1.0 }     // negative charge
    )
);
