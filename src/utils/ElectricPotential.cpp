#include "ElectricPotential.hpp"
#include "ChemicalDataHelper.hpp"
#include "GlobalEnumClass.hpp"
#include "Logger.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

const std::unordered_map<Element, std::array<double, 5>> ElectricPotential::a_neutral_par_map
{
    {Element::HYDROGEN,  {   0.0088,   0.0449,   0.1481,   0.2356,   0.0914 }}, // s up to 6.0 A^-1
    {Element::CARBON,    {   0.0489,   0.2091,   0.7537,   1.1420,   0.3555 }},
    {Element::NITROGEN,  {   0.0267,   0.1328,   0.5301,   1.1020,   0.4215 }},
    {Element::OXYGEN,    {   0.0365,   0.1729,   0.5805,   0.8814,   0.3121 }},
    {Element::SODIUM,    {   0.1260,   0.6442,   0.8893,   1.8197,   1.2988 }},
    {Element::MAGNESIUM, {   0.1130,   0.5575,   0.9046,   2.1580,   1.4735 }},
    {Element::PHOSPHORUS,{   0.1005,   0.4615,   1.0663,   2.5854,   1.2725 }},
    {Element::SULFUR,    {   0.0915,   0.4312,   1.0847,   2.4671,   1.0852 }},
    {Element::CHLORINE,  {   0.0799,   0.3891,   1.0037,   2.3332,   1.0507 }},
    {Element::CALCIUM,   {   0.2355,   0.9916,   2.3959,   3.7252,   2.5647 }},
    {Element::IRON,      {   0.1929,   0.8239,   1.8689,   2.3694,   1.9060 }},
    {Element::COPPER,    {   0.3501,   1.6558,   1.9582,   0.2134,   1.4109 }},
    //{Element::COPPER,    { 0.312E+0, 0.812E+0, 0.111E+1, 0.794E+0, 0.257E+0 }}, // Cu1+
    //{Element::COPPER,    { 0.224E+0, 0.544E+0, 0.970E+0, 0.727E+0, 0.182E+0 }}, // Cu2+
    {Element::ZINC,      {   0.1780,   0.8096,   1.6744,   1.9499,   1.4495 }}
};
const std::unordered_map<Element, std::array<double, 5>> ElectricPotential::b_neutral_par_map
{
    {Element::HYDROGEN,  {   0.1152,   1.0867,   4.9755,  16.5591,  43.2743 }}, // s up to 6.0 A^-1
    {Element::CARBON,    {   0.1140,   1.0825,   5.4281,  17.8811,  51.1341 }},
    {Element::NITROGEN,  {   0.0541,   0.5165,   2.8207,  10.6297,  34.3764 }},
    {Element::OXYGEN,    {   0.0652,   0.6184,   2.9449,   9.6298,  28.2194 }},
    {Element::SODIUM,    {   0.1684,   1.7150,   8.8386,  50.8265, 147.2073 }},
    {Element::MAGNESIUM, {   0.1356,   1.3579,   6.9255,  32.3165,  92.1138 }},
    {Element::PHOSPHORUS,{   0.0977,   0.9084,   4.9654,  18.5471,  54.3648 }},
    {Element::SULFUR,    {   0.0838,   0.7788,   4.3462,  15.5846,  44.6365 }},
    {Element::CHLORINE,  {   0.0694,   0.6443,   3.5351,  12.5058,  35.8633 }},
    {Element::CALCIUM,   {   0.1742,   1.8329,   8.8407,  47.4583, 134.9613 }},
    {Element::IRON,      {   0.1087,   1.0806,   4.7637,  22.8500,  76.7309 }},
    {Element::COPPER,    {   0.1867,   1.9917,  11.3396,  53.2619,  63.2520 }},
    //{Element::COPPER,    { 0.201E+0, 0.131E+1, 0.380E+1, 0.105E+2, 0.282E+2 }}, // Cu1+
    //{Element::COPPER,    { 0.145E+0, 0.933E+0, 0.269E+1, 0.711E+1, 0.194E+2 }}, // Cu2+
    {Element::ZINC,      {   0.0876,   0.8650,   3.8612,  18.8726,  64.7016 }}
};

const std::unordered_map<Element, std::array<double, 5>> ElectricPotential::a_positive_par_map
{
    {Element::HYDROGEN,  {      0.0,      0.0,      0.0,      0.0,      0.0 }},
    {Element::CARBON,    { 2.079e-2, 9.266e-2, 2.949e-1, 6.812e-1, 3.304e-1 }},
    {Element::NITROGEN,  { 2.296e-2, 1.004e-1, 3.289e-1, 6.546e-1, 2.733e-1 }},
    {Element::OXYGEN,    { 2.439e-2, 1.036e-1, 3.360e-1, 6.112e-1, 2.447e-1 }},
    {Element::SODIUM,    {   0.1260,   0.6442,   0.8893,   1.8197,   1.2988 }}, // neutral
    {Element::MAGNESIUM, {   0.1130,   0.5575,   0.9046,   2.1580,   1.4735 }}, // neutral
    {Element::PHOSPHORUS,{   0.1005,   0.4615,   1.0663,   2.5854,   1.2725 }}, // neutral
    {Element::SULFUR,    { 6.232e-2, 3.129e-1, 6.541e-1, 1.742e+0, 9.377e-1 }},
    {Element::CHLORINE,  {   0.0799,   0.3891,   1.0037,   2.3332,   1.0507 }}, // neutral
    {Element::CALCIUM,   {   0.2355,   0.9916,   2.3959,   3.7252,   2.5647 }}, // neutral
    {Element::IRON,      {   0.1929,   0.8239,   1.8689,   2.3694,   1.9060 }}, // neutral
    {Element::COPPER,    {   0.3501,   1.6558,   1.9582,   0.2134,   1.4109 }}, // neutral
    {Element::ZINC,      {   0.1780,   0.8096,   1.6744,   1.9499,   1.4495 }}  // neutral
};

const std::unordered_map<Element, std::array<double, 5>> ElectricPotential::b_positive_par_map
{
    {Element::HYDROGEN,  {      0.0,      0.0,      0.0,      0.0,      0.0 }},
    {Element::CARBON,    { 5.950e-2, 5.359e-1, 2.760e+0, 9.283e+0, 2.442e+1 }},
    {Element::NITROGEN,  { 5.522e-2, 4.910e-1, 2.402e+0, 7.751e+0, 2.051e+1 }},
    {Element::OXYGEN,    { 5.082e-2, 4.390e-1, 2.036e+0, 6.407e+0, 1.710e+1 }},
    {Element::SODIUM,    {   0.1684,   1.7150,   8.8386,  50.8265, 147.2073 }}, // neutral
    {Element::MAGNESIUM, {   0.1356,   1.3579,   6.9255,  32.3165,  92.1138 }}, // neutral
    {Element::PHOSPHORUS,{   0.0977,   0.9084,   4.9654,  18.5471,  54.3648 }}, // neutral
    {Element::SULFUR,    { 6.149e-2, 5.785e-1, 2.848e+0, 1.107e+1, 2.978e+1 }},
    {Element::CHLORINE,  {   0.0694,   0.6443,   3.5351,  12.5058,  35.8633 }}, // neutral
    {Element::CALCIUM,   {   0.1742,   1.8329,   8.8407,  47.4583, 134.9613 }}, // neutral
    {Element::IRON,      {   0.1087,   1.0806,   4.7637,  22.8500,  76.7309 }}, // neutral
    {Element::COPPER,    {   0.1867,   1.9917,  11.3396,  53.2619,  63.2520 }}, // neutral
    {Element::ZINC,      {   0.0876,   0.8650,   3.8612,  18.8726,  64.7016 }}  // neutral
};

const std::unordered_map<Element, std::array<double, 5>> ElectricPotential::a_negative_par_map
{
    {Element::HYDROGEN,  {      0.0,      0.0,      0.0,      0.0,      0.0 }},
    {Element::CARBON,    { 2.248e-1, 8.254e-1, 1.796e+0, 1.690e+0, 6.994e-1 }},
    {Element::NITROGEN,  { 2.192e-1, 7.256e-1, 1.398e+0, 1.245e+0, 4.381e-1 }},
    {Element::OXYGEN,    { 2.236e-1, 6.923e-1, 1.176e+0, 9.354e-1, 2.821e-1 }},
    {Element::SODIUM,    {   0.1260,   0.6442,   0.8893,   1.8197,   1.2988 }}, // neutral
    {Element::MAGNESIUM, {   0.1130,   0.5575,   0.9046,   2.1580,   1.4735 }}, // neutral
    {Element::PHOSPHORUS,{   0.1005,   0.4615,   1.0663,   2.5854,   1.2725 }}, // neutral
    {Element::SULFUR,    { 4.496e-1, 9.810e-1, 2.598e+0, 2.717e+0, 8.614e-1 }},
    {Element::CHLORINE,  {   0.0799,   0.3891,   1.0037,   2.3332,   1.0507 }}, // neutral
    {Element::CALCIUM,   {   0.2355,   0.9916,   2.3959,   3.7252,   2.5647 }}, // neutral
    {Element::IRON,      {   0.1929,   0.8239,   1.8689,   2.3694,   1.9060 }}, // neutral
    {Element::COPPER,    {   0.3501,   1.6558,   1.9582,   0.2134,   1.4109 }}, // neutral
    {Element::ZINC,      {   0.1780,   0.8096,   1.6744,   1.9499,   1.4495 }}  // neutral
};

const std::unordered_map<Element, std::array<double, 5>> ElectricPotential::b_negative_par_map
{
    {Element::HYDROGEN,  {      0.0,      0.0,      0.0,      0.0,      0.0 }},
    {Element::CARBON,    { 5.518e-1, 4.308e+0, 1.600e+1, 5.196e+1, 1.708e+2 }},
    {Element::NITROGEN,  { 4.784e-1, 3.389e+0, 1.171e+1, 3.604e+1, 1.125e+2 }},
    {Element::OXYGEN,    { 4.372e-1, 2.918e+0, 9.670e+0, 2.868e+1, 8.489e+1 }},
    {Element::SODIUM,    {   0.1684,   1.7150,   8.8386,  50.8265, 147.2073 }}, // neutral
    {Element::MAGNESIUM, {   0.1356,   1.3579,   6.9255,  32.3165,  92.1138 }}, // neutral
    {Element::PHOSPHORUS,{   0.0977,   0.9084,   4.9654,  18.5471,  54.3648 }}, // neutral
    {Element::SULFUR,    { 4.656e-1, 3.259e+0, 1.233e+1, 3.583e+1, 1.055e+2 }},
    {Element::CHLORINE,  {   0.0694,   0.6443,   3.5351,  12.5058,  35.8633 }}, // neutral
    {Element::CALCIUM,   {   0.1742,   1.8329,   8.8407,  47.4583, 134.9613 }}, // neutral
    {Element::IRON,      {   0.1087,   1.0806,   4.7637,  22.8500,  76.7309 }}, // neutral
    {Element::COPPER,    {   0.1867,   1.9917,  11.3396,  53.2619,  63.2520 }}, // neutral
    {Element::ZINC,      {   0.0876,   0.8650,   3.8612,  18.8726,  64.7016 }}  // neutral
};

ElectricPotential::ElectricPotential() :
    m_model_choice{ ModelChoice::SINGLE_GAUS }, m_blurring_width{ 0.5 }
{

}

void ElectricPotential::SetModelChoice(int value)
{
    m_model_choice = CheckModelChoice(value);
}

void ElectricPotential::SetBlurringWidth(double value)
{
    m_blurring_width = value;
}

double ElectricPotential::GetPotentialValue(
    Element element, double distance, double charge, double amplitude, double width) const
{
    switch (m_model_choice)
    {
        case ModelChoice::SINGLE_GAUS:
            return CalculateSingleGausModel(element, distance);
        case ModelChoice::FIVE_GAUS_CHARGE:
            return CalculateFiveGausChargeModel(element, distance, charge);
        case ModelChoice::SINGLE_GAUS_USER:
            return CalculateSingleGausUserModel(distance, amplitude, width);
        default:
            Logger::Log(LogLevel::Warning, "ElectricPotential::GetPotentialValue: ModelChoice not set.");
            return 0.0;
    }
}

const std::array<double, 5> & ElectricPotential::GetModelParameterAList(
    Element element, int delta_z) const
{
    switch(delta_z)
    {
        case  0:
            if (a_neutral_par_map.find(element) == a_neutral_par_map.end())
            {
                throw std::out_of_range(
                    "ElectricPotential::GetModelParameterAList element not found "
                    + std::to_string(static_cast<int>(element))
                );
            }
            return a_neutral_par_map.at(element);
        case  1:
            if (a_positive_par_map.find(element) == a_positive_par_map.end())
            {
                throw std::out_of_range(
                    "ElectricPotential::GetModelParameterAList element not found "
                    + std::to_string(static_cast<int>(element))
                );
            }
            return a_positive_par_map.at(element);
        case -1:
            if (a_negative_par_map.find(element) == a_negative_par_map.end())
            {
                throw std::out_of_range(
                    "ElectricPotential::GetModelParameterAList element not found "
                    + std::to_string(static_cast<int>(element))
                );
            }
            return a_negative_par_map.at(element);
        default:
            throw std::out_of_range(
                "ElectricPotential::GetModelParameterAList delta_z out of range "
                + std::to_string(delta_z)
            );
    }
}

const std::array<double, 5> & ElectricPotential::GetModelParameterBList(
    Element element, int delta_z) const
{
    switch(delta_z)
    {
        case  0:
            if (b_neutral_par_map.find(element) == b_neutral_par_map.end())
            {
                throw std::out_of_range(
                    "ElectricPotential::GetModelParameterBList element not found "
                    + std::to_string(static_cast<int>(element))
                );
            }
            return b_neutral_par_map.at(element);
        case  1:
            if (b_positive_par_map.find(element) == b_positive_par_map.end())
            {
                throw std::out_of_range(
                    "ElectricPotential::GetModelParameterBList element not found "
                    + std::to_string(static_cast<int>(element))
                );
            }
            return b_positive_par_map.at(element);
        case -1:
            if (b_negative_par_map.find(element) == b_negative_par_map.end())
            {
                throw std::out_of_range(
                    "ElectricPotential::GetModelParameterBList element not found "
                    + std::to_string(static_cast<int>(element))
                );
            }
            return b_negative_par_map.at(element);
        default:
            throw std::out_of_range(
                "ElectricPotential::GetModelParameterBList delta_z out of range "
                + std::to_string(delta_z)
            );
    }
}

ElectricPotential::ModelChoice ElectricPotential::CheckModelChoice(int value) const
{
    switch (value)
    {
        case static_cast<uint8_t>(ModelChoice::SINGLE_GAUS):
            return ModelChoice::SINGLE_GAUS;
        case static_cast<uint8_t>(ModelChoice::FIVE_GAUS_CHARGE):
            return ModelChoice::FIVE_GAUS_CHARGE;
        case static_cast<uint8_t>(ModelChoice::SINGLE_GAUS_USER):
            return ModelChoice::SINGLE_GAUS_USER;
        default:
            throw std::out_of_range(
                "ElectricPotential::ModelChoice out of range " + std::to_string(value)
            );
    }
}

double ElectricPotential::CalculateSingleGausModel(
    Element element, double distance) const
{
    auto atomic_number{ ChemicalDataHelper::GetAtomicNumber(element) };
    auto inv_atomic_number{ 1.0/static_cast<double>(atomic_number) };
    auto sigma_total_square{ inv_atomic_number * inv_atomic_number + m_blurring_width * m_blurring_width };
    auto distance_square{ distance * distance };
    auto exp_index{ -distance_square/(2.0 * sigma_total_square) };
    return std::pow(2.0 * M_PI * sigma_total_square, -1.5) * std::exp(exp_index);
}

double ElectricPotential::CalculateSingleGausUserModel(
    double distance, double amplitude, double width) const
{
    auto sigma_square{ width * width };
    auto distance_square{ distance * distance };
    auto exp_index{ -distance_square/(2.0 * sigma_square) };
    return amplitude * std::pow(2.0 * M_PI * sigma_square, -1.5) * std::exp(exp_index);
}


double ElectricPotential::CalculateFiveGausChargeModel(
    Element element, double distance, double charge) const
{
    auto potential_0{ CalculateFiveGausChargeIntrinsicTerm(element, distance, 0) };
    auto potential_p{ CalculateFiveGausChargeIntrinsicTerm(element, distance, 1) };
    auto potential_n{ CalculateFiveGausChargeIntrinsicTerm(element, distance,-1) };
    auto potential_delta{ CalculateFiveGausChargeDeltaTerm(distance, 1.0) };
    potential_p += potential_delta;
    potential_n -= potential_delta;
    if (charge >= 0.0)
    {
        return potential_0 + charge * (potential_p - potential_0);
    }
    else
    {
        return potential_0 + charge * (potential_0 - potential_n);
    }
}

double ElectricPotential::CalculateFiveGausChargeIntrinsicTerm(
    Element element, double distance, int delta_z) const
{
    const auto & a{ GetModelParameterAList(element, delta_z) };
    const auto & b{ GetModelParameterBList(element, delta_z) };
    auto potential_total{ 0.0 };
    auto scalar{ 8.0 * M_PI * M_PI };
    auto r_square{ distance * distance };
    auto blurring_width_square{ m_blurring_width * m_blurring_width };
    for (size_t i = 0; i < 5; i++)
    {
        auto factor{ b[i]/scalar + blurring_width_square };
        potential_total += a[i] * std::pow(2.0 * M_PI * factor, -1.5) * std::exp(-r_square/factor/2.0);
    }
    potential_total *= F_0;
    return potential_total;
}

double ElectricPotential::CalculateFiveGausChargeDeltaTerm(
    double distance, double charge) const
{
    auto blurring_width{ m_blurring_width };
    if (blurring_width == 0.00) return F_1 * charge/distance;
    if (blurring_width <= 1.0e-5) blurring_width = 1.0e-5;
    if (distance > 3.0) return 0.0;
    if (distance < 1.0e-5) return F_1 * charge * std::sqrt(2.0/M_PI) / blurring_width;
    return F_1 * charge/distance * std::erf(distance/blurring_width/std::sqrt(2.0));
}
