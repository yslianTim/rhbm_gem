#pragma once

#include <unordered_map>
#include <cstdint>
#include <array>

#include "GlobalEnumClass.hpp"

class ElectricPotential
{
    enum class ModelChoice : uint8_t
    {
        SINGLE_GAUS      = 0,
        FIVE_GAUS_CHARGE = 1,
        SINGLE_GAUS_USER = 2
    };

    static constexpr double F_0{ 47.87764193 };
    static constexpr double F_1{ 14.39964393 };
    static const std::unordered_map<Element, std::array<double, 5>> a_neutral_par_map;
    static const std::unordered_map<Element, std::array<double, 5>> b_neutral_par_map;
    static const std::unordered_map<Element, std::array<double, 5>> a_positive_par_map;
    static const std::unordered_map<Element, std::array<double, 5>> b_positive_par_map;
    static const std::unordered_map<Element, std::array<double, 5>> a_negative_par_map;
    static const std::unordered_map<Element, std::array<double, 5>> b_negative_par_map;

    ModelChoice m_model_choice;
    double m_blurring_width;

public:
    ElectricPotential(void);
    ~ElectricPotential() = default;
    
    void SetModelChoice(int value);
    void SetBlurringWidth(double value);
    double GetPotentialValue(Element element, double distance, double charge, double amplitude=0.0, double width=0.0) const;
    const std::array<double, 5> & GetModelParameterAList(Element element, int delta_z) const;
    const std::array<double, 5> & GetModelParameterBList(Element element, int delta_z) const;

private:
    ModelChoice CheckModelChoice(int value) const;
    double CalculateSingleGausModel(Element element, double distance) const;
    double CalculateSingleGausUserModel(double distance, double amplitude, double width) const;
    double CalculateFiveGausChargeModel(Element element, double distance, double charge) const;
    double CalculateFiveGausChargeIntrinsicTerm(Element element, double distance, int delta_z) const;
    double CalculateFiveGausChargeDeltaTerm(double distance, double charge) const;

};
