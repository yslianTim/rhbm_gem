#include "AtomicInfoHelper.hpp"

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

const std::unordered_map<Element, int> AtomicInfoHelper::atomic_number_map
{
    {Element::HYDROGEN, 1}, {Element::CARBON,    6}, {Element::NITROGEN,   7},
    {Element::OXYGEN,   8}, {Element::SODIUM,   11}, {Element::MAGNESIUM, 12},
    {Element::SULFUR,  16}, {Element::CHLORINE, 17}, {Element::CALCIUM,   20},
    {Element::IRON,    26}, {Element::ZINC,     30}
};

double AtomicInfoHelper::GausModelFunction(double * x, double * par)
{
    double tau_square{ par[1]*par[1] };
    double y{ par[0] * std::pow(2.0*M_PI*tau_square, -1.5) * std::exp(-x[0]*x[0]/(2.0*tau_square)) };
    return y;
}

