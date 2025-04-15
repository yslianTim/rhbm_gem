#include "AtomicInfoHelper.hpp"

#include <cmath>

double AtomicInfoHelper::GausModelFunction(double * x, double * par)
{
    double pi{ 3.14159265358979323846264338327950288 };
    double tau_square{ par[1]*par[1] };
    double y{ par[0] * std::pow(2.0*pi*tau_square, -1.5) * std::exp(-x[0]*x[0]/(2.0*tau_square)) };
    return y;
}

