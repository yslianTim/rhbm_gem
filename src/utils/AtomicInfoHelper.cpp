#include "AtomicInfoHelper.hpp"

double AtomicInfoHelper::GausModelFunction(double * x, double * par)
{
    double tau_square{ par[1]*par[1] };
    double y{ par[0] * pow(2.0*M_PI*tau_square, -1.5) * exp(-x[0]*x[0]/(2.0*tau_square)) };
    return y;
}

