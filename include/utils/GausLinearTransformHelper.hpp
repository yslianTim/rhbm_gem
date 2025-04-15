#pragma once

#if __INTELLISENSE__
#undef __ARM_NEON
#undef __ARM_NEON__
#endif

#define _USE_MATH_DEFINES // for MSVC compiler in Windows platform

#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <Eigen/Dense>

class GausLinearTransformHelper
{
    static constexpr double m_two_pi{ 2.0 * M_PI };

public:
    GausLinearTransformHelper(void) = default;
    ~GausLinearTransformHelper() = default;

    static Eigen::VectorXd BuildLinearModelDataVector(double gaus_x, double gaus_y);
    static Eigen::VectorXd BuildGausModel(const Eigen::VectorXd & linear_model);
    static std::tuple<Eigen::VectorXd, Eigen::VectorXd> BuildGausModelWithVariance(const Eigen::VectorXd & linear_model, const Eigen::MatrixXd & covariance_matrix);

private:

};