#pragma once

#if __INTELLISENSE__
#undef __ARM_NEON
#undef __ARM_NEON__
#endif

#include <vector>
#include <tuple>
#include <Eigen/Dense>

class GausLinearTransformHelper
{

public:
    GausLinearTransformHelper(void) = default;
    ~GausLinearTransformHelper() = default;

    static Eigen::VectorXd BuildLinearModelDataVector(double gaus_x, double gaus_y);
    static Eigen::VectorXd BuildGausModel(const Eigen::VectorXd & linear_model);
    static std::tuple<Eigen::VectorXd, Eigen::VectorXd> BuildGausModelWithVariance(
        const Eigen::VectorXd & linear_model, const Eigen::MatrixXd & covariance_matrix);

private:

};