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

    static double GetGaussianResponseAtDistance(double distance, double width, int dimension=3);
    static Eigen::VectorXd GetTylorSeriesBasisVector(double distance, const Eigen::VectorXd & model_par);
    static Eigen::VectorXd BuildLinearModelDataVector(double x, double y, const Eigen::VectorXd & model_par, int basis_dimension=2);
    static Eigen::VectorXd BuildGaus2DModel(const Eigen::VectorXd & linear_model);
    static Eigen::VectorXd BuildGaus3DModel(const Eigen::VectorXd & linear_model);
    static Eigen::VectorXd BuildGaus3DModel(const Eigen::VectorXd & linear_model, const Eigen::VectorXd & model_par);
    static std::tuple<Eigen::VectorXd, Eigen::VectorXd> BuildGaus2DModelWithVariance(
        const Eigen::VectorXd & linear_model, const Eigen::MatrixXd & covariance_matrix);
    static std::tuple<Eigen::VectorXd, Eigen::VectorXd> BuildGaus3DModelWithVariance(
        const Eigen::VectorXd & linear_model, const Eigen::MatrixXd & covariance_matrix);
    
};
