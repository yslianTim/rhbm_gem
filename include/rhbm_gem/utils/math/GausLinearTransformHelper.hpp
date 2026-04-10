#pragma once

#if __INTELLISENSE__
#undef __ARM_NEON
#undef __ARM_NEON__
#endif

#include <vector>
#include <Eigen/Dense>

#include <rhbm_gem/utils/math/SamplingTypes.hpp>

class GausLinearTransformHelper
{

public:
    GausLinearTransformHelper() = default;
    ~GausLinearTransformHelper() = default;

    static double GetGaussianResponseAtDistance(double distance, double width, int dimension=3);
    static double GetGaussianPesponseAtPoint(
        const Eigen::VectorXd & point, const Eigen::VectorXd & center, double width);
    static double GetGaussianPesponseAtPointWithNeighborhood(
        const Eigen::VectorXd & point, const Eigen::VectorXd & center,
        std::vector<Eigen::VectorXd> & neighbor_center_list, double width);
    static std::vector<Eigen::VectorXd> MapValueTransform(
        const LocalPotentialSampleList & sampling_entries,
        double x_min, double x_max, int basis_size=2
    );
    static Eigen::VectorXd GetTylorSeriesBasisVector(double distance, const Eigen::VectorXd & model_par);
    static Eigen::VectorXd BuildLinearModelDataVector(double x, double y, const Eigen::VectorXd & model_par, int basis_dimension=2);
    static Eigen::VectorXd BuildLinearModelCoefficentVector(double amplitude, double width);
    static Eigen::VectorXd BuildGaus2DModel(const Eigen::VectorXd & linear_model);
    static Eigen::VectorXd BuildGaus3DModel(const Eigen::VectorXd & linear_model);
    static Eigen::VectorXd BuildGaus3DModel(const Eigen::VectorXd & linear_model, const Eigen::VectorXd & model_par);
    static std::tuple<Eigen::VectorXd, Eigen::VectorXd> BuildGaus2DModelWithVariance(
        const Eigen::VectorXd & linear_model, const Eigen::MatrixXd & covariance_matrix);
    static std::tuple<Eigen::VectorXd, Eigen::VectorXd> BuildGaus3DModelWithVariance(
        const Eigen::VectorXd & linear_model, const Eigen::MatrixXd & covariance_matrix);
    
};
