#pragma once

#include <vector>
#include <string>
#include <tuple>
#include <Eigen/Dense>

class HRLModelTester
{
    const int m_basis_size;
    int m_member_size;
    double m_x_min, m_x_max;
    double m_data_error_sigma;
    Eigen::VectorXd m_model_par_prior;
    Eigen::VectorXd m_model_par_sigma;
    std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> m_data_array;

public:
    HRLModelTester(void) = delete;
    explicit HRLModelTester(int basis_size, int member_size);
    ~HRLModelTester() = default;

    void SetDataErrorSigma(double data_error_sigma) { m_data_error_sigma = data_error_sigma; }
    void SetModelParametersPrior(const Eigen::VectorXd & model_par_prior);
    void SetModelParametersSigma(const Eigen::VectorXd & model_par_sigma);
    void BuildDataArray(void);
    std::vector<std::tuple<double, double>> BuildRandomSamplingEntry(
        size_t sampling_entry_size, const Eigen::VectorXd & model_par);
    bool RunTest(double alpha_r, double alpha_g);

private:
    bool CheckModelParametersDimension(const Eigen::VectorXd & model_par);

};
