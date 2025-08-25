#pragma once

#if __INTELLISENSE__
#undef __ARM_NEON
#undef __ARM_NEON__
#endif

#include <cstddef>
#include <vector>
#include <string>
#include <tuple>
#include <Eigen/Dense>

class HRLModelHelper
{
    const int m_basis_size;
    int m_member_size; ///< The size of atomic member (number of members)  \f$ I \f$
    int m_maximum_iteration; ///< The upper limit of iteration (default = 100)
    double m_tolerance; ///< The tolerance for algorithms (default = \f$ 10^{-5} \f$)
    double m_omega_sum; ///< The sum of member weights \f$ \Omega^{(t)} = \sum_{0}^{I-1}\omega_{i}^{(t)} \f$
    double m_omega_h;   ///< The reciprocal of sum of inverse member weights
    double m_weight_data_min;
    double m_weight_member_min;

    /**
     * @brief   The list of selected data sample size \f$ = n_{i} \f$
     *          for each member \f$ i = 0\cdots I-1 \f$
     *          (length of list \f$ = I \f$)
     * @see     BuildDataArray()
     */
    std::vector<int> m_data_size_list;

    /**
     * @brief   The list of member information
     *          for each member \f$ i = 0\cdots I-1 \f$
     *          (length of list \f$ = I \f$)
     * @see     BuildDataArray()
     */
    std::vector<std::string> m_member_info_list;

    /**
     * @brief   The list of model basis matrix \f$ \boldsymbol{X}_{i} \f$
     *          for each member \f$ i = 0\cdots I-1 \f$
     *          (length of list \f$ = I \f$,
     *          dimension of matrix \f$ = \left[n_{i}\times b\right] \f$)
     * @see     BuildDataArray()
     */
    std::vector<Eigen::MatrixXd> m_X_list;

    /**
     * @brief   The list of data vector \f$ \boldsymbol{n}_{i} \f$
     *          for each member \f$ i = 0\cdots I-1 \f$
     *          (length of list \f$ = I \f$,
     *          dimension of vector \f$ = \left[n_{i}\times 1\right] \f$)
     * @see     BuildDataArray()
     */
    std::vector<Eigen::VectorXd> m_y_list;

    /**
     * @brief   The list of diagonal weight matrix \f$ \boldsymbol{W}_{i} \f$
     *          for each member \f$ i = 0\cdots I-1 \f$
     *          (length of list \f$ = I \f$,
     *          dimension of matrix \f$ = \left[n_{i}\times n_{i}\right] \f$)
     * @see     CalculateDataWeight()
     */
    std::vector<Eigen::DiagonalMatrix<double, Eigen::Dynamic>> m_W_list;

    /**
     * @brief   The list of data covariance matrix \f$ \boldsymbol{\Sigma}_{\boldsymbol{y}_{i}} \f$
     *          for each member \f$ i = 0\cdots I-1 \f$
     *          (length of list \f$ = I \f$,
     *          dimension of matrix \f$ = \left[n_{i}\times n_{i}\right] \f$)
     * @see     CalculateDataCovariance()
     */
    std::vector<Eigen::DiagonalMatrix<double, Eigen::Dynamic>> m_capital_sigma_list;

    /**
     * @brief   The list of posterior covariance matrix \f$ \tilde{\boldsymbol{\Sigma}}_{i} \f$
     *          for each member \f$ i = 0\cdots I-1 \f$
     *          (length of list \f$ = I \f$,
     *          dimension of matrix \f$ = \left[p\times p\right] \f$)
     * @see     AlgorithmWEB()
     */
    std::vector<Eigen::MatrixXd> m_capital_sigma_posterior_list;

    /**
     * @brief   The list of weighted prior covariance matrix \f$ \boldsymbol{\Lambda}_{i} \f$
     *          for each member \f$ i = 0\cdots I-1 \f$
     *          (length of list \f$ = I \f$,
     *          dimension of matrix \f$ = \left[p\times p\right] \f$)
     * @see     AlgorithmMuMDPDE()
     */
    std::vector<Eigen::MatrixXd> m_capital_lambda_list;

    /**
     * @brief   The double float array of square of data variance \f$ = \sigma^{2}_{i} \f$
     *          in data normal distribution \f$ = \mathcal{N}(0, \sigma^{2}_{i}) \f$
     *          (dimension of array \f$ = \left[I\times 1\right] \f$)
     * @see     CalculateDataVarianceSquare()
     */
    Eigen::ArrayXd m_sigma_square_array;

    /**
     * @brief   The double float array of member weight \f$ = \omega_{i} \f$
     *          (dimension of array \f$ = \left[I\times 1\right] \f$)
     * @see     CalculateMemberWeight()
     */
    Eigen::ArrayXd m_omega_array;

    /**
     * @brief   The double float array of statistical distance \f$ d_{i} \f$
     *          between estimated prior and posterior parameter
     * @see     CalculateStatisticalDistance()
     */
    Eigen::ArrayXd m_statistical_distance_array; // [member_size x 1]

    /**
     * @brief   The boolean flag array of member outlier
     *          (dimension of array \f$ = \left[I\times 1\right] \f$)
     * @see     LabelOutlierMember()
     */
    Eigen::Array<bool, Eigen::Dynamic, 1> m_outlier_flag_array;

    /**
     * @brief   The prior covariance matrix \f$ \boldsymbol{\Lambda} \f$
     *          (dimension of matrix \f$ = \left[p\times p\right] \f$)
     * @see     CalculateMemberCovariance()
     */
    Eigen::MatrixXd m_capital_lambda;


    Eigen::VectorXd m_mu_iter, m_mu_MDPDE, m_mu_prior, m_mu_mean; // [basis_size x 1]
    Eigen::MatrixXd m_beta_iter_array;
    Eigen::MatrixXd m_beta_OLS_array, m_beta_MDPDE_array, m_beta_posterior_array; // [basis_size x member_size]

public:
    static constexpr int DEFAULT_MAXIMUM_ITERATION{ 100 };
    static constexpr double DEFAULT_TOLERANCE{ 1.0e-5 };
    static constexpr double DEFAULT_WEIGHT_DATA_MIN{ 1.0e-8 };
    static constexpr double DEFAULT_WEIGHT_MEMBER_MIN{ 1.0e-2 };

    HRLModelHelper(void) = delete;
    explicit HRLModelHelper(int basis_size, int member_size);
    ~HRLModelHelper() = default;

    void SetThreadSize(int thread_size);
    void SetDataArray(const std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> & data_array);
    void RunEstimation(double alpha_r, double alpha_g);
    void SetMaximumIteration(int size);
    void SetTolerance(double value);
    bool GetOutlierFlag(int id) const;
    double GetStatisticalDistance(int id) const;
    double GetSigmaSquare(int id) const;
    double GetMemberWeight(int id) const;
    const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & GetDataWeightMatrix(int id) const;
    const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & GetDataCovarianceMatrix(int id) const;
    const Eigen::MatrixXd & GetCapitalSigmaMatrixPosterior(int id) const;
    const Eigen::MatrixXd & GetCapitalLambdaMatrix(void) const { return m_capital_lambda; }
    const Eigen::VectorXd & GetMuVectorPrior(void) const { return m_mu_prior; }
    const Eigen::VectorXd & GetMuVectorMDPDE(void) const { return m_mu_MDPDE; }
    const Eigen::VectorXd & GetMuVectorMean(void) const { return m_mu_mean; }
    Eigen::Ref<const Eigen::VectorXd> GetBetaMatrixPosterior(int id) const;
    Eigen::Ref<const Eigen::VectorXd> GetBetaMatrixMDPDE(int id) const;
    Eigen::Ref<const Eigen::VectorXd> GetBetaMatrixOLS(int id) const;

private:
    void ValidateMemberId(int id) const;
    void Initialization(void);
    void AlgorithmBetaMDPDE(double alpha_r);
    void AlgorithmMuMDPDE(double alpha_g);
    void AlgorithmWEB(void);
    void CalculateBetaByMDPDE(int member_id);
    void CalculateMuByMDPDE(void);
    void CalculateDataVarianceSquare(int member_id, double alpha_r);
    void CalculateDataCovariance(int member_id);
    void CalculateMemberCovariance(double alpha_g);
    void CalculateMemberWeight(double alpha_g);
    void CalculateDataWeight(int member_id, double alpha_r);
    void CalculateStatisticalDistance(void);
    void LabelOutlierMember(void);

};
