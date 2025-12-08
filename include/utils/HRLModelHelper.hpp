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
    bool m_quiet_mode;
    const int m_basis_size;
    int m_member_size; ///< The size of atomic member (number of members)  \f$ I \f$
    int m_maximum_iteration; ///< The upper limit of iteration (default = 100)
    double m_tolerance; ///< The tolerance for algorithms (default = \f$ 10^{-5} \f$)
    double m_omega_sum; ///< The sum of member weights \f$ \Omega^{(t)} = \sum_{0}^{I-1}\omega_{i}^{(t)} \f$

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


    Eigen::VectorXd m_mu_MDPDE, m_mu_prior, m_mu_mean; // [basis_size x 1]
    Eigen::MatrixXd m_beta_MDPDE_array, m_beta_posterior_array; // [basis_size x member_size]

public:
    static constexpr int DEFAULT_MAXIMUM_ITERATION{ 100 };
    static constexpr double DEFAULT_TOLERANCE{ 1.0e-5 };
    static constexpr double DEFAULT_WEIGHT_DATA_MIN{ 1.0e-8 };
    static constexpr double DEFAULT_WEIGHT_MEMBER_MIN{ 1.0e-2 };

    HRLModelHelper(void) = delete;
    explicit HRLModelHelper(int basis_size, int member_size);
    ~HRLModelHelper();

    void SetQuietMode(void) { m_quiet_mode = true; }
    void SetThreadSize(int thread_size);
    void SetMaximumIteration(int size);
    void SetTolerance(double value);
    void SetMemberDataEntriesList(const std::vector<std::vector<Eigen::VectorXd>> & data_list);
    void SetMemberBetaMDPDEList(
        const std::vector<Eigen::VectorXd> & beta_list,
        const std::vector<double> & sigma_square_list,
        const std::vector<Eigen::DiagonalMatrix<double, Eigen::Dynamic>> & W_list,
        const std::vector<Eigen::DiagonalMatrix<double, Eigen::Dynamic>> & capital_sigma_list
    );

    std::tuple<Eigen::MatrixXd, Eigen::VectorXd> BuildBasisVectorAndResponseArray(
        const std::vector<Eigen::VectorXd> & data_vector
    );
    
    void RunGroupEstimation(double alpha_g);
    void RunBetaMDPDE(
        const std::vector<Eigen::VectorXd> & data_vector,
        double alpha_r,
        Eigen::VectorXd & beta_OLS,
        Eigen::VectorXd & beta_MDPDE,
        double & sigma_square,
        Eigen::DiagonalMatrix<double, Eigen::Dynamic> & W,
        Eigen::DiagonalMatrix<double, Eigen::Dynamic> & capital_sigma
    );
    void RunMuMDPDE(
        const std::vector<Eigen::VectorXd> & beta_list,
        double alpha_g,
        Eigen::VectorXd & mu_median,
        Eigen::VectorXd & mu_MDPDE,
        Eigen::ArrayXd & omega_array,
        double & omega_sum,
        Eigen::MatrixXd & capital_lambda,
        std::vector<Eigen::MatrixXd> & member_capital_lambda_list
    );
    Eigen::VectorXd RunAlphaRTraining(
        const std::vector<Eigen::VectorXd> & data_list,
        const size_t subset_size,
        const std::vector<double> & alpha_list
    );
    Eigen::VectorXd RunAlphaGTraining(
        const std::vector<Eigen::VectorXd> & data_list,
        const size_t subset_size,
        const std::vector<double> & alpha_list
    );
    static bool AlgorithmBetaMDPDE(
        double alpha_r,
        const Eigen::MatrixXd & X,
        const Eigen::VectorXd & y,
        Eigen::VectorXd & beta_OLS,
        Eigen::VectorXd & beta,
        double & sigma_square,
        Eigen::DiagonalMatrix<double, Eigen::Dynamic> & W,
        Eigen::DiagonalMatrix<double, Eigen::Dynamic> & capital_sigma,
        int iteration = DEFAULT_MAXIMUM_ITERATION,
        double tolerance = DEFAULT_TOLERANCE,
        double weight_min = DEFAULT_WEIGHT_DATA_MIN,
        bool quite_mode = false
    );
    static bool AlgorithmMuMDPDE(
        double alpha_g,
        const Eigen::MatrixXd & beta_array,
        Eigen::VectorXd & mu_median,
        Eigen::VectorXd & mu,
        Eigen::ArrayXd & omega_array,
        double & omega_sum,
        Eigen::MatrixXd & capital_lambda,
        std::vector<Eigen::MatrixXd> & member_capital_lambda_list,
        int iteration = DEFAULT_MAXIMUM_ITERATION,
        double tolerance = DEFAULT_TOLERANCE,
        double weight_min = DEFAULT_WEIGHT_MEMBER_MIN,
        bool quite_mode = false
    );
    static bool AlgorithmWEB(
        const std::vector<Eigen::MatrixXd> & X_list,
        const std::vector<Eigen::VectorXd> & y_list,
        const std::vector<Eigen::DiagonalMatrix<double, Eigen::Dynamic>> & capital_sigma_list,
        const Eigen::VectorXd & mu_MDPDE,
        const std::vector<Eigen::MatrixXd> & member_capital_lambda_list,
        Eigen::VectorXd & mu_prior,
        Eigen::MatrixXd & beta_posterior_array,
        std::vector<Eigen::MatrixXd> & capital_sigma_posterior_list,
        bool quite_mode = false
    );
    
    bool GetOutlierFlag(int id) const;
    double GetStatisticalDistance(int id) const;
    Eigen::Ref<const Eigen::VectorXd> GetBetaPosterior(int id) const;
    const Eigen::MatrixXd & GetCapitalSigmaMatrixPosterior(int id) const;
    const Eigen::MatrixXd & GetCapitalLambdaMatrix(void) const { return m_capital_lambda; }
    const Eigen::VectorXd & GetMuVectorPrior(void) const { return m_mu_prior; }
    const Eigen::VectorXd & GetMuVectorMDPDE(void) const { return m_mu_MDPDE; }
    const Eigen::VectorXd & GetMuVectorMean(void) const { return m_mu_mean; }

private:
    void ValidateBasisSize(int size) const;
    void ValidateMemberSize(int size) const;
    void ValidateMemberId(int id) const;
    static void ValidateAlpha(double alpha);
    void PrepareDataSubset(
        const std::vector<Eigen::VectorXd> & data_list,
        size_t subset_size,
        std::vector<std::vector<Eigen::VectorXd>> & data_subset_list
    ) const;
    void PrepareTestAndTrainingDataSet(
        const std::vector<std::vector<Eigen::VectorXd>> & data_subset_list,
        size_t subset_size,
        size_t total_entries_size,
        std::vector<std::vector<Eigen::VectorXd>> & data_set_test,
        std::vector<std::vector<Eigen::VectorXd>> & data_set_training
    ) const;
    Eigen::MatrixXd ConvertBetaListToMatrix(const std::vector<Eigen::VectorXd> & beta_list) const;

    static Eigen::VectorXd CalculateBetaByOLS(
        const Eigen::MatrixXd & X,
        const Eigen::VectorXd & y
    );
    static Eigen::VectorXd CalculateBetaByMDPDE(
        const Eigen::MatrixXd & X,
        const Eigen::VectorXd & y,
        const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & W
    );
    static Eigen::VectorXd CalculateMuByMedian(
        const Eigen::MatrixXd & beta_array
    );
    static Eigen::VectorXd CalculateMuByMDPDE(
        const Eigen::MatrixXd & beta_array,
        const Eigen::ArrayXd & omega_array,
        double omega_sum
    );
    static Eigen::DiagonalMatrix<double, Eigen::Dynamic> CalculateDataWeight(
        double alpha,
        const Eigen::MatrixXd & X,
        const Eigen::VectorXd & y,
        const Eigen::VectorXd & beta,
        double sigma_square,
        double weight_min = DEFAULT_WEIGHT_DATA_MIN
    );
    static double CalculateDataVarianceSquare(
        double alpha,
        const Eigen::MatrixXd & X,
        const Eigen::VectorXd & y,
        const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & W,
        const Eigen::VectorXd & beta
    );
    static Eigen::DiagonalMatrix<double, Eigen::Dynamic> CalculateDataCovariance(
        double sigma_square,
        const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & W
    );
    static Eigen::ArrayXd CalculateMemberWeight(
        double alpha,
        const Eigen::MatrixXd & beta_array,
        const Eigen::VectorXd & mu,
        const Eigen::MatrixXd & capital_lambda,
        double weight_min = DEFAULT_WEIGHT_MEMBER_MIN
    );
    static Eigen::MatrixXd CalculateMemberCovariance(
        double alpha,
        const Eigen::MatrixXd & beta_array,
        const Eigen::VectorXd & mu,
        const Eigen::ArrayXd & omega_array,
        double omega_sum
    );
    static std::vector<Eigen::MatrixXd> CalculateWeightedMemberCovariance(
        const Eigen::MatrixXd & capital_lambda,
        const Eigen::ArrayXd & omega_array
    );
    
    static Eigen::ArrayXd CalculateMemberStatisticalDistance(
        const Eigen::VectorXd & mu_prior,
        const Eigen::MatrixXd & capital_lambda,
        const Eigen::MatrixXd & beta_posterior_array
    );
    static Eigen::Array<bool, Eigen::Dynamic, 1> CalculateOutlierMemberFlag(
        int basis_size,
        const Eigen::ArrayXd & statistical_distance_array
    );

};
