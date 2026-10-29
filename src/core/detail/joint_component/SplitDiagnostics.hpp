#pragma once
#include "Numerics.hpp"

namespace rhbm_gem::core::joint_component {

struct SplitDirectionDiagnostic
{
    std::string direction,status{"NotAssessed"},reason;
    std::optional<double> seed_prediction_robustness,seed_log_width_difference_inf,seed_scaled_ac_difference;
    std::string prediction_endpoint,log_width_endpoint,scaled_ac_endpoint;
};

struct LocalSplitDiagnostics
{
    std::string status{"Unavailable"},reason;
    std::optional<double> reoptimization_prediction_drift,reoptimization_log_width_difference_inf,
        reoptimization_scaled_ac_difference,seed_prediction_robustness,seed_log_width_difference_inf,
        seed_scaled_ac_difference,vector_decomposition_max_error;
    std::string reoptimization_prediction_endpoint{"null-reference"},
        reoptimization_log_width_endpoint{"null-reference"},reoptimization_scaled_ac_endpoint{"null-reference"},
        seed_prediction_endpoint,seed_log_width_endpoint,seed_scaled_ac_endpoint;
    std::size_t qualified_seed_count{},assessed_direction_count{},incomplete_direction_count{};
    std::vector<SplitDirectionDiagnostic> directions;
};

class SplitDiagnosticsAccumulator
{
public:
    SplitDiagnosticsAccumulator(const Vector &,const Vector &,const Vector &,const Vector &,const Vector &,
        const Vector &,const Vector &,double,bool,std::string,std::size_t);
    void AddEndpoint(const std::string &,const std::string &,const Vector &,const Vector &,const Vector &,
        bool,const std::string & = {});
    void FinishDirection(const std::string &,bool,const std::string & = {});
    LocalSplitDiagnostics Finish() const;

private:
    const Vector & eta0_,&beta0_,&prediction0_,&eta_null_,&beta_null_,&prediction_null_,&reference_column_norms_;
    double scale_{};
    bool null_qualified_{},inputs_valid_{};
    std::size_t expected_directions_{};
    LocalSplitDiagnostics result_;
};

} // namespace rhbm_gem::core::joint_component
