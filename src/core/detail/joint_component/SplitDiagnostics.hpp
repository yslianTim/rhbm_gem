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
    std::string status{"Unavailable"},reason,null_reference_status{"NotAssessed"},
        null_reference_failure_reason,high_budget_reference_status{"NotRun"},
        independent_replay_status{"NotRun"},width_stationarity_status{"NotRun"};
    std::optional<double> reoptimization_prediction_drift,reoptimization_log_width_difference_inf,
        reoptimization_scaled_ac_difference,seed_prediction_robustness,seed_log_width_difference_inf,
        seed_scaled_ac_difference,vector_decomposition_max_error,high_budget_width_difference_inf,
        high_budget_scaled_ac_difference,high_budget_prediction_difference_relative;
    std::string reoptimization_prediction_endpoint{"null-reference"},
        reoptimization_log_width_endpoint{"null-reference"},reoptimization_scaled_ac_endpoint{"null-reference"},
        seed_prediction_endpoint,seed_log_width_endpoint,seed_scaled_ac_endpoint;
    std::size_t qualified_seed_count{},assessed_direction_count{},incomplete_direction_count{};
    double null_search_seconds{},null_search_cpu_seconds{};
    std::size_t null_profile_evaluations{},null_accepted_updates{},null_replay_checks{},
        workspace_estimate_bytes{};
    std::vector<SplitDirectionDiagnostic> directions;
};

struct SplitDiagnosticBlock
{
    std::string block_id;
    LocalSplitDiagnostics diagnostics;
};

enum class SplitDiagnosticsMode {Disabled,Minimal,ReuseWorkspace};

class SplitDiagnosticsAccumulator
{
public:
    SplitDiagnosticsAccumulator(const Vector &,const Vector &,const Vector &,const Vector &,const Vector &,
        const Vector &,const Vector &,double,bool,std::string,std::size_t,bool=false);
    void AddEndpoint(const std::string &,const std::string &,const Vector &,const Vector &,const Vector &,
        bool,const std::string & = {});
    void FinishDirection(const std::string &,bool,const std::string & = {});
    LocalSplitDiagnostics Finish() const;

private:
    const Vector & eta0_,&beta0_,&prediction0_,&eta_null_,&beta_null_,&prediction_null_,&reference_column_norms_;
    double scale_{};
    bool null_qualified_{},inputs_valid_{};
    bool reuse_workspace_{};
    std::size_t expected_directions_{};
    LocalSplitDiagnostics result_;
    Vector prediction_delta_workspace_,width_delta_workspace_,beta_delta_workspace_,
        decomposition_left_workspace_,decomposition_right_workspace_;
};

} // namespace rhbm_gem::core::joint_component
