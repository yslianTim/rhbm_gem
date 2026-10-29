#include "SplitDiagnostics.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace rhbm_gem::core::joint_component {
namespace {
bool SameSize(const Vector & left,const Vector & right) {return left.size()==right.size();}
bool Finite(const Vector & value) {return value.size()>0 && value.allFinite();}
double PredictionDenominator(const Vector & prediction) {return std::max(1.0,prediction.norm());}
}

SplitDiagnosticsAccumulator::SplitDiagnosticsAccumulator(const Vector & eta0,const Vector & beta0,
    const Vector & prediction0,const Vector & eta_null,const Vector & beta_null,const Vector & prediction_null,
    const Vector & reference_column_norms,double scale,bool null_qualified,std::string null_failure_reason,
    std::size_t expected_directions,bool reuse_workspace):eta0_(eta0),beta0_(beta0),prediction0_(prediction0),eta_null_(eta_null),
    beta_null_(beta_null),prediction_null_(prediction_null),reference_column_norms_(reference_column_norms),
    scale_(scale),null_qualified_(null_qualified),reuse_workspace_(reuse_workspace),
    expected_directions_(expected_directions)
{
    result_.reason=std::move(null_failure_reason);
    if(!null_qualified_) return;
    inputs_valid_=SameSize(eta0_,eta_null_) && SameSize(beta0_,beta_null_) &&
        SameSize(beta0_,reference_column_norms_) && SameSize(prediction0_,prediction_null_) &&
        Finite(eta0_) && Finite(eta_null_) && Finite(beta0_) && Finite(beta_null_) &&
        Finite(prediction0_) && Finite(prediction_null_) && Finite(reference_column_norms_) &&
        (reference_column_norms_.array()>=0).all() && scale_>0 && std::isfinite(scale_);
    if(!inputs_valid_)
    {
        result_.reason="invalid-or-nonfinite-split-input";
        return;
    }
    if(reuse_workspace_)
    {
        prediction_delta_workspace_.resize(prediction0_.size());
        width_delta_workspace_.resize(eta0_.size());
        beta_delta_workspace_.resize(beta0_.size());
        decomposition_left_workspace_.resize(prediction0_.size());
        decomposition_right_workspace_.resize(prediction0_.size());
    }
    const double denominator=PredictionDenominator(prediction0_);
    result_.reoptimization_prediction_drift=(prediction_null_-prediction0_).norm()/denominator;
    result_.reoptimization_log_width_difference_inf=(eta_null_-eta0_).lpNorm<Eigen::Infinity>();
    result_.reoptimization_scaled_ac_difference=
        (((beta_null_-beta0_).array()*reference_column_norms_.array()).abs()/scale_).maxCoeff();
    result_.vector_decomposition_max_error=0.;
    if(!std::isfinite(*result_.reoptimization_prediction_drift) ||
        !std::isfinite(*result_.reoptimization_log_width_difference_inf) ||
        !std::isfinite(*result_.reoptimization_scaled_ac_difference))
    {
        inputs_valid_=false;
        result_.reason="nonfinite-split-metric";
        result_.reoptimization_prediction_drift.reset();
        result_.reoptimization_log_width_difference_inf.reset();
        result_.reoptimization_scaled_ac_difference.reset();
    }
}

void SplitDiagnosticsAccumulator::AddEndpoint(const std::string & direction,const std::string & identity,
    const Vector & eta,const Vector & beta,const Vector & prediction,bool qualified,const std::string & failure_reason)
{
    if(!inputs_valid_) return;
    auto found=std::find_if(result_.directions.begin(),result_.directions.end(),[&](const auto & item) {
        return item.direction==direction;
    });
    if(found==result_.directions.end())
    {
        result_.directions.push_back({});
        found=std::prev(result_.directions.end());
        found->direction=direction;
        found->status="Pending";
    }
    if(!qualified)
    {
        if(found->reason.empty()) found->reason=failure_reason.empty() ? "endpoint-not-qualified" : failure_reason;
        return;
    }
    if(!SameSize(eta,eta_null_) || !SameSize(beta,beta_null_) || !SameSize(prediction,prediction_null_) ||
        !Finite(eta) || !Finite(beta) || !Finite(prediction))
    {
        if(found->reason.empty()) found->reason="invalid-or-nonfinite-seed-endpoint";
        return;
    }
    double prediction_difference{},width_difference{},ac_difference{},decomposition_error{},decomposition_scale{};
    if(reuse_workspace_)
    {
        prediction_delta_workspace_.noalias()=prediction-prediction_null_;
        width_delta_workspace_.noalias()=eta-eta_null_;
        beta_delta_workspace_.noalias()=beta-beta_null_;
        decomposition_left_workspace_.noalias()=prediction-prediction0_;
        decomposition_right_workspace_.noalias()=prediction_null_-prediction0_;
        decomposition_right_workspace_+=prediction_delta_workspace_;
        prediction_difference=prediction_delta_workspace_.norm()/PredictionDenominator(prediction0_);
        width_difference=width_delta_workspace_.lpNorm<Eigen::Infinity>();
        ac_difference=((beta_delta_workspace_.array()*reference_column_norms_.array()).abs()/scale_).maxCoeff();
        decomposition_error=(decomposition_left_workspace_-decomposition_right_workspace_).norm();
        decomposition_scale=std::max({1.0,decomposition_left_workspace_.norm(),decomposition_right_workspace_.norm()});
    }
    else
    {
        const double denominator=PredictionDenominator(prediction0_);
        prediction_difference=(prediction-prediction_null_).norm()/denominator;
        width_difference=(eta-eta_null_).lpNorm<Eigen::Infinity>();
        ac_difference=(((beta-beta_null_).array()*reference_column_norms_.array()).abs()/scale_).maxCoeff();
        const Vector lhs=prediction-prediction0_;
        const Vector rhs=(prediction_null_-prediction0_)+(prediction-prediction_null_);
        decomposition_error=(lhs-rhs).norm();
        decomposition_scale=std::max({1.0,lhs.norm(),rhs.norm()});
    }
    const double decomposition_tolerance=64.*std::numeric_limits<double>::epsilon()*decomposition_scale;
    if(!std::isfinite(prediction_difference) || !std::isfinite(width_difference) ||
        !std::isfinite(ac_difference) || !std::isfinite(decomposition_error) ||
        decomposition_error>decomposition_tolerance)
    {
        if(found->reason.empty()) found->reason="nonfinite-metric-or-vector-decomposition-failed";
        return;
    }
    ++result_.qualified_seed_count;
    result_.vector_decomposition_max_error=std::max(*result_.vector_decomposition_max_error,decomposition_error);
    if(!found->seed_prediction_robustness || prediction_difference>*found->seed_prediction_robustness)
    {
        found->seed_prediction_robustness=prediction_difference;
        found->prediction_endpoint=identity;
    }
    if(!found->seed_log_width_difference_inf || width_difference>*found->seed_log_width_difference_inf)
    {
        found->seed_log_width_difference_inf=width_difference;
        found->log_width_endpoint=identity;
    }
    if(!found->seed_scaled_ac_difference || ac_difference>*found->seed_scaled_ac_difference)
    {
        found->seed_scaled_ac_difference=ac_difference;
        found->scaled_ac_endpoint=identity;
    }
}

void SplitDiagnosticsAccumulator::FinishDirection(const std::string & direction,bool complete,const std::string & reason)
{
    if(!inputs_valid_) return;
    auto found=std::find_if(result_.directions.begin(),result_.directions.end(),[&](const auto & item) {
        return item.direction==direction;
    });
    if(found==result_.directions.end())
    {
        result_.directions.push_back({});
        found=std::prev(result_.directions.end());
        found->direction=direction;
    }
    if(complete)
    {
        found->status="Assessed";
        ++result_.assessed_direction_count;
    }
    else
    {
        found->status="Incomplete";
        found->reason=reason.empty() ? "direction-incomplete" : reason;
        ++result_.incomplete_direction_count;
    }
}

LocalSplitDiagnostics SplitDiagnosticsAccumulator::Finish() const
{
    auto result=result_;
    if(!null_qualified_)
    {
        result.status="Unavailable";
        if(result.reason.empty()) result.reason="null-reference-not-qualified";
        return result;
    }
    if(!inputs_valid_)
    {
        result.status="Unavailable";
        if(result.reason.empty()) result.reason="invalid-or-nonfinite-split-input";
        return result;
    }
    for(const auto & direction:result.directions)
    {
        auto & aggregate=result;
        if(direction.seed_prediction_robustness && (!aggregate.seed_prediction_robustness ||
            *direction.seed_prediction_robustness>*aggregate.seed_prediction_robustness))
        {
            aggregate.seed_prediction_robustness=direction.seed_prediction_robustness;
            aggregate.seed_prediction_endpoint=direction.prediction_endpoint;
        }
        if(direction.seed_log_width_difference_inf && (!aggregate.seed_log_width_difference_inf ||
            *direction.seed_log_width_difference_inf>*aggregate.seed_log_width_difference_inf))
        {
            aggregate.seed_log_width_difference_inf=direction.seed_log_width_difference_inf;
            aggregate.seed_log_width_endpoint=direction.log_width_endpoint;
        }
        if(direction.seed_scaled_ac_difference && (!aggregate.seed_scaled_ac_difference ||
            *direction.seed_scaled_ac_difference>*aggregate.seed_scaled_ac_difference))
        {
            aggregate.seed_scaled_ac_difference=direction.seed_scaled_ac_difference;
            aggregate.seed_scaled_ac_endpoint=direction.scaled_ac_endpoint;
        }
    }
    if(result.incomplete_direction_count || result.assessed_direction_count<expected_directions_)
    {
        result.status="Incomplete";
        if(result.reason.empty()) result.reason="incomplete-perturbation-directions";
    }
    else
    {
        result.status="Qualified";
        result.reason.clear();
    }
    return result;
}

} // namespace rhbm_gem::core::joint_component
