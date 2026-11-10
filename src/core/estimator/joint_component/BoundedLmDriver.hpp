// Portions adapted from Eigen 5.0.1 NonLinearOptimization/LevenbergMarquardt.h.
// Copyright (C) 2009 Thomas Capricelli <orzel@freehackers.org>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once
#include <Eigen/QR>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>

namespace rhbm_gem::core::joint_component {

struct BoundedLmOptions
{
    double initial_radius_factor{100.};
    Eigen::Index evaluation_budget{400};
    double ftol{std::sqrt(std::numeric_limits<double>::epsilon())};
    double xtol{std::sqrt(std::numeric_limits<double>::epsilon())};
    double gtol{};
};

enum class BoundedLmStatus
{
    Ready, Running, InvalidInput, ReductionTolerance, StepTolerance, ReductionAndStepTolerance,
    Stationary, EvaluationBudget, ReductionPrecisionLimit, StepPrecisionLimit, GradientPrecisionLimit,
    ProblemStopped
};

struct BoundedLmState
{
    Eigen::VectorXd diagonal,column_norms,reduced_response;
    Eigen::VectorXd trial_eta,model_step;
    Eigen::MatrixXd reduced_jacobian;
    double accepted_residual_norm{},scaled_eta_norm{},scaled_gradient_norm{},trust_radius{},damping{};
    Eigen::Index accepted_iterations{},function_evaluations{};
};

// The problem supplies full residual norms, a reduced linear model, a feasible damped
// step, and stationarity and trial checks. AssessTrial records the candidate and
// returns its trustworthiness. The caller observes accepted updates after Advance.
template<class Problem>
class BoundedLmDriver
{
public:
    explicit BoundedLmDriver(Problem & problem,BoundedLmOptions options={})
        :problem_(problem),options_(options) {}

    Eigen::Index AcceptedIterations() const {return state_.accepted_iterations;}
    // Includes the initial evaluation and evaluated trials, excluding failed trials.
    Eigen::Index FunctionEvaluations() const {return state_.function_evaluations;}

    BoundedLmStatus Initialize(Eigen::VectorXd & eta)
    {
        state_.function_evaluations=state_.accepted_iterations=0;
        state_.damping=0.;
        const auto rows=problem_.ResidualCount();
        if(eta.size()<=0 || rows<eta.size() || options_.ftol<0. || options_.xtol<0. ||
            options_.gtol<0. || options_.evaluation_budget<=0 || options_.initial_radius_factor<=0.)
            return BoundedLmStatus::InvalidInput;

        state_.function_evaluations=1;
        double residual_norm{};
        if(!problem_.Evaluate(eta,residual_norm)) return BoundedLmStatus::ProblemStopped;
        state_.accepted_residual_norm=residual_norm;
        return BoundedLmStatus::Ready;
    }

    BoundedLmStatus Advance(Eigen::VectorXd & eta)
    {
        auto & s=state_;
        if(!problem_.Linearize(eta,s.reduced_jacobian,s.reduced_response,s.column_norms))
            return BoundedLmStatus::ProblemStopped;

        // Preserve the reduced model for the bounded solve. The pivoted QR is
        // shared by the scaled gradient and the prediction for each retry.
        const Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(s.reduced_jacobian);
        ComputeScalingAndGradient(eta,qr);
        if(s.scaled_gradient_norm<=options_.gtol && problem_.IsStationary(eta))
            return BoundedLmStatus::Stationary;
        s.diagonal=s.diagonal.cwiseMax(s.column_norms);

        bool rejected_untrusted=false;
        for(;;)
        {
            if(!std::isfinite(s.damping) || s.damping<0. ||
                !std::isfinite(s.trust_radius) || !(s.trust_radius>0.))
                return Stop("invalid-bounded-lm-trust-state");

            const auto bounded=problem_.SolveBoundedStep(s.reduced_jacobian,s.reduced_response,eta,
                s.diagonal,s.trust_radius,s.damping);
            if(!bounded.valid) return Stop("bounded-subproblem:"+bounded.status);
            const auto & step=bounded.step;
            s.damping=bounded.damping;
            if(!step.allFinite() || !std::isfinite(s.damping) || s.damping<0. ||
                !std::isfinite(bounded.scaled_step_norm))
                return Stop("non-finite-bounded-step");

            // A small diagonal can make a meaningful parameter step look tiny in the trust metric.
            if(step.template lpNorm<Eigen::Infinity>()<=64*epsilon_)
            {
                if(problem_.IsStationary(eta))
                {
                    problem_.failure="bounded-stationary";
                    return BoundedLmStatus::Stationary;
                }
                return Stop("bounded-tiny-step-nonstationary");
            }

            s.trial_eta=eta+step;
            const double scaled_step_norm=s.diagonal.cwiseProduct(step).stableNorm();
            if(!s.trial_eta.allFinite() || !std::isfinite(scaled_step_norm))
                return Stop("non-finite-lm-step");
            if(s.accepted_iterations==0) s.trust_radius=std::min(s.trust_radius,scaled_step_norm);
            if((s.trial_eta.array()==eta.array()).all()) return Stop("unrepresentable-step");

            Reduction reduction{};
            double trial_residual_norm{};
            const bool evaluated=problem_.Evaluate(s.trial_eta,trial_residual_norm);
            if(evaluated)
            {
                ++s.function_evaluations;
                if(!std::isfinite(trial_residual_norm)) return Stop("non-finite-trial-residual");

                // Actual reduction compares full residuals. Prediction uses the
                // compact model and the damping associated with the feasible step.
                reduction.actual=-1.;
                if(.1*trial_residual_norm<s.accepted_residual_norm)
                    reduction.actual=1.-Eigen::numext::abs2(trial_residual_norm/s.accepted_residual_norm);

                s.model_step.noalias()=qr.matrixQR().template triangularView<Eigen::Upper>()*
                    (qr.colsPermutation().inverse()*step);
                const double model_term=Eigen::numext::abs2(s.model_step.stableNorm()/s.accepted_residual_norm);
                const double damping_term=Eigen::numext::abs2(std::sqrt(s.damping)*scaled_step_norm/s.accepted_residual_norm);
                reduction.predicted=model_term+damping_term/.5;
                reduction.directional_derivative=-(model_term+damping_term);
                if(!std::isfinite(reduction.predicted) || !(reduction.predicted>0.) ||
                    !std::isfinite(reduction.directional_derivative))
                    return Stop("invalid-bounded-predicted-reduction");

                reduction.ratio=reduction.actual/reduction.predicted;
                if(!std::isfinite(reduction.ratio)) return Stop("non-finite-trust-ratio");
            }

            const bool proposed=reduction.ratio>=acceptance_threshold_;
            const bool trusted=problem_.AssessTrial(eta,step,s.diagonal,s.trust_radius,s.damping,
                reduction.actual,reduction.predicted,reduction.ratio,proposed);
            if(!evaluated || !trusted)
            {
                if(!problem_.CanRetry()) return BoundedLmStatus::ProblemStopped;
                rejected_untrusted=true;
                s.trust_radius*=.25;
                s.damping*=4.;
                continue;
            }

            UpdateTrustRegion(reduction,trial_residual_norm,scaled_step_norm);
            if(proposed)
            {
                eta=s.trial_eta;
                s.scaled_eta_norm=s.diagonal.cwiseProduct(eta).stableNorm();
                s.accepted_residual_norm=trial_residual_norm;
                ++s.accepted_iterations;
                rejected_untrusted=false;
            }

            const auto status=CheckTermination(reduction,rejected_untrusted);
            if(proposed || status!=BoundedLmStatus::Running) return status;
        }
    }

private:
    static constexpr double epsilon_=std::numeric_limits<double>::epsilon();
    static constexpr double acceptance_threshold_=1e-4;

    struct Reduction
    {
        double actual{},predicted{},directional_derivative{},ratio{};
    };

    void ComputeScalingAndGradient(const Eigen::VectorXd & eta,
        const Eigen::ColPivHouseholderQR<Eigen::MatrixXd> & qr)
    {
        auto & s=state_;
        if(s.accepted_iterations==0)
        {
            s.diagonal.resize(eta.size());
            for(Eigen::Index j=0;j<eta.size();++j)
                s.diagonal[j]=s.column_norms[j]==0. ? 1. : s.column_norms[j];
            s.scaled_eta_norm=s.diagonal.cwiseProduct(eta).stableNorm();
            // The starting radius measures a log-width increment from any origin.
            s.trust_radius=options_.initial_radius_factor*s.diagonal.stableNorm();
        }

        Eigen::VectorXd projected=s.reduced_response;
        projected.applyOnTheLeft(qr.householderQ().adjoint());
        s.scaled_gradient_norm=0.;
        if(s.accepted_residual_norm==0.) return;
        for(Eigen::Index j=0;j<eta.size();++j)
        {
            const double column_norm=s.column_norms[qr.colsPermutation().indices()[j]];
            if(column_norm==0.) continue;
            const double cosine=qr.matrixQR().col(j).head(j+1).dot(projected.head(j+1)/s.accepted_residual_norm);
            s.scaled_gradient_norm=std::max(s.scaled_gradient_norm,std::abs(cosine/column_norm));
        }
    }

    void UpdateTrustRegion(const Reduction & reduction,double trial_residual_norm,double scaled_step_norm)
    {
        auto & s=state_;
        if(reduction.ratio<=.25)
        {
            double shrink=reduction.actual>=0. ? .5 :
                .5*reduction.directional_derivative/(reduction.directional_derivative+.5*reduction.actual);
            if(.1*trial_residual_norm>=s.accepted_residual_norm || shrink<.1) shrink=.1;
            s.trust_radius=shrink*std::min(s.trust_radius,scaled_step_norm/.1);
            s.damping/=shrink;
        }
        else if(!(s.damping!=0. && reduction.ratio<.75))
        {
            s.trust_radius=scaled_step_norm/.5;
            s.damping=.5*s.damping;
        }
    }

    BoundedLmStatus CheckTermination(const Reduction & reduction,bool rejected_untrusted)
    {
        const auto & s=state_;
        const bool small_reduction=std::abs(reduction.actual)<=options_.ftol && reduction.predicted<=options_.ftol;
        const bool small_radius=s.trust_radius<=options_.xtol*s.scaled_eta_norm;
        const bool budget_reached=s.function_evaluations>=options_.evaluation_budget;
        if(rejected_untrusted && (small_reduction || small_radius || budget_reached))
            return Stop("no-trustworthy-descent-step");
        if(small_reduction && .5*reduction.ratio<=1.)
            return small_radius ? BoundedLmStatus::ReductionAndStepTolerance : BoundedLmStatus::ReductionTolerance;
        if(small_radius) return BoundedLmStatus::StepTolerance;
        if(budget_reached) return BoundedLmStatus::EvaluationBudget;
        if(std::abs(reduction.actual)<=epsilon_ && reduction.predicted<=epsilon_ && .5*reduction.ratio<=1.)
            return BoundedLmStatus::ReductionPrecisionLimit;
        if(s.trust_radius<=epsilon_*s.scaled_eta_norm) return BoundedLmStatus::StepPrecisionLimit;
        if(s.scaled_gradient_norm<=epsilon_) return BoundedLmStatus::GradientPrecisionLimit;
        return BoundedLmStatus::Running;
    }

    BoundedLmStatus Stop(std::string reason)
    {
        problem_.failure=std::move(reason);
        return BoundedLmStatus::ProblemStopped;
    }

    Problem & problem_;
    BoundedLmOptions options_;
    BoundedLmState state_;
};
}
