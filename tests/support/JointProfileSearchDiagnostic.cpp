#include "support/JointProfileSearchDiagnostic.hpp"
#include <algorithm>
#include <chrono>

namespace rhbm_gem::core::joint_component {
ProfileSearchDiagnosticResult AssessProfileSearchForDiagnostics(const Domain & domain,VectorRef y,const EvaluationContext & context,SearchResult search,
    const JointProgressObserver & observer,const JointProgressComponent * progress_component)
{
    if(progress_component)
        NotifyJointProgress(observer,JointProgressPhase::CertificationStarted,*progress_component,
            search.evaluations,context.profile_budget,search.accepted,context.update_budget,search.seconds,search.stop_reason,
            false,search.accepted_objective,search.accepted_gradient_inf_norm);
    ProfileSearchDiagnosticResult out; out.search=std::move(search);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    AssessmentStageTimerForTesting total_stage("total-assessment",y.size(),out.search.eta.size());
    AssessmentStageTimerForTesting primary_stage("endpoint-primary-evaluation",y.size(),out.search.eta.size());
#endif
    const auto audit_start=std::chrono::steady_clock::now();
    const auto endpoint=EvaluateProfile(domain,y,out.search.eta,false,&context);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    primary_stage.Finish();
    AssessmentStageTimerForTesting reference_stage("reference-evaluation",y.size(),out.search.eta.size());
#endif
    const auto reference=EvaluateProfile(domain,y,out.search.eta,true,&context);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    reference_stage.Finish();
#endif
    out.assessment=AssessEvaluated(domain,y,endpoint,reference,context);
    const bool accepted=std::any_of(out.search.trials.begin(),out.search.trials.end(),[](const Trial & trial) {
        return trial.accepted && trial.trust && trial.trust->passed;
    });
    if(accepted)
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        AssessmentStageTimerForTesting trust_stage("endpoint-trust-replay",y.size(),out.search.eta.size());
#endif
        out.endpoint_trust=CheckTrust(domain,y,endpoint,context,reference);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        trust_stage.Finish();
#endif
        if(out.endpoint_trust->passed)
        {out.trusted_state=endpoint; out.trusted_assessment=out.assessment;}
        else
        {
            if(!out.search.stopped) out.search.stop_reason="endpoint-certification-failed";
            out.search.stopped=true;
            auto fallback_reference=reference;
            for(std::size_t k=out.search.trials.size();k>0;--k)
            {
                const auto & trial=out.search.trials[k-1];
                if(!trial.accepted || !trial.trust || !trial.trust->passed) continue;
                const auto & state=trial.endpoint;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                AssessmentStageTimerForTesting candidate_stage("endpoint-state-evaluation",y.size(),state.eta.size());
#endif
                const auto candidate=EvaluateState(domain,y,state.eta,state.beta,context);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                candidate_stage.Finish();
                AssessmentStageTimerForTesting replay_stage("endpoint-trust-replay",y.size(),state.eta.size());
#endif
                const bool replay_passed=CheckReplay(domain,y,candidate,context).passed;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                replay_stage.Finish();
#endif
                if(!replay_passed) continue;
                if(state.eta.size()!=fallback_reference.eta.size() ||
                    !(state.eta.array()==fallback_reference.eta.array()).all())
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                {
                    AssessmentStageTimerForTesting fallback_reference_stage("reference-evaluation",y.size(),state.eta.size());
                    fallback_reference=EvaluateProfile(domain,y,state.eta,true,&context);
                    fallback_reference_stage.Finish();
                }
#else
                    fallback_reference=EvaluateProfile(domain,y,state.eta,true,&context);
#endif
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                AssessmentStageTimerForTesting fallback_trust_stage("endpoint-trust-replay",y.size(),state.eta.size());
#endif
                const bool fallback_trust_passed=CheckTrust(domain,y,candidate,context,fallback_reference).passed;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                fallback_trust_stage.Finish();
#endif
                if(!fallback_trust_passed) continue;
                out.trusted_state=state; out.trusted_trial=k-1;
                out.trusted_assessment=AssessEvaluated(domain,y,candidate,fallback_reference,context,true);
                break;
            }
        }
    }
    out.assessment_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-audit_start).count();
    out.search_success=out.trusted_state.has_value() && !out.search.stopped; return out;
}
}
