#include <gtest/gtest.h>
#include "support/SecondStageTestSupport.hpp"
#include "support/SecondStageNumericalProbe.hpp"
#include "core/detail/second_stage/CandidateEvaluation.hpp"
#include "core/detail/second_stage/CandidateTransaction.hpp"
#include "core/detail/second_stage/DependencyPolish.hpp"
#include "core/detail/second_stage/IterationProcess.hpp"
#include "core/detail/second_stage/IterationResult.hpp"
#include "core/detail/second_stage/observation/SecondStageLogging.hpp"
#include "core/detail/second_stage/observation/SecondStageObservation.hpp"
#include "core/detail/second_stage/observation/PerformanceCounters.hpp"
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <future>
#include <limits>
#include <sstream>
#include <type_traits>

namespace {
namespace detail=rhbm_gem::core::detail;
namespace rg=rhbm_gem;
using namespace second_stage_test;
struct LogScope { LogLevel previous{Logger::GetLogLevel()}; LogScope(){Logger::SetLogLevel(LogLevel::Debug);} ~LogScope(){Logger::SetLogLevel(previous);SetAuditFault(AuditFault::None);} };
std::vector<std::string> written;
void Write(std::string_view text) { written.emplace_back(text); }
void FailWrite(std::string_view) { throw std::runtime_error("writer unavailable"); }
constexpr bool CompiledAudit() {
#ifdef RHBM_GEM_ENABLE_SECOND_STAGE_AUDIT
    return true;
#else
    return false;
#endif
}
}

TEST(SecondStageObservationTest, BoundedWorkerMergeKeepsDeterministicDetailsAndAllCounts)
{
    static_assert(std::is_trivially_copyable_v<detail::AuditEvent>);
    std::array<detail::AuditBatch,4> batches;
    std::vector<std::future<void>> workers;
    for(std::size_t worker=0;worker<4;++worker)
        workers.push_back(std::async(std::launch::async,[&,worker] {
            for(std::size_t i=100;i-- > 0;)
                batches[worker].Add({.stage=detail::AuditStage::LocalSearch, .category=detail::AuditCategory::Rejected,
                    .first_atom=worker, .atom_count=1, .trial=i, .outcome="rejected", .reason="previous-gate"});
        }));
    for(auto & worker:workers) worker.get();
    detail::AuditBatch forward, reverse;
    for(const auto & batch:batches) { EXPECT_EQ(batch.detail_count,5U); forward.Merge(batch); }
    for(std::size_t i=batches.size();i-- > 0;) reverse.Merge(batches[i]);
    EXPECT_EQ(forward.abnormal_count,400U); EXPECT_EQ(forward.detail_count,5U);
    EXPECT_EQ(forward.categories[static_cast<std::size_t>(detail::AuditCategory::Rejected)],400U);
    for(std::size_t i=0;i<5;++i) {
        EXPECT_EQ(forward.details[i].first_atom,0U); EXPECT_EQ(forward.details[i].trial,i);
        EXPECT_EQ(forward.details[i].trial,reverse.details[i].trial);
    }
}

TEST(SecondStageObservationTest, SingleSwitchQuietAndVerbosityGateAllocationAndOutput)
{
    LogScope logs;
    for(auto level:{LogLevel::Info,LogLevel::Debug})
    for(bool quiet:{false,true})
    {
        Logger::SetLogLevel(level); written.clear();
        detail::SecondStageObservationSession session(quiet,Write);
        EXPECT_EQ(session.Enabled(),CompiledAudit() && !quiet && level==LogLevel::Debug);
        auto options=MakeSecondStageOptions(); detail::LogDecisionAuditStart(session,options);
        EXPECT_EQ(written.size(),session.Enabled() ? 1U : 0U);
    }
}

TEST(SecondStageObservationTest, AllocationCollectionAndWriterFailuresDisableTheSession)
{
    LogScope logs;
    SetAuditFault(AuditFault::Allocation);
    detail::SecondStageObservationSession failed_allocation(false,Write);
    EXPECT_FALSE(failed_allocation.Enabled());
    SetAuditFault(AuditFault::None);
    detail::SecondStageObservationSession failed_collection(false,Write);
    SetAuditFault(AuditFault::Collection);
    failed_collection.Record({.category=detail::AuditCategory::Rejected});
    EXPECT_FALSE(failed_collection.Enabled());
    SetAuditFault(AuditFault::None);
    detail::SecondStageObservationSession failed_writer(false,FailWrite);
    EXPECT_NO_THROW(detail::LogDecisionAuditStart(failed_writer,MakeSecondStageOptions()));
    EXPECT_FALSE(failed_writer.Enabled());
}

TEST(SecondStageObservationTest, SelectionExecutionAndMissingCertificateRemainExplicit)
{
    LogScope logs; written.clear();
    detail::SecondStageObservationSession session(false,Write);
    session.BeginAttempt(1,1,1,false,false);
    session.ObserveSelectionAudit(false,false,"unavailable","previous-objective-unavailable");
    session.ObserveSelectionAudit(true,true,"empty_after_salvage","no-selection-remains",2);
    detail::IterationResult result; result.attempt_number=1;
    result.stop_reason=detail::SecondStageStopReason::AllRejectedBacktrackingExhausted;
    detail::LogDecisionAuditIteration(session,result);
    if(!CompiledAudit()) { EXPECT_TRUE(written.empty()); return; }
    ASSERT_EQ(written.size(),1U);
    EXPECT_NE(written[0].find("\"ordinary\":{\"executed\":false"),std::string::npos);
    EXPECT_NE(written[0].find("\"after_rescue\":{\"executed\":true"),std::string::npos);
    EXPECT_NE(written[0].find("\"result\":\"empty_after_salvage\""),std::string::npos);
    EXPECT_NE(written[0].find("\"reference\":\"iteration_previous\",\"status\":\"not_evaluated\""),std::string::npos);
}

TEST(SecondStageObservationTest, GlobalGateRecordsActualReferencesAndDoesNotAddEvaluation)
{
    LogScope logs;
    auto fixture=BuildJointPolishFixture({{6.0,0.5,0.0}},{{6.2,0.5,0.0}});
    const detail::ClusterKey key{0}; const auto baseline=detail::BuildResidualBaseline(fixture.context,fixture.state);
    const auto domain=detail::BuildObjectiveDomain(fixture.context,baseline.model_snapshot,{key});
    const auto previous=detail::EvaluateAuditObjective(domain,baseline); ASSERT_TRUE(previous);
    auto improved=fixture.state; improved[0]=MakeGaussianResult({6.1,0.5,0.0});
    const auto patch=detail::FitStatePatch::FromState(improved,key);
    const detail::CandidateEvaluationOverlay overlay{fixture.context,baseline,fixture.state,patch};
    detail::ClusterSolverWorkspaceMap workspaces; detail::BoundaryJointCorrectionWorkspaceMap corrections;
    detail::PerformanceCounters counters(true,fixture.context,workspaces,corrections);
    detail::SecondStageObservationSession session(false,Write);
    session.BeginAttempt(1,1,1,false,false);
    const detail::ObjectiveBreakdown best{0,0,0};
    std::array<std::size_t,4> unobserved_work{};
    for(bool observe:{false,true})
    {
        BeginNumericalCapture();
        auto result=detail::EvaluateGlobalCandidate(overlay,detail::GlobalCandidateReference{
            fixture.sample_ref_list,domain,&best,&*previous,counters},observe ? &session : nullptr,false);
        const auto captured=EndNumericalCapture();
        EXPECT_FALSE(result);
        if(!observe) unobserved_work=captured.work; else EXPECT_EQ(captured.work,unobserved_work);
    }
    if(session.Enabled()) {
        ASSERT_TRUE(session.Audit()->selection[0].candidate);
        EXPECT_EQ(session.Audit()->selection[0].evaluations,1U);
        EXPECT_DOUBLE_EQ(session.Audit()->selection[0].best->GetTotalObjective(),0.0);
        EXPECT_EQ(session.Audit()->batch.detail_count,1U);
    }
}

TEST(SecondStageObservationTest, CollectionAndWriterFailurePreserveSmallProductionRun)
{
    LogScope logs;
    std::vector<double> baseline_values; NumericalCapture baseline_capture;
    for(auto level:{LogLevel::Debug,LogLevel::Info})
    for(auto fault:{AuditFault::None,AuditFault::Allocation,AuditFault::Collection,AuditFault::Writer})
    {
        Logger::SetLogLevel(level);
        auto model=BuildJointPolishDefenseModel();
        model->EditAnalysis().CopyLocalFittingStageResult(rg::FittingStage::Second,rg::FittingStage::First);
        auto options=MakeSecondStageOptions(); options.quiet_mode=false;
        SetAuditFault(fault); testing::internal::CaptureStdout(); BeginNumericalCapture();
        detail::RunSecondStageIterations(*model,options);
        const auto capture=EndNumericalCapture(); testing::internal::GetCapturedStdout();
        std::vector<double> values;
        for(auto * atom:model->GetSelectedAtoms()) {
            const auto view=rg::AtomLocalPotentialView::For(*atom); const auto fit=view.GetGaussianResult(rg::FittingStage::Second);
            for(int par=0;par<3;++par) values.push_back(fit.mdpde.GetModelParameter(par));
            for(const auto & sample:view.GetPeelingSamplingEntries(false)) values.push_back(sample.response);
        }
        if(baseline_values.empty()) { baseline_values=values; baseline_capture=capture; }
        else { EXPECT_EQ(values,baseline_values); EXPECT_EQ(capture.commits,baseline_capture.commits); EXPECT_EQ(capture.work,baseline_capture.work); EXPECT_EQ(capture.terminal,baseline_capture.terminal); EXPECT_EQ(capture.backgrounds,baseline_capture.backgrounds); }
    }
}

TEST(SecondStageObservationTest, BasicPerformanceSummaryRemainsAvailableWithoutAudit)
{
    LogScope logs;
    for(bool quiet:{true,false}) {
        detail::SecondStageContext context; detail::ClusterSolverWorkspaceMap workspaces;
        detail::BoundaryJointCorrectionWorkspaceMap corrections;
        testing::internal::CaptureStdout();
        { detail::PerformanceCounters counters(quiet,context,workspaces,corrections);
          counters.RecordBoundaryReconciliation(4.25); counters.RecordBoundaryJointCorrection(3.0);
          counters.RecordDependencyPolish(8.5); }
        const auto output=testing::internal::GetCapturedStdout();
        if(quiet) EXPECT_TRUE(output.empty());
        else { EXPECT_NE(output.find("boundary_reconciliation_ms = 4.250"),std::string::npos);
          EXPECT_NE(output.find("dependency_polish_ms = 8.500"),std::string::npos); }
        EXPECT_EQ(output.find("full_state_materializations"),std::string::npos);
    }
}

TEST(SecondStageObservationTest, LaterCorrectionDoesNotRelabelSelectedEndpoint)
{
    LogScope logs;
    detail::SecondStageObservationSession session(false,Write);
    session.BeginAttempt(1,1,1,false,false);
    const detail::ClusterKey key{2,3};
    detail::JointCandidateObservation trial(&session,key);
    trial.BeginBoundary(detail::BoundaryAcceptancePolicy::Ordinary,detail::BoundaryObservationStage::Endpoint,1.0);
    trial.BeginBoundary(detail::BoundaryAcceptancePolicy::Ordinary,detail::BoundaryObservationStage::Correction,0.5);
    trial.RejectMemberSamplesUnavailable({3});
    trial.BeginBoundary(detail::BoundaryAcceptancePolicy::Ordinary,detail::BoundaryObservationStage::Backtracking,0.25);
    trial.RejectSuspicious();
    trial.Flush();
    if (!session.Enabled()) return;
    const auto & batch=session.Audit()->batch;
    EXPECT_EQ(batch.stages[static_cast<std::size_t>(detail::AuditStage::BoundaryEndpoint)].accepted,1U);
    EXPECT_EQ(batch.stages[static_cast<std::size_t>(detail::AuditStage::BoundaryCorrection)].rejected,1U);
    EXPECT_EQ(batch.stages[static_cast<std::size_t>(detail::AuditStage::BoundaryBacktracking)].rejected,1U);
    ASSERT_EQ(batch.detail_count,2U);
    EXPECT_EQ(batch.details[0].stage,detail::AuditStage::BoundaryCorrection);
    EXPECT_EQ(batch.details[0].first_atom,3U); EXPECT_EQ(batch.details[0].atom_count,1U);
    EXPECT_EQ(batch.details[1].first_atom,2U); EXPECT_EQ(batch.details[1].atom_count,2U);
}

TEST(SecondStageObservationTest, BoundaryMissingQuietAndEnabledObserverPreserveDecisionAndWork)
{
    LogScope logs;
    auto fixture=BuildJointPolishFixture({{6.0,0.5,0.0}},{{6.4,0.5,0.0}});
    const detail::ClusterKey key{0};
    detail::CouplingGraphPartition partition; partition.sample_id_list_by_key[key]=fixture.sample_ref_list;
    const auto baseline=detail::BuildResidualBaseline(fixture.context,fixture.state);
    const auto domain=detail::BuildObjectiveDomain(fixture.context,baseline.model_snapshot,{key});
    const auto previous_by_key=detail::BuildObjectiveByKey(partition,domain,baseline);
    const auto previous=detail::EvaluateAuditObjective(domain,baseline); ASSERT_TRUE(previous);
    const detail::BoundaryReconciliationComponent component{.key_list={key},
        .affected_sample_ref_list=fixture.sample_ref_list,.halo_atom_index_list=key};
    const detail::FitState improved{MakeGaussianResult({6.2,0.5,0.0})};
    const auto patch=detail::FitStatePatch::FromState(improved,key);
    const detail::CandidateEvaluationOverlay overlay{fixture.context,baseline,fixture.state,patch};
    detail::ClusterSolverWorkspaceMap workspaces; detail::BoundaryJointCorrectionWorkspaceMap corrections;
    detail::PerformanceCounters counters(true,fixture.context,workspaces,corrections);
    const detail::BoundaryCandidateReference reference{
        .policy=detail::BoundaryAcceptancePolicy::Ordinary,.samples_by_key=partition.sample_id_list_by_key,
        .domain=domain,.previous_objective_by_key=previous_by_key,.best_audit=nullptr,.counters=counters,
        .component=component};
    std::array<std::size_t,4> work{}; double objective=0;
    for (int mode=0;mode<3;++mode) {
        detail::SecondStageObservationSession session(mode==1,Write);
        detail::JointCandidateObservation trial(mode ? &session : nullptr,key);
        trial.BeginBoundary(reference.policy,detail::BoundaryObservationStage::Endpoint,1.0);
        BeginNumericalCapture();
        const auto candidate=detail::EvaluateBoundaryCandidate(overlay,reference,&*previous,mode ? &trial : nullptr);
        const auto capture=EndNumericalCapture();
        ASSERT_TRUE(candidate);
        if (mode==0) { work=capture.work; objective=candidate->GetTotalObjective(); }
        else { EXPECT_EQ(capture.work,work); EXPECT_DOUBLE_EQ(candidate->GetTotalObjective(),objective); }
    }
}

TEST(SecondStageObservationTest, FinalPolishAcceptanceDoesNotImplyCertificationOrApplication)
{
    LogScope logs; written.clear();
    detail::SecondStageObservationSession session(false,Write);
    if (!session.Enabled()) return;
    session.ObserveFinalPolishAttempt();
    detail::FinalDependencyPolishResult result; result.accepted=true;
    detail::ConvergenceAssessment certificate; certificate.certificate.operator_complete=false;
    session.ObserveFinalCertification(result,detail::FinalPolishResidualSafetyStatus::Failed,certificate,false);
    detail::LogDecisionAuditTerminal(session,"converged","latest-validated",{});
    ASSERT_EQ(written.size(),1U);
    EXPECT_NE(written[0].find("\"objective_accepted\":true,\"operator_certified\":false"),std::string::npos);
    EXPECT_NE(written[0].find("\"applied\":false"),std::string::npos);
    EXPECT_EQ(written[0].find("accepted_active_p99"),std::string::npos);
    EXPECT_NE(written[0].find("\"value\":null,\"reason\":\"unavailable\""),std::string::npos);
}

TEST(SecondStageObservationTest, CorrectionKeepsGuardShortCircuitAndStrictImprovementReference)
{
    LogScope logs;
    auto fixture=BuildJointPolishFixture({{6.0,0.5,0.0}},{{6.4,0.5,0.0}});
    const detail::ClusterKey key{0};
    detail::CouplingGraphPartition partition; partition.sample_id_list_by_key[key]=fixture.sample_ref_list;
    const auto baseline=detail::BuildResidualBaseline(fixture.context,fixture.state);
    const auto domain=detail::BuildObjectiveDomain(fixture.context,baseline.model_snapshot,{key});
    const auto previous_by_key=detail::BuildObjectiveByKey(partition,domain,baseline);
    const auto previous=detail::EvaluateAuditObjective(domain,baseline); ASSERT_TRUE(previous);
    const detail::BoundaryReconciliationComponent component{.key_list={key},
        .affected_sample_ref_list=fixture.sample_ref_list,.halo_atom_index_list=key};
    const detail::FitStatePatch empty_patch;
    const detail::FitStateView endpoint{fixture.state,empty_patch};
    detail::ClusterSolverWorkspaceMap workspaces; detail::BoundaryJointCorrectionWorkspaceMap corrections;
    detail::PerformanceCounters counters(true,fixture.context,workspaces,corrections);
    const detail::BoundaryCandidateReference reference{
        .policy=detail::BoundaryAcceptancePolicy::Ordinary,.samples_by_key=partition.sample_id_list_by_key,
        .domain=domain,.previous_objective_by_key=previous_by_key,.best_audit=nullptr,.counters=counters,.component=component};
    for (bool suspicious : {false,true})
    {
        const detail::FitState state{MakeGaussianResult({6.2,suspicious ? 100.0 : 0.5,0.0})};
        const auto patch=detail::FitStatePatch::FromState(state,key);
        const detail::CandidateEvaluationOverlay overlay{fixture.context,baseline,fixture.state,patch};
        const auto improvement=detail::EvaluateObjectiveDelta(overlay,fixture.sample_ref_list,domain,*previous,counters);
        ASSERT_TRUE(improvement);
        std::array<std::size_t,4> unobserved_work{};
        for (int mode=0;mode<3;++mode)
        {
            detail::SecondStageObservationSession session(mode==1,Write);
            detail::JointCandidateObservation trial(mode ? &session : nullptr,key);
            trial.BeginBoundary(reference.policy,detail::BoundaryObservationStage::Correction,0.5);
            BeginNumericalCapture();
            EXPECT_FALSE(detail::EvaluateBoundaryCorrection(overlay,reference,endpoint,*previous,*improvement,
                mode ? &trial : nullptr));
            const auto capture=EndNumericalCapture();
            EXPECT_EQ(capture.work[2],suspicious ? 0U : 4U);
            if (mode==0) unobserved_work=capture.work; else EXPECT_EQ(capture.work,unobserved_work);
            trial.Flush();
            if (mode!=2 || !session.Enabled()) continue;
            ASSERT_EQ(session.Audit()->batch.detail_count,1U);
            const auto & event=session.Audit()->batch.details[0];
            EXPECT_EQ(event.stage,detail::AuditStage::BoundaryCorrection);
            EXPECT_EQ(event.trial,0U); EXPECT_EQ(event.factor,0.5);
            EXPECT_EQ(event.reason,suspicious ? "suspicious" : "strict-improvement");
            EXPECT_EQ(event.reference,suspicious ? "iteration_previous" : "best_boundary_candidate");
            EXPECT_EQ(event.previous_checked,!suspicious); EXPECT_FALSE(event.best_checked);
            if (!suspicious) EXPECT_DOUBLE_EQ(event.previous->GetTotalObjective(),improvement->GetTotalObjective());
        }
    }
}

TEST(SecondStageObservationTest, PerformanceRecordingStopsWithItsSession)
{
    LogScope logs;
    detail::SecondStageContext context; detail::ClusterSolverWorkspaceMap workspaces;
    detail::BoundaryJointCorrectionWorkspaceMap corrections;
    detail::SecondStageObservationSession session(false,Write);
    detail::PerformanceCounters counters(true,context,workspaces,corrections,&session);
    counters.RecordFullStateMaterialization();
    const auto before=counters.AuditCounts();
    session.Disable(); counters.RecordFullStateMaterialization();
    EXPECT_FALSE(counters.AuditEnabled());
    EXPECT_EQ(before[0],CompiledAudit() ? 1U : 0U);
    EXPECT_EQ(counters.AuditCounts(),(std::array<std::size_t,6>{}));
    detail::PerformanceCounters missing(false,context,workspaces,corrections);
    EXPECT_FALSE(missing.AuditEnabled());
}

TEST(SecondStageObservationTest, ScoresAndFinalizationKeepOriginalReferencesAndChosenState)
{
    LogScope logs;
    detail::SecondStageObservationSession session(false,Write);
    static_assert(std::is_same_v<decltype(session.Audit()), const detail::SecondStageAuditData *>);
    session.BeginAttempt(1,1,1,false,false);
    const detail::ObjectiveBreakdown previous{3.0,0,0}, candidate{2.0,0,0};
    detail::ObjectiveBreakdown best{2.5,0,0};
    session.ObserveScoreReferences(previous,&best);
    session.ObserveCandidateScoreSource(true);
    session.ObserveScores(previous,candidate,&best);
    best.fit_range_residual_objective=1.0;
    if (!session.Enabled()) { EXPECT_EQ(session.Audit(),nullptr); return; }
    EXPECT_DOUBLE_EQ(session.Audit()->best->GetTotalObjective(),2.5);
    EXPECT_EQ(session.Audit()->score_source,"selection_audit");
    session.BeginFinalization(nullptr);
    EXPECT_DOUBLE_EQ(session.Audit()->final_objective->GetTotalObjective(),2.0);
    session.BeginFinalization(&best);
    EXPECT_DOUBLE_EQ(session.Audit()->final_objective->GetTotalObjective(),1.0);

    session.BeginAttempt(2,1,1,false,false);
    session.ObserveScoreReferences(previous,&best);
    session.ObserveCandidateScoreSource(false);
    EXPECT_EQ(session.Audit()->score_source,"post_commit_evaluation");
    detail::CandidateCommitResult rejected;
    detail::TrustRegionStateSet radii;
    session.ObserveCommit({},rejected,radii);
    EXPECT_EQ(session.Audit()->score_source,"restored_previous");
    EXPECT_DOUBLE_EQ(session.Audit()->candidate->GetTotalObjective(),3.0);
    EXPECT_FALSE(session.Audit()->convergence);
    session.BeginFinalization(nullptr);
    for (const auto status : {detail::FinalPolishResidualSafetyStatus::NotEvaluated,
        detail::FinalPolishResidualSafetyStatus::Error, detail::FinalPolishResidualSafetyStatus::Failed,
        detail::FinalPolishResidualSafetyStatus::AbsolutePassed})
    {
        detail::FinalDependencyPolishResult result;
        result.accepted=status!=detail::FinalPolishResidualSafetyStatus::NotEvaluated;
        if (result.accepted) result.objective=candidate;
        std::optional<detail::ConvergenceAssessment> certificate;
        const bool applied=status==detail::FinalPolishResidualSafetyStatus::AbsolutePassed;
        if (status==detail::FinalPolishResidualSafetyStatus::Failed || applied) certificate.emplace();
        session.ObserveFinalPolishAttempt();
        session.ObserveFinalCertification(result,status,certificate,applied);
        EXPECT_EQ(session.Audit()->polish_status,status);
        EXPECT_EQ(session.Audit()->polish_applied,applied);
        EXPECT_EQ(session.Audit()->polish_certificate.has_value(),certificate.has_value());
        EXPECT_DOUBLE_EQ(session.Audit()->final_objective->GetTotalObjective(),applied ? 2.0 : 3.0);
    }
    session.Disable();
    EXPECT_NO_THROW(session.ObserveScores(previous,candidate,&best));
    EXPECT_NO_THROW(session.BeginFinalization(&best));
    EXPECT_NO_THROW(session.ObserveFinalPolishAttempt());
    EXPECT_EQ(session.Audit(),nullptr);
}

TEST(SecondStageObservationTest, SelectionAuditReportsActualConditionalAndSalvageBranches)
{
    LogScope logs;
    for (int mode=0;mode<4;++mode) {
        SCOPED_TRACE(mode);
        const std::vector<rg::GaussianModel3D> models{{6.0,0.5,0.0},{7.0,0.5,0.0}};
        auto fixture=BuildJointPolishFixture(models,models);
        const std::vector<detail::ClusterKey> keys{{0},{1}};
        detail::CouplingGraphPartition partition;
        for (const auto & key:keys) partition.sample_id_list_by_key[key]=fixture.sample_ref_list;
        if (mode!=0) partition.boundary_sample_dependency_list={{{0,0},keys,{0,1}}};
        const auto baseline=detail::BuildResidualBaseline(fixture.context,fixture.state);
        auto domain=detail::BuildObjectiveDomain(fixture.context,baseline.model_snapshot,keys);
        if (mode==2) domain.cluster_by_key.at(keys[0]).scale.reset();
        const auto previous=detail::BuildObjectiveByKey(partition,domain,baseline);
        const detail::PolishProvenance provenance(2,0);
        const detail::SuspiciousBlockActivity fixed{{1,1},{1,1},{1,1}};
        const std::vector<double> ridge(2,1.0); const detail::ClusterHealthMap health;
        detail::BestAuditState best;
        if (mode==3) best=detail::AuditedState{detail::ObjectiveBreakdown{-1.0,0,0},fixture.state,false,0};
        detail::TrustRegionStateSet trust; trust.Reconcile(keys);
        detail::ClusterSolverWorkspaceMap workspaces; detail::BoundaryJointCorrectionWorkspaceMap corrections;
        detail::SecondStageObservationSession session(false,Write);
        session.BeginAttempt(1,1,1,false,false);
        detail::PerformanceCounters counters(true,fixture.context,workspaces,corrections);
        const auto options=MakeSecondStageOptions();
        const detail::CandidateSelectionInputs inputs{
            fixture.context,options,baseline,partition,health,fixture.state,provenance,fixture.state,fixed,
            ridge,domain,previous,best,trust,workspaces,corrections,counters,&session};
        detail::CandidateSelection selected; selected.block_activity=fixed;
        selected.assembled_state=fixture.state; selected.assembled_polish_provenance=provenance;
        selected.accepted_key_list=keys;
        detail::CandidateTransactionBuilder builder(std::move(selected));
        BeginNumericalCapture();
        builder.ReconcileSelectedBoundaries(inputs);
        const auto capture=EndNumericalCapture();
        EXPECT_EQ(builder.View().accepted_key_list.empty(),mode==2 || mode==3);
        if (mode==1 || mode==3) {
            // Previous audit + two member objectives + boundary and selection
            // deltas. Each delta includes two contribution evaluations. Only a
            // failed selection adds the ranking delta; empty selection adds none.
            EXPECT_EQ(capture.work[static_cast<std::size_t>(Work::Objective)],mode==1 ? 9U : 12U);
            EXPECT_EQ(builder.View().final_audit_objective.has_value(),mode==1);
        }
        if (!session.Enabled()) continue;
        const auto & ordinary=session.Audit()->selection[0];
        EXPECT_EQ(ordinary.executed,mode==1 || mode==3);
        EXPECT_EQ(ordinary.result,mode==0 ? "skipped" : mode==1 ? "passed" : mode==2 ? "unavailable" : "empty_after_salvage");
        if (mode==3) EXPECT_EQ(ordinary.removed_clusters,2U);
        EXPECT_FALSE(session.Audit()->selection[1].executed);
    }
}
