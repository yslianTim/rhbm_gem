#include <gtest/gtest.h>
#include "support/SecondStageTestSupport.hpp"
#include "support/SecondStageNumericalProbe.hpp"
#include "core/detail/second_stage/CandidateEvaluation.hpp"
#include "core/detail/second_stage/CandidateTransaction.hpp"
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
    detail::ObserveSelectionAudit(&session,false,false,"unavailable","previous-objective-unavailable");
    detail::ObserveSelectionAudit(&session,true,true,"empty_after_salvage","no-selection-remains",2);
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
        auto result=detail::EvaluateCandidate(overlay,detail::GlobalCandidateReference{
            fixture.sample_ref_list,domain,&best,&*previous,counters,observe ? &session : nullptr,false});
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
    if (auto * record=trial.Record()) {
        record->outcome="rejected"; record->category=detail::AuditCategory::Guard;
        record->reason="suspicious"; record->first_atom=3; record->atom_count=1;
    }
    trial.BeginBoundary(detail::BoundaryAcceptancePolicy::Ordinary,detail::BoundaryObservationStage::Backtracking,0.25);
    if (auto * record=trial.Record()) {
        EXPECT_EQ(record->first_atom,2U); EXPECT_EQ(record->atom_count,2U);
    }
    trial.Flush();
    if (!session.Enabled()) return;
    const auto & batch=session.Audit()->batch;
    EXPECT_EQ(batch.stages[static_cast<std::size_t>(detail::AuditStage::BoundaryEndpoint)].accepted,1U);
    EXPECT_EQ(batch.stages[static_cast<std::size_t>(detail::AuditStage::BoundaryCorrection)].rejected,1U);
    EXPECT_EQ(batch.stages[static_cast<std::size_t>(detail::AuditStage::BoundaryBacktracking)].accepted,1U);
    ASSERT_EQ(batch.detail_count,1U);
    EXPECT_EQ(batch.details[0].stage,detail::AuditStage::BoundaryCorrection);
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
        .component=component,.previous_audit=&*previous};
    std::array<std::size_t,4> work{}; double objective=0;
    for (int mode=0;mode<3;++mode) {
        detail::SecondStageObservationSession session(mode==1,Write);
        detail::JointCandidateObservation trial(mode ? &session : nullptr,key);
        trial.BeginBoundary(reference.policy,detail::BoundaryObservationStage::Endpoint,1.0);
        BeginNumericalCapture();
        const auto candidate=detail::EvaluateCandidate(overlay,reference,mode ? &trial : nullptr);
        const auto capture=EndNumericalCapture();
        ASSERT_TRUE(candidate);
        if (mode==0) { work=capture.work; objective=candidate->audit_objective.GetTotalObjective(); }
        else { EXPECT_EQ(capture.work,work); EXPECT_DOUBLE_EQ(candidate->audit_objective.GetTotalObjective(),objective); }
    }
}

TEST(SecondStageObservationTest, FinalPolishAcceptanceDoesNotImplyCertificationOrApplication)
{
    LogScope logs; written.clear();
    detail::SecondStageObservationSession session(false,Write);
    if (!session.Enabled()) return;
    auto & data=*session.Audit();
    data.polish_attempted=true; data.polish_accepted=true;
    data.polish_status=detail::FinalPolishResidualSafetyStatus::Failed;
    data.polish_certificate.emplace();
    data.polish_certificate->certificate.operator_complete=false;
    detail::LogDecisionAuditTerminal(session,"converged","latest-validated",{},{});
    ASSERT_EQ(written.size(),1U);
    EXPECT_NE(written[0].find("\"objective_accepted\":true,\"operator_certified\":false"),std::string::npos);
    EXPECT_NE(written[0].find("\"applied\":false"),std::string::npos);
    EXPECT_EQ(written[0].find("accepted_active_p99"),std::string::npos);
    EXPECT_NE(written[0].find("\"value\":null,\"reason\":\"unavailable\""),std::string::npos);
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
        builder.ReconcileSelectedBoundaries(inputs);
        EXPECT_EQ(builder.View().accepted_key_list.empty(),mode==2 || mode==3);
        if (!session.Enabled()) continue;
        const auto & ordinary=session.Audit()->selection[0];
        EXPECT_EQ(ordinary.executed,mode==1 || mode==3);
        EXPECT_EQ(ordinary.result,mode==0 ? "skipped" : mode==1 ? "passed" : mode==2 ? "unavailable" : "empty_after_salvage");
        if (mode==3) EXPECT_EQ(ordinary.removed_clusters,2U);
        EXPECT_FALSE(session.Audit()->selection[1].executed);
    }
}
