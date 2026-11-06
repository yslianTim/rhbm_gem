#include <gtest/gtest.h>
#include "support/FirstStageSamplingFixture.hpp"
#include "support/FibonacciNeighborSelection.hpp"
#include "support/CommandTestHelpers.hpp"
#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/FittingWorkset.hpp"
#include "core/detail/PotentialFittingWorkflow.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/WidthCalibration.hpp"
#include <rhbm_gem/core/MapSampler.hpp>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>

namespace {
using namespace rhbm_gem;
using namespace rhbm_gem::core;
namespace d = rhbm_gem::core::detail;
namespace n = rhbm_gem::core::joint_component;
constexpr std::array<d::JointSeedSamplingPolicy,4> policies{{{},
    {d::JointSeedSamplingMethod::FibonacciShells,d::MapInterpolationMethod::Trilinear},
    {d::JointSeedSamplingMethod::QScoreSpiralShells,d::MapInterpolationMethod::Tricubic},
    {d::JointSeedSamplingMethod::QScoreSpiralShells,d::MapInterpolationMethod::Trilinear}}};

void CheckEndpoint(const JointProblem & problem,const JointFitResult & fit,double bound)
{
    ASSERT_EQ(fit.RuntimeConvergence(),JointCheckStatus::Passed);
    ASSERT_EQ(fit.TargetRuntimeConvergence(),JointCheckStatus::Passed);
    ASSERT_TRUE(fit.assembled_state); ASSERT_TRUE(fit.objective);
    const auto & state=*fit.assembled_state;
    std::vector<double> prediction(problem.Input().observations.size(),0.);
    for (std::size_t atom=0;atom<state.tau.size();++atom)
    {
        EXPECT_GT(state.tau[atom],0.); EXPECT_LE(state.tau[atom],bound*(1.+1e-14));
        EXPECT_GE(state.ac[2*atom],0.);
        for (const auto & support:problem.Input().support[atom])
        {
            const auto kernel=n::EvaluateKernel(support.squared_distance,state.tau[atom],2.5);
            prediction[support.row]+=state.ac[2*atom]*kernel.gaussian+state.ac[2*atom+1]*kernel.charge;
        }
    }
    double objective=0.;
    for (std::size_t row=0;row<prediction.size();++row)
    {
        const double residual=(prediction[row]-problem.Input().observations[row])/problem.ObservationScale();
        objective+=residual*residual/2.;
    }
    EXPECT_NEAR(objective,*fit.objective,1e-18);
    for (const auto & component:fit.components)
    {
        EXPECT_LE(component.search_sweeps,3u);
        EXPECT_NE(component.stop_reason,"block-sweep-budget");
        for (const auto * name:{"conditional-ac-optimality","width-stationarity",
            "independent-replay","endpoint-verification","width-bound-feasibility"})
        {
            const auto check=std::find_if(component.evidence.begin(),component.evidence.end(),
                [&](const auto & c){return c.name==name;});
            ASSERT_NE(check,component.evidence.end())<<name;
            EXPECT_EQ(check->status,JointCheckStatus::Passed)<<name<<": "<<check->reason;
        }
    }
}
}

TEST(JointSeedSensitivityTest, WidthTruthRecoveryAndQualificationAtTwoVoxelSpacings)
{
    for (double spacing:{.2,.4})
    {
        auto source=second_stage_test::FirstStageStudyModel(5.5,0);
        auto map=second_stage_test::FirstStageStudyMap(*source,spacing,1e-4,0);
        const auto problem=BuildPureHaloSupportExclusionProblem(*map,*source);
        for (const auto & policy:policies)
        {
            ModelObject model(*source); model.EditAnalysis().InitializeFromSelection();
            FitOptions options; options.quiet_mode=true;
            const auto workset=d::MakeJointFittingWorkset(model,problem);
            const auto init=d::PrepareJointWidthInitialization(*map,model,workset,options,policy);
            const auto calibration=n::CalibrateGlobalWidth(init,"synthetic-truth-regression");
            EXPECT_EQ(calibration.status,"Available"); EXPECT_EQ(calibration.qualified_donor_count,24u);
            EXPECT_EQ(calibration.qualified_donor_coverage,1.); EXPECT_EQ(calibration.fallback_count,0u);
            double square=0.;
            for (std::size_t i=0;i<init.atoms.size();++i)
            {
                const double error=init.atoms[i].mdpde[1]-(.35+.005*static_cast<double>(i));
                square+=error*error;
            }
            // Truth-based bounds cover interpolation error, not just solver self-consistency.
            const bool cubic=policy.interpolation==d::MapInterpolationMethod::Tricubic;
            EXPECT_LT(std::sqrt(square/24.),spacing==.2 ? (cubic ? .002 : .012) : (cubic ? .02 : .045));
        }
    }
}

TEST(JointSeedSensitivityTest, ProductionAndCommonBoundsRetainCertifiedEquivalentEndpoints)
{
    auto source=second_stage_test::FirstStageStudyModel(5.5,0);
    auto map=second_stage_test::FirstStageStudyMap(*source,.4,1e-4,0);
    const auto problem=BuildPureHaloSupportExclusionProblem(*map,*source);
    std::vector<JointInitialization> initializations;
    std::vector<n::WidthBoundCalibration> calibrations;
    double common_bound=0.;
    for (const auto & policy:policies)
    {
        ModelObject model(*source); model.EditAnalysis().InitializeFromSelection();
        FitOptions options; options.quiet_mode=true;
        initializations.push_back(d::PrepareJointWidthInitialization(*map,model,
            d::MakeJointFittingWorkset(model,problem),options,policy));
        calibrations.push_back(n::CalibrateGlobalWidth(initializations.back(),"production-calibration"));
        ASSERT_TRUE(calibrations.back().tau_max4);
        common_bound=std::max(common_bound,*calibrations.back().tau_max4);
    }
    std::optional<JointState> reference;
    for (std::size_t i=0;i<policies.size();++i)
    {
        const auto production=FitJointComponentsInitialized(problem,initializations[i]);
        CheckEndpoint(problem,production,*calibrations[i].tau_max4);
        auto calibration=calibrations[i]; calibration.tau_max4=common_bound;
        const auto projected=n::MakeBoundedInitialization(initializations[i],calibration);
        n::BcdExecutionPolicy execution; execution.eta_upper_bound=std::log(common_bound);
        const auto controlled=n::FitBcdComponents(problem,projected.tau_values,execution);
        CheckEndpoint(problem,controlled,common_bound);
        ASSERT_TRUE(controlled.assembled_state); ASSERT_TRUE(production.assembled_state);
        if (!reference) reference=controlled.assembled_state;
        EXPECT_NEAR(*production.objective,*controlled.objective,1e-12);
        for (std::size_t atom=0;atom<24;++atom)
        {
            EXPECT_NEAR(controlled.assembled_state->tau[atom],reference->tau[atom],1e-6);
            EXPECT_NEAR(production.assembled_state->tau[atom],reference->tau[atom],1e-6);
        }
    }
}

TEST(JointSeedSensitivityTest, ExperimentalSamplesPersistWithoutBecomingGridIdentities)
{
    auto model=second_stage_test::FirstStageStudyModel(5.5,0);
    auto map=second_stage_test::FirstStageStudyMap(*model,.4,1e-4,0);
    FitOptions options; options.quiet_mode=true;
    d::JointSeedSamplingReport report;
    d::RunJointComponentWorkflow(*map,*model,options,nullptr,policies[3],&report);
    ASSERT_EQ(report.policy.sampling,policies[3].sampling);
    ASSERT_TRUE(model->GetAnalysisView().GetJointResult());
    const auto original=AtomLocalPotentialView::For(*model->FindAtomPtr(1));
    const auto raw=original.GetRawSamplingEntries(false);
    ASSERT_EQ(raw.size(),161u); ASSERT_TRUE(original.GetPostFitPeeling());
    ASSERT_TRUE(original.GetPostJointLocalResult());
    const auto peeling_count=original.GetPostFitPeeling()->samples.size();
    const auto retained=original.GetPostJointLocalResult()->effective_sample_count;
    command_test::ScopedTempDir directory("joint-seed-policy");
    DataRepository repository(directory.path()/"samples.sqlite");
    repository.SaveModel(*model,"Q1"); auto restored=repository.LoadModel("Q1");
    const auto view=AtomLocalPotentialView::For(*restored->FindAtomPtr(1));
    const auto saved=view.GetRawSamplingEntries(false);
    ASSERT_EQ(saved.size(),raw.size());
    for (std::size_t i=0;i<raw.size();++i)
    {
        EXPECT_EQ(saved[i].point.position,raw[i].point.position);
        EXPECT_DOUBLE_EQ(saved[i].point.distance,raw[i].point.distance);
        EXPECT_DOUBLE_EQ(saved[i].response,raw[i].response);
        EXPECT_EQ(saved[i].point.is_selected,raw[i].point.is_selected);
    }
    restored->EditAnalysis().SetAtomLocalRawSamplingEntries(*restored->FindAtomPtr(1),
        SampleAtomMapValues(*map,*restored->FindAtomPtr(1)));
    ASSERT_TRUE(view.GetPostFitPeeling()); ASSERT_TRUE(view.GetPostJointLocalResult());
    EXPECT_EQ(view.GetPostFitPeeling()->samples.size(),peeling_count);
    EXPECT_EQ(view.GetPostJointLocalResult()->effective_sample_count,retained);
    EXPECT_EQ(view.GetRawSamplingEntries(false).size(),200u);
}


TEST(JointSeedSensitivityTest, FixedGeometryRelaxationIncreasesCrowdedWidthErrorAndContamination)
{
    auto source=second_stage_test::FirstStageStudyModel(1.4,0);
    auto map=second_stage_test::FirstStageStudyMap(*source,.4,1e-4,0);
    FitOptions options; options.quiet_mode=true;
    std::vector<double> rmse,contamination;
    for (int control=0;control<3;++control)
    {
        ModelObject model(*source); model.EditAnalysis().InitializeFromSelection();
        const auto workset=d::MakeTwoStageFittingWorkset(model);
        d::RunJointSeedSamplingWorkflow(*map,model,workset.contributors,1,{},true,{});
        double neighbor_sum=0.; std::size_t retained=0;
        for (const auto * atom:workset.contributors)
        {
            auto raw=AtomLocalPotentialView::For(*atom).GetRawSamplingEntries(false);
            if (control)
                second_stage_test::ApplyFibonacciNeighborSelection(raw,
                    second_stage_test::FibonacciNeighborMasks(*atom,model,raw),control==2);
            // Truth is used only to assess selected observations, never to fit them.
            for (const auto & sample:raw)
                if (sample.point.is_selected && sample.point.distance<=1.)
                {
                    ++retained;
                    for (std::size_t i=0;i<model.GetAtomList().size();++i)
                    {
                        const auto & other=*model.GetAtomList()[i];
                        if (other.GetSerialID()==atom->GetSerialID()) continue;
                        double square=0.;
                        for (std::size_t axis=0;axis<3;++axis)
                        {
                            const double delta=sample.point.position[axis]-other.GetPositionRef()[axis];
                            square+=delta*delta;
                        }
                        neighbor_sum+=(1.5+.02*static_cast<double>(i))*
                            n::EvaluateKernel(square,.35+.005*static_cast<double>(i),2.5).gaussian;
                    }
                }
            model.EditAnalysis().SetAtomLocalRawSamplingEntries(*atom,std::move(raw));
        }
        const auto init=d::RunJointFirstStageInitializationFromPreparedSamples(model,workset,options);
        const auto calibration=n::CalibrateGlobalWidth(init,"neighbor-selection-regression");
        EXPECT_EQ(calibration.qualified_donor_count,24u); EXPECT_EQ(calibration.fallback_count,0u);
        double square=0.;
        for (std::size_t i=0;i<init.atoms.size();++i)
        {
            const double error=init.atoms[i].mdpde[1]-(.35+.005*static_cast<double>(i));
            square+=error*error;
        }
        rmse.push_back(std::sqrt(square/24.)); contamination.push_back(neighbor_sum/static_cast<double>(retained));
    }
    EXPECT_LT(rmse[0],.03); EXPECT_GT(rmse[1],.045); EXPECT_GT(rmse[2],.06);
    EXPECT_LT(rmse[0],rmse[1]); EXPECT_LT(rmse[1],rmse[2]);
    EXPECT_LT(contamination[0],contamination[1]); EXPECT_LT(contamination[1],contamination[2]);
}

TEST(JointSeedSensitivityTest, NeighborControlEndpointsUseOriginalGridDownstreamAndPersistMasks)
{
    for (int control=0;control<3;++control)
    {
        // Keep the Production donor minimum; crowd only the first pair.
        auto source=second_stage_test::FirstStageStudyModel(5.5,0);
        source->FindAtomPtr(2)->SetPosition(1.53,.17,.19);
        ModelObject model(*source); model.EditAnalysis().InitializeFromSelection();
        auto map=second_stage_test::FirstStageStudyMap(model,.4,1e-4,0);
        const auto problem=BuildPureHaloSupportExclusionProblem(*map,model);
        const auto workset=d::MakeJointFittingWorkset(model,problem);
        d::RunJointSeedSamplingWorkflow(*map,model,workset.contributors,1,{},true,{});
        if (control)
            for (const auto * atom:workset.contributors)
            {
                auto raw=AtomLocalPotentialView::For(*atom).GetRawSamplingEntries(false);
                second_stage_test::ApplyFibonacciNeighborSelection(raw,
                    second_stage_test::FibonacciNeighborMasks(*atom,model,raw),control==2);
                model.EditAnalysis().SetAtomLocalRawSamplingEntries(*atom,std::move(raw));
            }
        FitOptions options; options.quiet_mode=true;
        const auto init=d::RunJointFirstStageInitializationFromPreparedSamples(model,workset,options);
        const auto fit=FitJointComponentsInitialized(problem,init);
        ASSERT_EQ(fit.RuntimeConvergence(),JointCheckStatus::Passed);
        ASSERT_EQ(fit.TargetRuntimeConvergence(),JointCheckStatus::Passed);
        const auto raw=AtomLocalPotentialView::For(*workset.contributors[0]).GetRawSamplingEntries(false);
        d::JointPostprocessingMetrics metrics;
        d::RunJointPostprocessing(*map,model,problem,CaptureJointAnalysisResult(fit),options,metrics);
        EXPECT_EQ(metrics.available_peeling_sample_count,metrics.grid_peeling_sample_count);
        const auto view=AtomLocalPotentialView::For(*workset.contributors[0]);
        ASSERT_TRUE(view.GetPostFitPeeling()); ASSERT_TRUE(view.GetPostJointLocalResult());
        EXPECT_EQ(view.GetPostJointLocalResult()->status,EvidenceStatus::Available);
        command_test::ScopedTempDir directory("neighbor-selection-persistence");
        DataRepository repository(directory.path()/"samples.sqlite");
        repository.SaveModel(model,"control"); auto restored=repository.LoadModel("control");
        const auto saved=AtomLocalPotentialView::For(*restored->FindAtomPtr(1)).GetRawSamplingEntries(false);
        ASSERT_EQ(saved.size(),raw.size());
        for (std::size_t k=0;k<raw.size();++k)
        {
            EXPECT_EQ(saved[k].point.position,raw[k].point.position);
            EXPECT_EQ(saved[k].point.distance,raw[k].point.distance);
            EXPECT_EQ(saved[k].response,raw[k].response);
            EXPECT_EQ(saved[k].point.is_selected,raw[k].point.is_selected);
        }
    }
}
