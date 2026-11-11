#include <gtest/gtest.h>
#include "support/JointQScoreFixture.hpp"
#include "core/painter/detail/PotentialPlotBuilder.hpp"
#include "data/detail/LocalPotentialEntry.hpp"
#include "data/detail/ModelAnalysisData.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

#ifdef HAVE_ROOT
#include <TGraphErrors.h>
#endif

namespace {
using namespace qscore_test;
// Scalar oracle uses the mathematical kernel, separate from GaussianModel3D.
double Response(const GaussianModel3D & m,double r)
{
    const double w=m.GetWidth();
    return m.GetAmplitude()*std::exp(-r*r/(2*w*w))/std::pow(2*std::numbers::pi*w*w,1.5)+
        m.GetOffset()*(r==0 ? std::sqrt(2/std::numbers::pi)/w : std::erf(r/(std::sqrt(2.)*w))/r);
}
void ExpectOracle(Fixture & f)
{
    const auto & result=f.results.front().second;
    ASSERT_TRUE(result.scores) << result.reason;
    const auto & target=*f.model->FindAtomPtr(1);
    const auto samples=PrepareQScoreSamples(target,*f.map,*f.model);
    std::vector<double> y,standard,joint;
    std::size_t evaluations{};
    const auto & target_point=AtomLocalPotentialView::For(target).GetJointStageModel();
    for (const auto & sample:samples)
    {
        long double response=sample.response;
        for (const auto & atom:f.model->GetAtomList())
        {
            if (atom.get()==&target || atom->GetElement()==Element::HYDROGEN) continue;
            double square{};
            for (std::size_t axis=0;axis<3;++axis)
            {const double delta=sample.point.position[axis]-atom->GetPositionRef()[axis];square+=delta*delta;}
            if (square>6.25) continue;
            GaussianModel3D point{0,1,0};
            const auto view=AtomLocalPotentialView::For(*atom);
            if (view.IsAvailable() && view.HasJointStageModel()) point=view.GetJointStageModel();
            else
            {
                const GeometryLocalTargetDonorIndex index(f.halos.target_donors);
                const auto median=GeometryLocalTargetMedian(index.Query(atom->GetPositionRef()));
                point=GaussianModel3D{median[0],median[1],median[2]};
            }
            response-=Response(point,std::sqrt(square));
            ++evaluations;
        }
        y.push_back(static_cast<double>(response));
        const double r=sample.point.distance;
        standard.push_back(f.model->GetReferenceHeight()*std::exp(-r*r/(2*.6*.6))+f.model->GetReferenceOffset());
        joint.push_back(Response(target_point,r));
    }
    const auto green=CalculateQScoreCorrelation(y,standard),red=CalculateQScoreCorrelation(y,joint);
    ASSERT_TRUE(green.value); ASSERT_TRUE(red.value);
    EXPECT_NEAR(result.scores->standard,*green.value,2e-14);
    EXPECT_NEAR(result.scores->joint,*red.value,2e-14);
    EXPECT_EQ(result.sample_count,samples.size());
    EXPECT_EQ(result.evaluated_contributions,evaluations);
    for (const auto score:{result.scores->standard,result.scores->joint})
    { EXPECT_TRUE(std::isfinite(score));EXPECT_GE(score,-1.-1e-12);EXPECT_LE(score,1.+1e-12); }
}
}

TEST(PostJointQScoreTest, IsolatedAtomPreservesBlueAndEquivalentReferences)
{
    auto f=Make();
    ASSERT_TRUE(f.results.front().second.scores);
    const auto scores=*f.results.front().second.scores;
    EXPECT_EQ(scores.standard,f.model->FindAtomPtr(1)->GetStandardQScore());
    EXPECT_NEAR(scores.standard,scores.joint,1e-12);
    EXPECT_EQ(f.metrics.qscore_contributions,0u);
    EXPECT_EQ(f.metrics.qscore_available,1u);
    ExpectOracle(f);
}

TEST(PostJointQScoreTest, OverlappingGaussianAndChargeUseOnePeeledDataset)
{
    for (const double charge:{0.,.35})
    {
        auto f=Make({{{0,0,0},{2,.6,charge}},{{1.2,0,0},{2.3,.55,-charge}}});
        ExpectOracle(f);
        EXPECT_EQ(f.metrics.qscore_available,2u);
        EXPECT_GT(f.metrics.qscore_contributions,0u);
        EXPECT_NE(f.results.front().second.scores->standard,f.model->FindAtomPtr(1)->GetStandardQScore());
    }
}

TEST(PostJointQScoreTest, ChargeReferenceIsKernelRatherThanConstantOffset)
{
    auto f=Make({{{0,0,0},{2,.6,.8}}});
    ExpectOracle(f);
    const auto original=f.results.front().second;
    ASSERT_TRUE(original.scores);
    EXPECT_GT(std::abs(original.scores->standard-original.scores->joint),1e-4);
    // Change only target reference. With no neighbors the peeled data and count stay fixed.
    auto changed=f.snapshot;
    changed.components.front().state->ac[1]=0;
    f.model->EditAnalysis().ApplyJointResult(changed,"qscore-test-run");
    const auto comparisons=BuildPostJointQScores(*f.map,*f.model,*f.problem,changed,&f.halos,f.metrics);
    ASSERT_TRUE(comparisons.front().second.scores);
    EXPECT_EQ(comparisons.front().second.sample_count,original.sample_count);
    EXPECT_EQ(comparisons.front().second.scores->standard,original.scores->standard);
    EXPECT_NE(comparisons.front().second.scores->joint,original.scores->joint);
}

TEST(PostJointQScoreTest, InclusiveSampleSupportAndDistantNeighborsAndHydrogen)
{
    for (const double z:{2.5,3.5})
    {
        auto f=Make({{{0,0,0},{2,.6,.1}},{{0,0,z},{2,.6,.1},false},
            {{0,0,-3.2},{2,.6,.1},false},{{.8,0,0},{2,.6,10},false,true}});
        ExpectOracle(f);
        EXPECT_EQ(f.results.front().second.contributor_count,2u);
        EXPECT_EQ(f.results.front().second.halo_count,2u);
        EXPECT_TRUE(f.results.front().second.extra_halos.empty());
    }
}

TEST(PostJointQScoreTest, ExtraDomainHaloUsesQualifiedLocalDonorsWithoutChangingSnapshot)
{
    auto f=Make({{{0,0,0},{2,.6,.2}},{{3.5,0,0},{2,.6,.2},false}},true);
    ASSERT_EQ(f.snapshot.atom_ids.size(),1u);
    ASSERT_TRUE(f.results.front().second.scores);
    const auto & extra=f.results.front().second.extra_halos;
    ASSERT_EQ(extra.size(),1u); EXPECT_EQ(extra.front().atom_id,2);
    EXPECT_EQ(extra.front().donor_atom_ids,std::vector<std::string>{"1"});
    EXPECT_FALSE(f.snapshot.halo_median_reconstruction);
    EXPECT_FALSE(AtomLocalPotentialView::For(*f.model->FindAtomPtr(2)).IsAvailable());
    ExpectOracle(f);
    auto stale=f.halos;
    stale.target_donors.front().point[2]+=1;
    const auto missing=BuildPostJointQScores(*f.map,*f.model,*f.problem,f.snapshot,&stale,f.metrics);
    EXPECT_FALSE(missing.front().second.scores);
    EXPECT_EQ(missing.front().second.reason,"contributor-unavailable");
}

TEST(PostJointQScoreTest, MissingContributorAndNonfinitePredictionAreUnavailable)
{
    auto f=Make({{{0,0,0},{2,.6,.1}},{{1.2,0,0},{2,.6,.1}}});
    auto & entry=*ModelAnalysisData::Of(*f.model).FindAtomLocalEntry(*f.model->FindAtomPtr(2));
    const auto original=entry.JointStageEstimate();
    auto missing=original;missing.point.reset();entry.SetJointStageEstimate(missing);
    f.Calculate();
    EXPECT_FALSE(f.results.front().second.scores);
    EXPECT_EQ(f.results.front().second.reason,"contributor-unavailable");
    EXPECT_EQ(f.results.front().second.missing_contributor_count,1u);
    missing=original;missing.source.run_id="other-run";entry.SetJointStageEstimate(missing);
    f.Calculate();EXPECT_FALSE(f.results.front().second.scores);
    entry.SetJointStageEstimate(original);
    auto nonfinite=f.snapshot;
    nonfinite.components.front().state->ac[2]=std::numeric_limits<double>::max();
    nonfinite.components.front().state->tau[1]=.01;
    f.model->EditAnalysis().ApplyJointResult(nonfinite,"qscore-test-run");
    const auto results=BuildPostJointQScores(*f.map,*f.model,*f.problem,nonfinite,&f.halos,f.metrics);
    EXPECT_FALSE(results.front().second.scores);
    EXPECT_EQ(results.front().second.reason,"nonfinite-contribution");
}

TEST(PostJointQScoreTest, BoundaryAndDegenerateDataHaveExplicitReasons)
{
    auto f=Make();
    f.model->SetReferenceHeight(0);
    f.Calculate();
    EXPECT_FALSE(f.results.front().second.scores);
    EXPECT_EQ(f.results.front().second.reason,"standard:reference-zero-variance");
    f.model->SetReferenceHeight(1);
    auto boundary=Make({Atom{}},false,true);
    EXPECT_EQ(boundary.results.front().second.reason,"sample-outside-map");
    EXPECT_FALSE(boundary.results.front().second.scores);
    auto coarse=Make({Atom{}},true);
    const auto size=coarse.map->GetMapValueArraySize();
    auto values=std::make_unique<double[]>(size);
    std::copy_n(coarse.map->GetMapValueArray(),size,values.get());
    // x=3 is outside the Original support, but enters QScore trilinear stencils.
    for (std::size_t grid=0;grid<size;++grid)
        if (coarse.map->GetGridPosition(grid)[0]==3) values[grid]=std::numeric_limits<double>::quiet_NaN();
    MapObject nonfinite(coarse.map->GetGridSize(),coarse.map->GetGridSpacing(),coarse.map->GetOrigin(),std::move(values));
    const auto invalid=BuildPostJointQScores(nonfinite,*coarse.model,*coarse.problem,coarse.snapshot,&coarse.halos,coarse.metrics);
    EXPECT_FALSE(invalid.front().second.scores);
    EXPECT_EQ(invalid.front().second.reason,"nonfinite-response");
    EXPECT_THROW(BuildPostJointQScores(nonfinite,*f.model,*f.problem,f.snapshot,&f.halos,f.metrics),std::invalid_argument);

}

TEST(QScorePainterTest, UsesPersistedSourcesSortedResiduesAndSharedCoverage)
{
    auto f=Make({{{0,0,0}},{{1.2,0,0}}});
    f.results.front().second.scores=PostJointQScoreValues{-.7,-.2};
    f.model->EditAnalysis().ApplyPostJointQScoresBatch(f.results);
    f.model->FindAtomPtr(1)->SetSequenceID(20);f.model->FindAtomPtr(2)->SetSequenceID(3);
    f.model->FindAtomPtr(1)->SetStandardQScore(.25);f.model->FindAtomPtr(2)->SetStandardQScore(.75);
    const auto blue=CollectResidueQScores(*f.model,QScoreMode::Standard);
    const auto green=CollectResidueQScores(*f.model,QScoreMode::JointPeeledStandard);
    const auto red=CollectResidueQScores(*f.model,QScoreMode::JointPeeledModel);
    EXPECT_EQ(blue.at("A").begin()->first,3);
    EXPECT_EQ(*blue.at("A").at(20).average,.25);
    EXPECT_EQ(green.at("A").at(20).average,f.results.front().second.scores->standard);
    EXPECT_EQ(red.at("A").at(20).average,f.results.front().second.scores->joint);
#ifdef HAVE_ROOT
    auto sorted=PotentialPlotBuilder(f.model.get()).CreateAverageQScoreToSequenceIDGraphMap(QScoreMode::Standard);
    ASSERT_EQ(sorted.at("A")->GetN(),2);
    EXPECT_EQ(sorted.at("A")->GetPointX(0),3);
    EXPECT_EQ(sorted.at("A")->GetPointX(1),20);
#endif
    f.model->FindAtomPtr(2)->SetSequenceID(20);
    for (const auto mode:{QScoreMode::JointPeeledStandard,QScoreMode::JointPeeledModel})
    {
        const auto paired=CollectResidueQScores(*f.model,mode);
        const auto & first=*f.results.front().second.scores;
        const auto & second=*f.results.back().second.scores;
        const double expected=mode==QScoreMode::JointPeeledStandard ?
            (first.standard+second.standard)/2 : (first.joint+second.joint)/2;
        EXPECT_EQ(*paired.at("A").at(20).average,expected);
        EXPECT_EQ(paired.at("A").at(20).available_atoms,2u);
    }
    const auto average=CollectResidueQScores(*f.model,QScoreMode::Standard);
    EXPECT_EQ(*average.at("A").at(20).average,.5);
    f.model->FindAtomPtr(2)->SetChainID("B");
    ModelAnalysisData::Of(*f.model).FindAtomLocalEntry(*f.model->FindAtomPtr(2))->ClearPostJointQScore();
    const auto missing=CollectResidueQScores(*f.model,QScoreMode::JointPeeledStandard);
    EXPECT_FALSE(missing.at("B").at(20).average);
    EXPECT_EQ(missing.at("B").at(20).available_atoms,0u);
#ifdef HAVE_ROOT
    PotentialPlotBuilder builder(f.model.get());
    auto graphs=builder.CreateAverageQScoreToSequenceIDGraphMap(QScoreMode::JointPeeledModel);
    EXPECT_FALSE(graphs.contains("B"));
    EXPECT_EQ(graphs.at("A")->GetPointY(0),f.results.front().second.scores->joint);
#endif
    const auto range=QScoreDisplayRange({-.7,.1,1});
    EXPECT_LE(range.first,-.7);EXPECT_EQ(range.second,1);
    EXPECT_EQ(QScoreDisplayRange({}),std::make_pair(-1.,1.));
    EXPECT_EQ(QScoreMethodLabels.size(),3u);
    EXPECT_NE(std::string(QScoreMethodLabels[2]).find("Joint model"),std::string::npos);
    f.model->SelectAtoms([](const auto &){return false;});
    EXPECT_TRUE(CollectResidueQScores(*f.model,QScoreMode::Standard).empty());
}

TEST(PostJointQScoreTest, SyntheticWorkloadRecordsCostAndPreservesEndpoints)
{
    std::vector<qscore_test::Atom> atoms;
    for (int z=-1;z<=1;++z) for (int y=-1;y<=1;++y) for (int x=-1;x<=1;++x)
        atoms.push_back({{2.*x,2.*y,2.*z},{2.,.55,.1}});
    auto f=qscore_test::Make(std::move(atoms));
    const auto before=f.snapshot.components;
    std::vector<double> blue;
    for (const auto & atom:f.model->GetAtomList()) blue.push_back(atom->GetStandardQScore());
    f.Calculate();
    EXPECT_EQ(f.metrics.qscore_targets,27u);
    EXPECT_EQ(f.metrics.qscore_available,27u);
    EXPECT_EQ(f.metrics.qscore_unavailable,0u);
    EXPECT_GT(f.metrics.qscore_contributions,f.metrics.qscore_samples);
    for (std::size_t c=0;c<before.size();++c)
    {
        EXPECT_EQ(f.snapshot.components[c].state->ac,before[c].state->ac);
        EXPECT_EQ(f.snapshot.components[c].state->tau,before[c].state->tau);
    }
    for (std::size_t i=0;i<blue.size();++i) EXPECT_EQ(f.model->GetAtomList()[i]->GetStandardQScore(),blue[i]);
    RecordProperty("qscore_targets",static_cast<int>(f.metrics.qscore_targets));
    RecordProperty("qscore_available",static_cast<int>(f.metrics.qscore_available));
    RecordProperty("qscore_unavailable",static_cast<int>(f.metrics.qscore_unavailable));
    RecordProperty("qscore_samples",static_cast<int>(f.metrics.qscore_samples));
    RecordProperty("qscore_contributions",static_cast<int>(f.metrics.qscore_contributions));
    RecordProperty("qscore_seconds",std::to_string(f.metrics.qscore_seconds));
    RecordProperty("qscore_estimated_payload_bytes",static_cast<int>(f.metrics.estimated_qscore_peak_bytes));
}
