#include <gtest/gtest.h>
#include "support/JointComponentChecks.hpp"
#include "support/JointLocalEvidence.hpp"
#include "support/CommandTestHelpers.hpp"
#include <fstream>

namespace {
namespace p=second_stage_test::matched::joint_abc;
using Vector=Eigen::VectorXd;
}
namespace {
struct TwoBlocks
{
    p::Domain domain{81,{{},{}}};
    Vector y=Vector::Zero(81),eta=Vector::Constant(2,std::log(.55)),beta=Vector::Zero(4);
    TwoBlocks(double amplitude=2)
    {
        beta<<2,.2,amplitude,-.15;
        for(int block=0;block<2;++block) for(int k=0;k<40;++k)
        {
            const auto row=block*40+k; const double square=.003*k*k;
            domain.atoms[static_cast<std::size_t>(block)].push_back({row,square});
            const auto b=second_stage_test::matched::EvaluateBasis(square,.5,2.5);
            y(row)=beta(2*block)*b.gaussian+beta(2*block+1)*b.charge+.001*std::sin(k);
        }
        y(80)=3;
    }
};
}
namespace {
boost::json::value Science(boost::json::value value)
{
    if(value.is_array()) for(auto & entry:value.as_array()) entry=Science(entry);
    if(value.is_object())
    {
        std::vector<std::string> erase;
        for(auto & entry:value.as_object())
            if(entry.key().ends_with("seconds") || entry.key()=="process_peak_rss_bytes") erase.emplace_back(entry.key());
            else entry.value()=Science(entry.value());
        for(const auto & key:erase) value.as_object().erase(key);
    }
    return value;
}
}
TEST(JointComponentChecksTest, StructuralGraphPreservesBridgeChainConstantAndUnobserved)
{
    const p::Domain domain(6,{{{0,0},{1,6.25}},{{1,1},{2,1}},{{2,1}},{{3,1}}, {}});
    const auto c=p::MakeContext(Vector::Ones(6),5); const auto part=p::BuildPartition(domain,{"a","b","c","d","e"});
    ASSERT_EQ(part.components.size(),3); EXPECT_EQ(part.components[0].atoms,(std::vector<Eigen::Index>{0,1,2}));
    EXPECT_EQ(part.constant_rows,(std::vector<Eigen::Index>{4,5})); EXPECT_EQ(part.unobserved_atoms,(std::vector<Eigen::Index>{4}));
    EXPECT_EQ(part.components[0].domain.atoms[0][1].square,6.25);
    EXPECT_EQ(p::Census(domain,part,c).at("constant_objective"),1.0);
}

TEST(JointComponentChecksTest, SharedParentScaleAndDirectionsAreNotRenormalized)
{
    const p::Domain domain(4,{{{0,0},{1,1}},{{2,0},{3,1}}});
    auto c=p::MakeContext((Vector(4)<<100,200,.01,.02).finished(),2);
    c.audit.directions=Eigen::MatrixXd::Constant(2,1,1/std::sqrt(2.0));
    const auto part=p::BuildPartition(domain,c.atom_ids);
    const auto child=p::ComponentContext(c,part.components[1],true);
    EXPECT_EQ(child.scale,c.scale); EXPECT_EQ(child.observations,c.observations);
    EXPECT_EQ(child.audit.directions(0,0),c.audit.directions(1,0)); EXPECT_LT(child.audit.directions.col(0).norm(),1);
    EXPECT_EQ(child.rank.rows,2); EXPECT_EQ(child.rank.design_columns,2);
    const auto shared=p::ComponentContext(c,part.components[1],false); EXPECT_EQ(shared.rank.rows,4);
}

TEST(JointComponentChecksTest, PartitionIDsSurviveStoragePermutationAndRejectMalformedCSR)
{
    const p::Domain a(3,{{{0,0},{1,1}},{{1,2}},{{2,3}}}),b(3,{{{0,3}},{{1,2}},{{2,0},{1,1}}});
    const auto x=p::BuildPartition(a,{"A","B","C"}),y=p::BuildPartition(b,{"C","B","A"});
    ASSERT_EQ(x.components.size(),y.components.size());
    for(std::size_t k=0;k<x.components.size();++k) EXPECT_EQ(x.components[k].id,y.components[k].id);
    EXPECT_THROW(p::BuildPartition(p::Domain(1,{{{0,1},{0,1}}}),{"A"}),std::invalid_argument);
    EXPECT_THROW(p::BuildPartition(p::Domain(1,{{{1,1}}}),{"A"}),std::invalid_argument);
}

TEST(JointComponentChecksTest, RegisteredAuditDoesNotShrinkWithAComponent)
{
    const auto plan=p::RegisteredAudit(168,"heterogeneous-168","first-stage-double");
    const p::Domain domain(2,{{{0,0}},{{1,0}}});
    const auto c=p::MakeContext(Vector::Ones(2),2,"snapshot",&plan);
    const auto part=p::BuildPartition(domain,c.atom_ids);
    const auto child=p::ComponentContext(c,part.components[0],true);
    EXPECT_FALSE(child.audit.trial_details); EXPECT_FALSE(child.audit.expanded_if_unverified);
    EXPECT_TRUE(p::RegisteredAudit(12,"near-0.02","narrower-double").precision);
}

TEST(JointComponentChecksTest, EvaluationUsesParentScaleForBothKKTAndWidthGradient)
{
    std::vector<p::Support> support;
    for(int k=0;k<40;++k) support.push_back({k,.003*k*k});
    const p::Domain domain(40,{support}); const Vector eta=Vector::Constant(1,std::log(.5));
    Vector y(40);
    for(int k=0;k<40;++k)
    {
        const auto basis=second_stage_test::matched::EvaluateBasis(.003*k*k,.5,2.5);
        y(k)=2*basis.gaussian-.3*basis.charge+.01*std::sin(k);
    }
    Vector parent(41); parent.head(40)=y; parent(40)=1000;
    auto context=p::MakeContext(parent,1); context.independent_search=true;
    const auto local=p::Evaluate(domain,y,eta),global=p::Evaluate(domain,y,eta,false,&context);
    ASSERT_TRUE(local.valid && global.valid);
    EXPECT_EQ(local.beta,global.beta);
    const double local_scale=std::max(1.0,y.norm());
    EXPECT_NEAR(global.gradient(0),local.gradient(0)*std::pow(local_scale/context.scale,2),1e-18);
    EXPECT_NEAR(boost::json::value_to<double>(global.certificate.at("projected_kkt")),
        boost::json::value_to<double>(local.certificate.at("projected_kkt"))*local_scale/context.scale,1e-16);
    EXPECT_TRUE(p::Trust(domain,y,global,&context).at("passed").as_bool());
}

TEST(JointComponentChecksTest, ZeroObservationsAndUnderflowDoNotChangeStructuralPartition)
{
    const p::Domain domain(3,{{{0,0},{1,6.25}},{{1,6.25},{2,0}}});
    const auto c=p::MakeContext(Vector::Zero(3),2);
    const auto before=p::BuildPartition(domain,c.atom_ids);
    (void)p::Evaluate(domain,Vector::Zero(3),Vector::Constant(2,-100));
    const auto after=p::BuildPartition(domain,c.atom_ids);
    ASSERT_EQ(before.components.size(),1); EXPECT_EQ(before.atom_component,after.atom_component);
}

TEST(JointComponentChecksTest, SameStateIncludesResidualTermAndConstantRows)
{
    const TwoBlocks f; const auto c=p::MakeContext(f.y,2); const auto part=p::BuildPartition(f.domain,c.atom_ids);
    const auto parity=p::SameState(f.domain,f.y,f.eta,f.beta,part,c);
    EXPECT_TRUE(parity.at("passed").as_bool())<<parity;
    EXPECT_TRUE(parity.at("full_equivalence").as_bool());
    EXPECT_EQ(parity.at("raw").at("constant_objective"),4.5);
    RecordProperty("scientific_record",boost::json::serialize(parity));
    const auto e=p::EvaluatePartitioned(f.domain,f.y,f.eta,part,c);
    const auto d=p::DifferentiatePartitioned(e,part,c);
    ASSERT_TRUE(e.valid && d.valid); EXPECT_GT((d.projected-d.jacobian).norm(),1e-5);
    Vector direction(2); direction<<.3,-.7; constexpr double h=1e-5;
    const auto plus=p::EvaluatePartitioned(f.domain,f.y,f.eta+h*direction,part,c,true);
    const auto minus=p::EvaluatePartitioned(f.domain,f.y,f.eta-h*direction,part,c,true);
    ASSERT_TRUE(plus.valid && minus.valid);
    EXPECT_LT(((plus.residual-minus.residual)/(2*h*c.scale)-d.jacobian*direction).norm(),1e-8);
}

TEST(JointComponentChecksTest, FiniteUntrustedStateKeepsRawParityWithoutProfileClaims)
{
    TwoBlocks f; f.eta.setConstant(std::log(1e9));
    const auto c=p::MakeContext(f.y,2); const auto part=p::BuildPartition(f.domain,c.atom_ids);
    const auto parity=p::SameState(f.domain,f.y,f.eta,f.beta,part,c);
    RecordProperty("scientific_record",boost::json::serialize(parity));
    EXPECT_TRUE(parity.at("raw").at("passed").as_bool());
    EXPECT_TRUE(parity.at("passed").as_bool()); EXPECT_FALSE(parity.at("full_equivalence").as_bool());
    EXPECT_FALSE(parity.at("profile_solves").at(0).at("available").as_bool());
}

TEST(JointComponentChecksTest, BlockBackendSharesConstrainedActiveSet)
{
    const TwoBlocks f(-2); const auto c=p::MakeContext(f.y,2); const auto part=p::BuildPartition(f.domain,c.atom_ids);
    for(bool reference:{false,true})
    {
        const auto a=p::Evaluate(f.domain,f.y,f.eta,reference,&c),b=p::EvaluatePartitioned(f.domain,f.y,f.eta,part,c,reference);
        ASSERT_TRUE(a.valid && b.valid)<<a.reason<<' '<<b.reason;
        EXPECT_EQ(a.certificate.at("active_atoms"),b.certificate.at("active_atoms"));
        EXPECT_FALSE(a.certificate.at("active_atoms").as_array().empty()); EXPECT_LT((a.beta-b.beta).norm(),1e-12);
        EXPECT_EQ(a.certificate.at("linear_solves"),b.certificate.at("linear_solves"));
    }
}

TEST(JointComponentChecksTest, RankEvidenceUsesGlobalDimensionsAndCurrentSpectrum)
{
    const auto c=p::MakeContext(Vector::Ones(1000),2);
    Eigen::MatrixXd small(2,2); small<<1,0,0,1e-14;
    const auto global=p::MatrixSpectrum(small,c.rank,4,false);
    const p::RankPolicy local{2,2,1}; const auto block=p::MatrixSpectrum(small,local,2,false);
    EXPECT_EQ(global.at("rank"),1); EXPECT_EQ(block.at("rank"),2);
    const auto changed=p::MatrixSpectrum(Eigen::MatrixXd(2*small),c.rank,4,false);
    EXPECT_DOUBLE_EQ(boost::json::value_to<double>(changed.at("rank_threshold")),2*boost::json::value_to<double>(global.at("rank_threshold")));
}

TEST(JointComponentChecksTest, FailureIsolationDoesNotFillMissingRowsWithZero)
{
    TwoBlocks f; f.domain.atoms.push_back(f.domain.atoms[1]);
    auto c=p::MakeContext(f.y,3); const auto part=p::BuildPartition(f.domain,c.atom_ids);
    const Vector initial=Vector::Constant(3,.55);
    const auto result=p::FitComponents(f.domain,f.y,initial,part,c);
    ASSERT_EQ(result.at("components").as_array().size(),2);
    EXPECT_TRUE(result.at("components").at(0).at("search_success").as_bool());
    EXPECT_FALSE(result.at("components").at(1).at("usable_state").as_bool());
    EXPECT_FALSE(result.at("prediction_available").as_bool()); EXPECT_TRUE(result.at("prediction").is_null());
    EXPECT_FALSE(result.at("joint_qualified").as_bool());
    for(std::size_t r=0;r<81;++r) EXPECT_EQ(result.at("available_row_mask").at(r).as_bool(),r<40 || r==80);
    const auto failed_first=p::FitComponent(part.components[1],f.y,initial,c);
    const auto healthy_second=p::FitComponent(part.components[0],f.y,initial,c);
    EXPECT_EQ(Science(failed_first),Science(result.at("components").at(1)));
    EXPECT_EQ(Science(healthy_second),Science(result.at("components").at(0)));
}

TEST(JointComponentChecksTest, BudgetAndInvalidStartRetainHonestAvailability)
{
    const TwoBlocks f; auto c=p::MakeContext(f.y,2); const auto part=p::BuildPartition(f.domain,c.atom_ids);
    const Vector initial=Vector::Constant(2,.6);
    const auto healthy=p::FitComponent(part.components[0],f.y,initial,c);
    auto limited=c; limited.profile_budget=1;
    const auto exhausted=p::FitComponent(part.components[1],f.y,initial,limited);
    EXPECT_FALSE(exhausted.at("search_success").as_bool()); EXPECT_TRUE(exhausted.at("usable_state").as_bool());
    EXPECT_EQ(exhausted.at("stop_reason"),"profile-budget");
    const auto assembled=p::Assemble(f.domain,f.y,part,c,boost::json::array{healthy,exhausted});
    EXPECT_TRUE(assembled.at("prediction_available").as_bool()); EXPECT_FALSE(assembled.at("joint_qualified").as_bool());
    EXPECT_TRUE(assembled.at("search_stopped_without_convergence").as_bool());
    Vector invalid=initial; invalid(1)=0;
    const auto failure=p::FitComponent(part.components[1],f.y,invalid,c);
    EXPECT_FALSE(failure.at("usable_state").as_bool());
    EXPECT_EQ(Science(healthy),Science(p::FitComponent(part.components[0],f.y,initial,c)));
}

TEST(JointComponentChecksTest, UnobservedAtomsKeepRawPredictionButDisableProfileClaims)
{
    TwoBlocks f; f.domain.atoms[1].clear(); const auto c=p::MakeContext(f.y,2);
    const auto part=p::BuildPartition(f.domain,c.atom_ids);
    const auto state=p::SameState(f.domain,f.y,f.eta,f.beta,part,c);
    EXPECT_TRUE(state.at("raw").at("passed").as_bool()); EXPECT_TRUE(state.at("passed").as_bool());
    EXPECT_FALSE(state.at("full_equivalence").as_bool());
    EXPECT_FALSE(state.at("profile_solves").at(0).at("available").as_bool());
    EXPECT_EQ(part.unobserved_atoms,(std::vector<Eigen::Index>{1}));
    const auto fit=p::FitComponents(f.domain,f.y,Vector(f.eta.array().exp()),part,c);
    EXPECT_TRUE(fit.at("components").at(0).at("search_success").as_bool());
    EXPECT_FALSE(fit.at("components").at(1).at("usable_state").as_bool());
    EXPECT_FALSE(fit.at("prediction_available").as_bool());
    RecordProperty("fit_availability",boost::json::serialize(fit));
}

TEST(JointComponentChecksTest, LocalAuditUsesActualStateAndDiscardsSiblingDirections)
{
    const TwoBlocks f;
    auto parent=p::MakeContext(f.y,2);
    const auto part=p::BuildPartition(f.domain,parent.atom_ids);
    const auto & view=part.components[0];
    auto fit=p::FitComponent(view,f.y,Vector::Constant(2,.55),parent);
    ASSERT_TRUE(fit.at("usable_state").as_bool());
    const auto y=p::Select(f.y,view.rows);
    auto child=p::ComponentContext(parent,view,true);
    child.audit.directions=Eigen::MatrixXd::Zero(1,3);
    const auto before=second_stage_test::matched::certification::PrepareLocalAudit(view.domain,y,fit,child);
    (void)p::FitComponent(part.components[1],f.y,Vector::Constant(2,.55),parent);
    child.audit.directions=Eigen::MatrixXd::Constant(1,3,1e-100);
    const auto after=second_stage_test::matched::certification::PrepareLocalAudit(view.domain,y,fit,child);
    const auto rerun=p::FitComponent(view,f.y,Vector::Constant(2,.55),parent);
    const auto reverse=second_stage_test::matched::certification::PrepareLocalAudit(view.domain,y,rerun,child);
    EXPECT_EQ(Science(fit),Science(rerun));
    EXPECT_EQ(Science(before.fit),Science(reverse.fit));
    EXPECT_EQ(before.scope,reverse.scope);
    EXPECT_EQ(before.fit,after.fit);
    EXPECT_EQ(before.scope,after.scope);
    EXPECT_EQ(before.context.audit.directions,after.context.audit.directions);
    EXPECT_EQ(before.context.scale,parent.scale);
    EXPECT_EQ(before.context.rank.rows,view.domain.rows);
    EXPECT_EQ(before.fit.at("primary").at("beta"),fit.at("last_trusted_state").at("beta"));
    ASSERT_TRUE(before.scope.at("weak_direction_available").as_bool());
    EXPECT_DOUBLE_EQ(before.context.audit.directions.col(0).norm(),1);
    fit["usable_state"]=false;
    const auto missing=second_stage_test::matched::certification::PrepareLocalAudit(view.domain,y,fit,child);
    EXPECT_FALSE(missing.fit.at("primary").at("valid").as_bool());
    EXPECT_FALSE(missing.scope.at("weak_direction_available").as_bool());
}
