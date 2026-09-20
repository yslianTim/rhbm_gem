#include <gtest/gtest.h>
#include "support/JointComponentChecks.hpp"
#include "support/JointPrecisionAudit.hpp"
#include "support/JointLocalEvidence.hpp"
#include "support/JointOfflineTools.hpp"
#include "support/CommandTestHelpers.hpp"
#include <fstream>
#include <boost/multiprecision/cpp_dec_float.hpp>

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
        auto support=domain.CopySupport();
        for(int block=0;block<2;++block) for(int k=0;k<40;++k)
        {
            const auto row=block*40+k; const double square=.003*k*k;
            support[static_cast<std::size_t>(block)].push_back({row,square});
            const auto b=second_stage_test::matched::EvaluateBasis(square,.5,2.5);
            y(row)=beta(2*block)*b.gaussian+beta(2*block+1)*b.charge+.001*std::sin(k);
        }
        domain=p::Domain(81,std::move(support));
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
TEST(JointComponentChecksTest, PrecisionReferenceUsesFullParentObservations)
{
    std::vector<p::Support> support; Vector y(30);
    for(int k=0;k<30;++k)
    {
        const double square=.004*k*k; support.push_back({k,square});
        const auto basis=second_stage_test::matched::EvaluateBasis(square,.5,2.5);
        y(k)=2*basis.gaussian+.2*basis.charge+.02*std::sin(k);
    }
    const p::Domain domain(30,{support}); const Vector eta=Vector::Constant(1,std::log(.6));
    Vector all(31); all.head(30)=y; all(30)=1000;
    const auto context=p::MakeContext(all,1);
    const auto local=p::Evaluate(domain,y,eta),parent=p::Evaluate(domain,y,eta,false,&context);
    ASSERT_TRUE(local.valid && parent.valid);
    const auto directions=Eigen::MatrixXd::Identity(1,1);
    const auto a=second_stage_test::matched::certification::PrecisionAudit(domain,y,local,directions);
    const auto b=second_stage_test::matched::certification::PrecisionAudit(domain,y,parent,directions,&context);
    ASSERT_TRUE(a.at("agreement_passed").as_bool() && b.at("agreement_passed").as_bool());
    EXPECT_TRUE(b.at("derivative_passed").as_bool());
    const double ratio=std::pow(std::max(1.0,y.norm())/context.scale,2);
    EXPECT_NEAR(boost::json::value_to<double>(b.at("b_gradient").at(0)),
        boost::json::value_to<double>(a.at("b_gradient").at(0))*ratio,1e-17);
    EXPECT_NEAR(boost::json::value_to<double>(b.at("local_correction_inf")),
        boost::json::value_to<double>(a.at("local_correction_inf")),1e-14);
}

TEST(JointComponentChecksTest, BlockMultiprecisionReferenceMatchesDenseGlobalReference)
{
    const TwoBlocks f; const auto c=p::MakeContext(f.y,2); const auto e=p::Evaluate(f.domain,f.y,f.eta,false,&c);
    ASSERT_TRUE(e.valid); Eigen::MatrixXd directions(2,2); directions<<.3,-.2,.7,.5;
    const auto dense=second_stage_test::matched::certification::PrecisionAudit(f.domain,f.y,e,directions,&c);
    const auto block=second_stage_test::matched::certification::PrecisionAudit(f.domain,f.y,e,directions,&c,true);
    ASSERT_TRUE(dense.at("agreement_passed").as_bool() && block.at("agreement_passed").as_bool());
    EXPECT_EQ(dense.at("fixed_face_feasible"),block.at("fixed_face_feasible"));
    for(const char * key:{"b_gradient","local_correction","derivative_relative_errors"})
        for(std::size_t k=0;k<dense.at(key).as_array().size();++k)
            EXPECT_NEAR(boost::json::value_to<double>(dense.at(key).at(k)),boost::json::value_to<double>(block.at(key).at(k)),1e-14);
}

TEST(JointComponentChecksTest, GlobalBlockBoundaryScanIncludesOtherComponentsAndConstants)
{
    TwoBlocks f; f.beta(2)=0; f.beta(3)=.2;
    auto c=p::MakeContext(f.y,2); c.audit.boundary_atoms={0,1}; c.audit.boundary=true;
    const auto dense=second_stage_test::matched::certification::BoundaryAudit(f.domain,f.y,f.eta,f.beta,&c);
    c.audit.block_precision=true; c.audit.cache_precision=true;
    second_stage_test::matched::certification::ResetPrecisionCache();
    const auto block=second_stage_test::matched::certification::BoundaryAudit(f.domain,f.y,f.eta,f.beta,&c);
    ASSERT_TRUE(dense.at("agreement_passed").as_bool() && block.at("agreement_passed").as_bool());
    using P=boost::multiprecision::cpp_dec_float_100;
    for(const char * precision:{"precision50","precision100"})
    {
        const auto & a=dense.at(precision).at("rows").as_array(); const auto & b=block.at(precision).at("rows").as_array();
        ASSERT_EQ(a.size(),b.size());
        for(std::size_t k=0;k<a.size();++k)
        {
            EXPECT_EQ(a[k].at("atom"),b[k].at("atom")); EXPECT_EQ(a[k].at("profile_valid"),b[k].at("profile_valid"));
            for(const char * key:{"prediction_change","structure_loss","path_objective","fixed_face_objective","profile_objective"})
            {
                const P x(boost::json::value_to<std::string>(a[k].at(key))),y(boost::json::value_to<std::string>(b[k].at(key)));
                EXPECT_LE(abs(x-y)/(1+abs(x)),P("1e-20"));
            }
        }
    }
    EXPECT_EQ(Science(block),Science(second_stage_test::matched::certification::BoundaryAudit(f.domain,f.y,f.eta,f.beta,&c)));
    EXPECT_GT(boost::json::value_to<int>(second_stage_test::matched::certification::PrecisionCacheCosts().at("boundary_reference_reuses")),0);
    second_stage_test::matched::certification::ResetPrecisionCache();
    EXPECT_EQ(second_stage_test::matched::certification::PrecisionCacheCosts().at("boundary_reference_reuses"),0);
}

TEST(JointComponentChecksTest, StandaloneBundleNeedsNoDatasetOrSiblingFiles)
{
    namespace j=boost::json;
    const TwoBlocks f; const auto parent=p::MakeContext(f.y,2);
    const auto partition=p::BuildPartition(f.domain,parent.atom_ids); const auto & view=partition.components[0];
    j::array observations,atoms,rows,memberships;
    for(double value:f.y) observations.push_back(value);
    for(auto atom:view.atoms) atoms.push_back(atom);
    for(auto row:view.rows) rows.push_back(row);
    for(std::size_t atom=0;atom<view.domain.atoms.size();++atom)
        for(const auto & support:view.domain.atoms[atom]) memberships.push_back(j::array{support.row,atom,support.square});
    j::object bundle{{"search_context",p::ContextEvidence(parent)},{"parent_observations",observations},
        {"component_initial_b",j::array{.55}},{"dataset","unit"},{"case","unit"},
        {"component_input",j::object{{"id",view.id},{"atoms",atoms},{"rows",rows},{"memberships",memberships}}}};
    command_test::ScopedTempDir temporary("joint-local-bundle"); const auto input=temporary.path()/"input.json";
    const auto write=[&] {std::ofstream stream(input); stream<<j::serialize(bundle);};
    const auto read=[](const auto & path) {
        std::ifstream stream(path); j::parse_options options; options.numbers=j::number_precision::precise;
        return j::parse(std::string(std::istreambuf_iterator<char>(stream),{}),{},options);
    };
    write(); p::ComponentLocalBundleRerun(input.string(),(temporary.path()/"fit").string());
    auto expected=p::FitComponent(view,f.y,Vector::Constant(2,.55),parent);
    expected["dataset"]="unit"; expected["case"]="unit"; expected["initial_b"]=j::array{.55};
    expected["observation_snapshot_sha256"]=parent.snapshot_hash;
    EXPECT_EQ(Science(expected),Science(read(temporary.path()/"fit/fit.json")));
    EXPECT_EQ(read(temporary.path()/"fit/audit/scope.json").at("state").at("beta"),expected.at("last_trusted_state").at("beta"));
    bundle["parent_observations"].as_array()[0]=1e6; write();
    EXPECT_THROW(p::ComponentLocalBundleRerun(input.string(),(temporary.path()/"invalid-parent").string()),std::invalid_argument);
    EXPECT_FALSE(std::filesystem::exists(temporary.path()/"invalid-parent"));
    bundle["parent_observations"]=observations; bundle["component_initial_b"]=j::array{0}; write();
    EXPECT_THROW(p::ComponentLocalBundleRerun(input.string(),(temporary.path()/"invalid-width").string()),std::invalid_argument);
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
    for(const char * key:{"runtime_checks","runtime_convergence","last_trusted_state"})
        EXPECT_EQ(before.fit.at(key),fit.at(key));
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
