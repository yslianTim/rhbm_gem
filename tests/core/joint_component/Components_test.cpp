#include <gtest/gtest.h>
#include "support/JointComponentChecks.hpp"
#include "support/CommandTestHelpers.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include <algorithm>
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
    TwoBlocks(double amplitude=2,double noise=.001)
    {
        beta<<2,.2,amplitude,-.15;
        auto support=domain.CopySupport();
        for(int block=0;block<2;++block) for(int k=0;k<40;++k)
        {
            const auto row=block*40+k; const double square=.003*k*k;
            support[static_cast<std::size_t>(block)].push_back({row,square});
            const auto b=second_stage_test::matched::EvaluateBasis(square,.5,2.5);
            y(row)=beta(2*block)*b.gaussian+beta(2*block+1)*b.charge+noise*std::sin(k);
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
TEST(JointComponentChecksTest, StructuralGraphPreservesBridgeChainConstantAndUnobserved)
{
    const p::Domain domain(6,{{{0,0},{1,6.25}},{{1,1},{2,1}},{{2,1}},{{3,1}}, {}});
    const auto c=p::MakeContext(Vector::Ones(6),5); const auto part=p::BuildPartition(domain,{"a","b","c","d","e"});
    ASSERT_EQ(part.components.size(),3); EXPECT_EQ(part.components[0].atoms,(std::vector<Eigen::Index>{0,1,2}));
    EXPECT_EQ(part.constant_rows,(std::vector<Eigen::Index>{4,5})); EXPECT_EQ(part.unobserved_atoms,(std::vector<Eigen::Index>{4}));
    EXPECT_EQ(part.components[0].domain.atoms[0][1].square,6.25);
    EXPECT_EQ(p::Census(domain,part,c).at("constant_objective"),1.0);
}

TEST(JointComponentChecksTest, SharedParentScaleAndComponentRankPolicy)
{
    const p::Domain domain(4,{{{0,0},{1,1}},{{2,0},{3,1}}});
    auto c=p::MakeContext((Vector(4)<<100,200,.01,.02).finished(),2);
    const auto part=p::BuildPartition(domain,c.atom_ids);
    const auto child=p::ComponentContext(c,part.components[1],true);
    EXPECT_EQ(child.scale,c.scale); EXPECT_EQ(child.observations,c.observations);
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

TEST(JointComponentChecksTest, ImmutableComponentViewsMapDomainsWithoutCopyingProblemData)
{
    using rhbm_gem::core::JointProblem;
    using rhbm_gem::core::JointProblemInput;
    namespace n=rhbm_gem::core::joint_component;
    const auto component_layout=[](const auto & parent,const auto & view) {
        rhbm_gem::JointParameterLayout out;
        for(const auto atom:parent.full_atoms) if(view.LocalAtom(static_cast<Eigen::Index>(atom))>=0) out.full_atoms.push_back(atom);
        for(const auto row:parent.informative_rows) if(view.LocalRow(static_cast<Eigen::Index>(row))>=0) out.informative_rows.push_back(row);
        return out;
    };
    const auto check=[&](JointProblem problem) {
        const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
        for(const auto & view:data.partition.components)
        {
            const auto layout=component_layout(data.layout,view);
            const auto prepared=n::MakeComponentProblemView(layout,view);
            ASSERT_EQ(prepared.local_to_parent_atoms.size(),prepared.local_layout.full_atoms.size());
            ASSERT_EQ(prepared.local_to_parent_rows.size(),prepared.local_layout.informative_rows.size());
            EXPECT_EQ(prepared.domain.atoms.Storage(),data.input.get());
            EXPECT_EQ(prepared.domain.atoms.size(),prepared.local_layout.full_atoms.size());
            EXPECT_EQ(prepared.domain.rows,static_cast<Eigen::Index>(prepared.local_layout.informative_rows.size()));
            for(std::size_t local=0;local<prepared.local_to_parent_atoms.size();++local)
            {
                const auto parent=prepared.local_to_parent_atoms[local];
                EXPECT_EQ(prepared.local_layout.full_atoms[local],local);
                EXPECT_GE(view.LocalAtom(parent),0);
                EXPECT_EQ(prepared.domain.atoms.ParentAtom(local),static_cast<std::size_t>(parent));
                EXPECT_EQ(prepared.atom_ids[local],data.input->atom_ids[static_cast<std::size_t>(parent)]);
                const auto supports=prepared.domain.atoms[local];
                for(std::size_t k=0;k<supports.size();++k)
                {
                    const auto support=supports[k];
                    const auto parent_row=prepared.local_to_parent_rows.at(static_cast<std::size_t>(support.row));
                    EXPECT_EQ(supports.ParentRow(k),parent_row);
                    const auto & source=data.input->support.at(static_cast<std::size_t>(parent));
                    EXPECT_TRUE(std::any_of(source.begin(),source.end(),[&](const auto & item) {
                        return item.row==static_cast<std::size_t>(parent_row) && item.squared_distance==support.square;
                    }));
                }
            }
            for(std::size_t local=0;local<prepared.local_to_parent_rows.size();++local)
            {
                const auto parent=prepared.local_to_parent_rows[local];
                EXPECT_EQ(prepared.local_layout.informative_rows[local],local);
                EXPECT_GE(view.LocalRow(parent),0);
                EXPECT_EQ(prepared.row_ids[local],data.input->row_ids[static_cast<std::size_t>(parent)]);
            }
            for(const auto atom:view.atoms)
            {
                const auto found=std::find(prepared.local_to_parent_atoms.begin(),prepared.local_to_parent_atoms.end(),atom);
                const bool included=std::find(layout.full_atoms.begin(),layout.full_atoms.end(),static_cast<std::size_t>(atom))!=layout.full_atoms.end();
                EXPECT_EQ(found!=prepared.local_to_parent_atoms.end(),included);
            }
            for(const auto row:view.rows)
            {
                const auto found=std::find(prepared.local_to_parent_rows.begin(),prepared.local_to_parent_rows.end(),row);
                const bool included=std::find(layout.informative_rows.begin(),layout.informative_rows.end(),static_cast<std::size_t>(row))!=layout.informative_rows.end();
                EXPECT_EQ(found!=prepared.local_to_parent_rows.end(),included);
            }
        }
    };

    JointProblemInput disconnected;
    disconnected.atom_ids={"a","b","c"}; disconnected.row_ids={"r0","r1","r2"};
    disconnected.observations={1,2,3}; disconnected.support={{{0,0}},{{1,1}},{{2,2}}};
    check(JointProblem(std::move(disconnected)));

    JointProblemInput observable;
    observable.atom_ids={"target","halo","other"}; observable.row_ids={"r0","r1","r2"};
    observable.observations={1,2,3}; observable.support={{{0,0},{1,1}},{{0,1}},{{2,2}}};
    observable.selection_domain=rhbm_gem::JointSelectionDomain{};
    observable.selection_domain->target_indices={0,2};
    check(JointProblem(std::move(observable)));
}

TEST(JointComponentChecksTest, NestedCompactRowsPreserveNoncontiguousAndZeroValuedMemberships)
{
    using rhbm_gem::core::JointProblem;
    using rhbm_gem::core::JointProblemInput;
    namespace n=rhbm_gem::core::joint_component;
    JointProblemInput input;
    input.atom_ids={"a","b","c"}; input.row_ids={"r0","r1","r2","r3","r4","r5"};
    input.observations={1,0,-1,0,2,0};
    input.support={{{0,0},{4,1}},{{2,0},{4,2}},{{2,3}}};
    const JointProblem problem(std::move(input));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    ASSERT_EQ(data.partition.components.size(),1u);
    const auto & view=data.partition.components.front();
    EXPECT_EQ(view.rows,(std::vector<Eigen::Index>{0,2,4}));
    rhbm_gem::JointParameterLayout component_layout;
    for(const auto atom:data.layout.full_atoms) if(view.LocalAtom(static_cast<Eigen::Index>(atom))>=0)
        component_layout.full_atoms.push_back(atom);
    for(const auto row:data.layout.informative_rows) if(view.LocalRow(static_cast<Eigen::Index>(row))>=0)
        component_layout.informative_rows.push_back(row);
    const auto prepared=n::MakeComponentProblemView(component_layout,view);
    ASSERT_EQ(prepared.domain.rows,3);
    ASSERT_EQ(prepared.local_to_parent_rows,(n::Indices{0,2,4}));
    for(std::size_t atom=0;atom<prepared.domain.atoms.size();++atom)
    {
        const auto rows=prepared.domain.atoms[atom];
        const auto root_atom=prepared.local_to_parent_atoms[atom];
        const auto & root_support=data.input->support[static_cast<std::size_t>(root_atom)];
        ASSERT_EQ(rows.size(),root_support.size());
        for(std::size_t k=0;k<rows.size();++k)
        {
            const auto row=rows[k];
            EXPECT_EQ(rows.ParentRow(k),static_cast<Eigen::Index>(root_support[k].row));
            EXPECT_EQ(prepared.local_to_parent_rows[static_cast<std::size_t>(row.row)],
                static_cast<Eigen::Index>(root_support[k].row));
            EXPECT_DOUBLE_EQ(row.square,root_support[k].squared_distance);
        }
    }
    EXPECT_EQ(prepared.domain.atoms[0][0].row,0);
    EXPECT_EQ(prepared.domain.atoms[0][1].row,2);
    EXPECT_EQ(prepared.domain.atoms[1][0].row,1);
    EXPECT_LT(data.y(2),0);
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
    const auto plus=p::EvaluatePartitioned(f.domain,f.y,f.eta+h*direction,part,c);
    const auto minus=p::EvaluatePartitioned(f.domain,f.y,f.eta-h*direction,part,c);
    ASSERT_TRUE(plus.valid && minus.valid);
    EXPECT_LT(((plus.residual-minus.residual)/(2*h*c.scale)-d.jacobian*direction).norm(),1e-8);
}

TEST(JointComponentChecksTest, FailureIsolationDoesNotFillMissingRowsWithZero)
{
    TwoBlocks f(2,0); auto support=f.domain.CopySupport(); support.push_back(support[1]); f.domain=p::Domain(f.domain.rows,std::move(support));
    const auto input=p::MakeSyntheticJointProblemInput(f.domain,f.y); f.domain=p::Domain(input);
    auto c=p::MakeContext(input); const auto part=p::BuildPartition(f.domain,c.atom_ids);
    const Vector initial=Vector::Constant(3,.55);
    const auto result=p::FitComponents(input,initial,part,c);
    ASSERT_EQ(result.at("components").as_array().size(),2);
    EXPECT_TRUE(result.at("components").at(0).at("search_success").as_bool());
    EXPECT_FALSE(result.at("components").at(1).at("usable_state").as_bool());
    EXPECT_FALSE(result.at("prediction_available").as_bool()); EXPECT_TRUE(result.at("prediction").is_null());
    EXPECT_NE(result.at("runtime_convergence"),"passed");
    for(std::size_t r=0;r<81;++r) EXPECT_EQ(result.at("available_row_mask").at(r).as_bool(),r<40 || r==80);
    const auto failed_first=p::FitComponent(part.components[1],initial,c);
    const auto healthy_second=p::FitComponent(part.components[0],initial,c);
    EXPECT_EQ(Science(failed_first),Science(result.at("components").at(1)));
    EXPECT_EQ(Science(healthy_second),Science(result.at("components").at(0)));
}

TEST(JointComponentChecksTest, UnobservedAtomsKeepRawPredictionButDisableProfileClaims)
{
    TwoBlocks f(2,0); auto support=f.domain.CopySupport(); support[1].clear(); f.domain=p::Domain(f.domain.rows,std::move(support));
    const auto input=p::MakeSyntheticJointProblemInput(f.domain,f.y); f.domain=p::Domain(input); const auto c=p::MakeContext(input);
    const auto part=p::BuildPartition(f.domain,c.atom_ids);
    const auto state=p::SameState(f.domain,f.y,f.eta,f.beta,part,c);
    EXPECT_TRUE(state.at("raw").at("passed").as_bool()); EXPECT_TRUE(state.at("passed").as_bool());
    EXPECT_FALSE(state.at("full_equivalence").as_bool());
    EXPECT_FALSE(state.at("profile_solves").at(0).at("available").as_bool());
    EXPECT_EQ(part.unobserved_atoms,(std::vector<Eigen::Index>{1}));
    const auto fit=p::FitComponents(input,Vector(f.eta.array().exp()),part,c);
    EXPECT_TRUE(fit.at("components").at(0).at("search_success").as_bool());
    EXPECT_FALSE(fit.at("components").at(1).at("usable_state").as_bool());
    EXPECT_FALSE(fit.at("prediction_available").as_bool());
    RecordProperty("fit_availability",boost::json::serialize(fit));
}
