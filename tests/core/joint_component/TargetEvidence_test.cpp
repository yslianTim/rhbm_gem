#include <gtest/gtest.h>
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/TargetEvidence.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>

namespace {
namespace c=rhbm_gem::core;
namespace n=c::joint_component;
using Status=rhbm_gem::JointCheckStatus;

c::JointProblem MakeProblem(std::size_t atoms)
{
    c::JointProblemInput input; input.observations={0.0}; input.row_ids={"row-0"};
    input.support.resize(atoms); input.atom_ids.reserve(atoms);
    input.selection_domain.emplace(); input.selection_domain->target_indices.resize(atoms);
    std::iota(input.selection_domain->target_indices.begin(),input.selection_domain->target_indices.end(),0);
    for(std::size_t atom=0;atom<atoms;++atom)
    {
        input.atom_ids.push_back("atom-"+std::to_string(atom));
        input.support[atom].push_back({0,.1});
    }
    return c::JointProblem(std::move(input));
}

rhbm_gem::JointState MakeState(std::size_t atoms)
{
    rhbm_gem::JointState state;
    state.ac.resize(2*atoms); state.b.assign(atoms,.7); state.log_b.assign(atoms,std::log(.7));
    for(std::size_t atom=0;atom<atoms;++atom)
    {state.ac[2*atom]=.5; state.ac[2*atom+1]=.1;}
    state.objective=0.;
    return state;
}

rhbm_gem::JointBlockConditioning ConditioningBlock(std::string id,std::vector<std::size_t> atoms)
{
    rhbm_gem::JointBlockConditioning block; block.block_id=std::move(id); block.status="assessed";
    block.conditioning_classification="WellConditioned"; block.core_atoms=std::move(atoms);
    return block;
}

rhbm_gem::JointBlockParameterStability StabilityBlock(std::string id,std::vector<std::size_t> atoms)
{
    rhbm_gem::JointBlockParameterStability block; block.block_id=std::move(id); block.status="Stable";
    block.core_atoms=std::move(atoms); return block;
}

TEST(JointTargetEvidenceTest, DuplicateMembershipAndBlockIdentityAreExplicitWarnings)
{
    const auto problem=MakeProblem(3); const auto state=MakeState(3);
    const std::vector<rhbm_gem::JointBlockConditioning> conditioning{
        ConditioningBlock("block-a",{0,1}),ConditioningBlock("block-b",{0,2})};
    const auto evidence=n::AssessTargetPointState(problem,problem.ParameterLayout(),state,Status::Passed,
        conditioning,{},rhbm_gem::JointEvidenceScope::AssembledGlobal,1);
    ASSERT_EQ(evidence.point_atoms.size(),3u);
    EXPECT_EQ(evidence.point_atoms[0].point_quality,"Warning");
    EXPECT_EQ(evidence.point_atoms[0].reason,"invalid-diagnostic-block-mapping");
    const auto mapping=std::find_if(evidence.checks.begin(),evidence.checks.end(),[](const auto & check) {
        return check.name=="diagnostic-block-mapping";
    });
    ASSERT_NE(mapping,evidence.checks.end());
    EXPECT_EQ(mapping->status,Status::Unavailable);
    EXPECT_EQ(mapping->reason,"duplicate-diagnostic-core-membership");
    EXPECT_EQ(n::TargetConvergence(evidence,true),Status::Passed);

    const std::vector<rhbm_gem::JointBlockConditioning> duplicate_identity{
        ConditioningBlock("block-a",{0}),ConditioningBlock("block-a",{1,2})};
    const auto invalid_id=n::AssessTargetPointState(problem,problem.ParameterLayout(),state,Status::Passed,
        duplicate_identity,{},rhbm_gem::JointEvidenceScope::AssembledGlobal,1);
    EXPECT_EQ(invalid_id.point_atoms[0].local_conditioning,"Unavailable");
    const auto identity_check=std::find_if(invalid_id.checks.begin(),invalid_id.checks.end(),[](const auto & check) {
        return check.name=="diagnostic-block-mapping";
    });
    ASSERT_NE(identity_check,invalid_id.checks.end());
    EXPECT_EQ(identity_check->reason,"duplicate-diagnostic-block-identity");
}

TEST(JointTargetEvidenceTest, MissingMembershipIsNotSilentlyIgnored)
{
    const auto problem=MakeProblem(3); const auto state=MakeState(3);
    const std::vector<rhbm_gem::JointBlockConditioning> conditioning{ConditioningBlock("block-a",{0,1})};
    const auto evidence=n::AssessTargetPointState(problem,problem.ParameterLayout(),state,Status::Passed,
        conditioning,{},rhbm_gem::JointEvidenceScope::AssembledGlobal,1);
    EXPECT_EQ(evidence.point_atoms[2].local_conditioning,"Unavailable");
    EXPECT_EQ(evidence.point_atoms[2].reason,"invalid-diagnostic-block-mapping");
    const auto mapping=std::find_if(evidence.checks.begin(),evidence.checks.end(),[](const auto & check) {
        return check.name=="diagnostic-block-mapping";
    });
    ASSERT_NE(mapping,evidence.checks.end());
    EXPECT_EQ(mapping->reason,"missing-diagnostic-core-membership");
}

TEST(JointTargetEvidenceTest, ComponentLocalAndAssembledEvidenceAgreeForNoncontiguousAtoms)
{
    c::JointProblemInput input;
    input.atom_ids={"a","b","c","d"}; input.row_ids={"r0","r1","r2"};
    input.observations={0.,0.,0.};
    input.support={{{0,.1}},{{1,.1}},{{2,.1}},{{1,.2}}};
    input.selection_domain.emplace(); input.selection_domain->target_indices={1,3};
    const c::JointProblem problem(std::move(input));
    const auto & data=c::JointProblemAccess::Get(problem);
    const auto view=std::find_if(data.partition.components.begin(),data.partition.components.end(),[](const auto & component) {
        return component.LocalAtom(1)>=0 && component.LocalAtom(3)>=0;
    });
    ASSERT_NE(view,data.partition.components.end());
    rhbm_gem::JointParameterLayout local_layout;
    for(const auto atom:data.layout.full_atoms) if(view->LocalAtom(static_cast<Eigen::Index>(atom))>=0)
        local_layout.full_atoms.push_back(atom);
    for(const auto row:data.layout.informative_rows) if(view->LocalRow(static_cast<Eigen::Index>(row))>=0)
        local_layout.informative_rows.push_back(row);
    ASSERT_EQ(local_layout.full_atoms,(std::vector<std::size_t>{1,3}));

    const auto local_state=MakeState(local_layout.full_atoms.size());
    const std::vector<rhbm_gem::JointBlockConditioning> local_conditioning{
        ConditioningBlock("component/block-1",{1,3})};
    const std::vector<rhbm_gem::JointBlockParameterStability> local_stability{
        StabilityBlock("component/block-1",{1,3})};
    const auto local=n::AssessTargetPointState(problem,local_layout,local_state,Status::Passed,local_conditioning,
        local_stability,rhbm_gem::JointEvidenceScope::ComponentLocal,view->rows.size());

    const auto & global_layout=problem.ParameterLayout(); const auto global_state=MakeState(global_layout.full_atoms.size());
    std::vector<rhbm_gem::JointBlockConditioning> global_conditioning;
    std::vector<rhbm_gem::JointBlockParameterStability> global_stability;
    for(const auto atom:global_layout.full_atoms)
    {
        global_conditioning.push_back(ConditioningBlock("global/"+std::to_string(atom),{atom}));
        global_stability.push_back(StabilityBlock("global/"+std::to_string(atom),{atom}));
    }
    const auto global=n::AssessTargetPointState(problem,global_layout,global_state,Status::Passed,global_conditioning,
        global_stability,rhbm_gem::JointEvidenceScope::AssembledGlobal,problem.Input().observations.size());
    ASSERT_EQ(local.point_atoms.size(),2u);
    for(std::size_t index=0;index<local.point_atoms.size();++index)
    {
        const auto found=std::find_if(global.point_atoms.begin(),global.point_atoms.end(),[&](const auto & point) {
            return point.atom==local.point_atoms[index].atom;
        });
        ASSERT_NE(found,global.point_atoms.end());
        EXPECT_EQ(local.point_atoms[index].local_conditioning,found->local_conditioning);
        EXPECT_EQ(local.point_atoms[index].parameter_stability,found->parameter_stability);
        EXPECT_EQ(local.point_atoms[index].point_quality,found->point_quality);
    }
    EXPECT_EQ(n::TargetConvergence(local,true),Status::Passed);
    EXPECT_EQ(n::TargetConvergence(global,true),Status::Passed);
}

#ifdef RHBM_GEM_TEST_INSTRUMENTATION
TEST(JointTargetEvidenceTest, LargeComponentMappingWorkIsLinearInPointsAndMemberships)
{
    constexpr std::size_t atom_count=37'406,core_size=12;
    const auto problem=MakeProblem(atom_count); const auto state=MakeState(atom_count);
    std::vector<rhbm_gem::JointBlockConditioning> conditioning;
    std::vector<rhbm_gem::JointBlockParameterStability> stability;
    for(std::size_t first=0;first<atom_count;first+=core_size)
    {
        std::vector<std::size_t> atoms;
        const auto last=std::min(atom_count,first+core_size);
        for(std::size_t atom=first;atom<last;++atom) atoms.push_back(atom);
        const auto identity=std::to_string(first/core_size);
        conditioning.push_back(ConditioningBlock("conditioning/"+identity,atoms));
        stability.push_back(StabilityBlock("stability/"+identity,std::move(atoms)));
    }
    const auto & data=c::JointProblemAccess::Get(problem);
    ASSERT_EQ(data.partition.components.size(),1u);
    const auto & local_component=data.partition.components.front();
    auto & work=n::TargetEvidenceWorkForTesting(); work={};
    auto started=std::chrono::steady_clock::now();
    const auto local=n::AssessTargetPointState(problem,problem.ParameterLayout(),state,Status::Passed,
        conditioning,stability,rhbm_gem::JointEvidenceScope::ComponentLocal,local_component.rows.size());
    const double local_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    ASSERT_EQ(local.point_atoms.size(),atom_count);
    EXPECT_EQ(work.memberships_indexed,4*atom_count);
    EXPECT_EQ(work.atom_lookups,3*atom_count);
    const auto local_work=work;

    work={};
    started=std::chrono::steady_clock::now();
    const auto global=n::AssessTargetPointState(problem,problem.ParameterLayout(),state,Status::Passed,
        conditioning,stability,rhbm_gem::JointEvidenceScope::AssembledGlobal,problem.Input().observations.size());
    const double global_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    ASSERT_EQ(global.point_atoms.size(),atom_count);
    EXPECT_EQ(work.memberships_indexed,4*atom_count);
    EXPECT_EQ(work.atom_lookups,3*atom_count);
    EXPECT_EQ(global.point_atoms.back().point_quality,"Available");
    ASSERT_EQ(local.point_atoms.size(),global.point_atoms.size());
    for(std::size_t atom=0;atom<atom_count;++atom)
    {
        EXPECT_EQ(local.point_atoms[atom].local_conditioning,global.point_atoms[atom].local_conditioning);
        EXPECT_EQ(local.point_atoms[atom].parameter_stability,global.point_atoms[atom].parameter_stability);
        EXPECT_EQ(local.point_atoms[atom].point_quality,global.point_atoms[atom].point_quality);
    }
    EXPECT_EQ(local_work.memberships_indexed,work.memberships_indexed);
    EXPECT_EQ(local_work.atom_lookups,work.atom_lookups);
    std::cout<<"large-target-evidence-points="<<atom_count<<" blocks="<<conditioning.size()
        <<" observations="<<problem.Input().observations.size()
        <<" component-local-seconds="<<local_seconds<<" assembled-global-seconds="<<global_seconds<<"\n";
}
#endif
}
