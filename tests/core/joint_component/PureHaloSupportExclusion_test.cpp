#include <gtest/gtest.h>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include "support/CommandTestHelpers.hpp"
#include <boost/json.hpp>

#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/gaussian_fit/PreparedLocalGaussianFit.hpp"
#include "core/detail/MapInterpolation.hpp"
#include "core/detail/joint_component/HaloMedianReconstruction.hpp"
#include "core/detail/GroupPotentialFitting.hpp"
#include "core/detail/joint_component/PostFitPeeling.hpp"
#include "support/JointPostprocessingTestHelpers.hpp"
#include "core/detail/PotentialFittingWorkflow.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/GaussSeidelBcd.hpp"
#include "core/detail/joint_component/WidthCalibration.hpp"
#include "core/command/detail/SimulationGeometry.hpp"
#include "data/detail/AtomClassifier.hpp"
#include "data/detail/JointStageAdapter.hpp"
#include "data/io/detail/JointResultJson.hpp"
#include "data/io/detail/StageResultJson.hpp"
#include "core/detail/joint_component/JointComponentEstimator.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

#include <Eigen/QR>
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace rhbm_gem;
using namespace rhbm_gem::core;
namespace core=rhbm_gem::core;
namespace postprocessing_test=second_stage_test::postprocessing;

JointProblemInput PureInput()
{
    JointProblemInput input;
    input.atom_ids={"1"};
    input.row_ids={"0","1","2","3","4","5","6","7","8"};
    input.support.resize(1);
    for (std::size_t row=0;row<input.row_ids.size();++row)
    {
        const auto basis=core::joint_component::EvaluateKernel(.04*static_cast<double>(row*row),.65,2.5);
        input.support[0].push_back({row,.04*static_cast<double>(row*row)});
        input.observations.push_back(2.0*basis.gaussian-.2*basis.charge);
    }
    input.selection_domain.emplace();
    input.selection_domain->contributor_policy=std::string(PureHaloSupportExclusionContributorPolicy);
    input.selection_domain->target_indices={0};
    input.selection_domain->pure_halo_support_exclusion=JointPureHaloSupportExclusionProvenance{
        input.observations.size(),input.observations.size(),0,0,0,0};
    return input;
}

JointInitialization TargetInitialization()
{
    JointInitialization initialization;
    initialization.valid=true;
    initialization.tau={.65};
    JointInitializationAtom atom;
    atom.id="1";
    atom.target=true;
    atom.reason="valid-width";
    atom.used_tau=.65;
    atom.original_tau=.65;
    initialization.atoms.push_back(std::move(atom));
    initialization.reason="valid-widths";
    return initialization;
}

JointAnalysisResult PurePeelingSnapshot(const JointProblem & problem,const ModelObject & model)
{
    (void)model;
    JointAnalysisResult result;
    result.parameterization_contract=std::string(JointPureHaloSupportExclusionParameterizationContract);
    result.atom_ids=problem.OriginalInput().atom_ids;
    result.row_ids=problem.OriginalInput().row_ids;
    result.fitting_row_ids=problem.Input().row_ids;
    result.fitting_to_original_contributors=problem.FittingToOriginalContributors();
    result.fitting_to_original_rows=problem.FittingToOriginalRows();
    result.selection_domain=problem.OriginalInput().selection_domain;
    result.available_row_mask.assign(result.row_ids.size(),false);
    for (const auto row:problem.FittingToOriginalRows()) result.available_row_mask[row]=true;
    result.runtime_convergence=JointCheckStatus::Passed;
    result.target_runtime_convergence=JointCheckStatus::Passed;
    JointAnalysisComponent component;
    component.id="pure-component";
    component.atoms={0};
    component.rows=problem.FittingToOriginalRows();
    component.runtime_convergence=JointCheckStatus::Passed;
    component.target_runtime_convergence=JointCheckStatus::Passed;
    component.state=JointState{{1.0,0.0},{1.0},{0.0},{0.0},0.0};
    result.components.push_back(component);
    result.assembled_state=component.state;
    JointHaloMedianSnapshot medians;
    medians.joint_run_id="pure-peeling-run";
    medians.valid_target_count=1;
    medians.global_median={.5,.8,0.0};
    if (problem.OriginalInput().atom_ids.size()>
        problem.OriginalInput().selection_domain->target_indices.size())
    {
        JointHaloMedianReconstruction reconstruction;
        reconstruction.halo_contributor_index=1;
        reconstruction.atom_id="2";
        reconstruction.amplitude=.5;
        reconstruction.width=.8;
        reconstruction.charge=0.0;
        reconstruction.source=HaloMedianSource::GeometryLocalTargetMedian;
        reconstruction.donor_count=1;
        reconstruction.joint_run_id="pure-peeling-run";
        reconstruction.donor_atom_ids={"1"};
        medians.halos.push_back(std::move(reconstruction));
        result.halo_median_reconstruction=std::move(medians);
    }
    return result;
}

struct MapFixture
{
    ModelObject model;
    MapObject map;
};

MapFixture BoundaryFixture()
{
    auto target=std::make_unique<AtomObject>();
    target->SetSerialID(1);
    target->SetElement(Element::CARBON);
    target->SetPosition(0,0,0);
    auto halo=std::make_unique<AtomObject>();
    halo->SetSerialID(2);
    halo->SetElement(Element::CARBON);
    halo->SetPosition(5,0,0);
    std::vector<std::unique_ptr<AtomObject>> atoms;
    atoms.push_back(std::move(target));
    atoms.push_back(std::move(halo));
    ModelObject model(std::move(atoms));
    model.SelectAtoms([](const auto & atom){return atom.GetSerialID()==1;});
    auto values=std::make_unique<double[]>(11);
    for (std::size_t row=0;row<11;++row) values[row]=static_cast<double>(row+1);
    MapObject map({11,1,1},{.5,.5,.5},{-2.5,0,0},std::move(values));
    return {std::move(model),std::move(map)};
}

MapFixture MultipleHaloBoundaryFixture()
{
    std::vector<std::unique_ptr<AtomObject>> atoms;
    const auto add_atom=[&](int serial,double x) {
        auto atom=std::make_unique<AtomObject>();
        atom->SetSerialID(serial);
        atom->SetElement(Element::CARBON);
        atom->SetPosition(x,0,0);
        atoms.push_back(std::move(atom));
    };
    add_atom(30,0.0);
    add_atom(10,4.5);
    add_atom(20,5.0);
    add_atom(40,20.0);
    ModelObject model(std::move(atoms));
    model.SelectAtoms([](const auto & atom){return atom.GetSerialID()==30;});
    auto values=std::make_unique<double[]>(11);
    for(std::size_t row=0;row<11;++row) values[row]=static_cast<double>(row+1);
    MapObject map({11,1,1},{.5,.5,.5},{-2.5,0,0},std::move(values));
    return {std::move(model),std::move(map)};
}

MapFixture WorkflowFixture(bool qualified_donors=false,bool coupled_targets=false)
{
    std::vector<std::unique_ptr<AtomObject>> atoms;
    const auto add_atom=[&](int serial,const std::array<double,3> & position,int component,int atom_key,
        const char * atom_id) {
        auto atom=std::make_unique<AtomObject>();
        atom->SetSerialID(serial);
        atom->SetElement(Element::CARBON);
        atom->SetComponentKey(static_cast<ComponentKey>(component));
        atom->SetAtomKey(static_cast<AtomKey>(atom_key));
        atom->SetChainID("A");
        atom->SetComponentID("ALA");
        atom->SetAtomID(atom_id);
        atom->SetPosition(position);
        atoms.push_back(std::move(atom));
    };
    struct Profile
    {
        std::array<double,3> position;
        double amplitude,width,charge;
    };
    std::vector<Profile> profiles{
        {{0,0,0},2.0,.65,-.2},{{4.5,0,0},.4,.8,.1},{{4.8,0,0},.3,.9,-.1}};
    if(coupled_targets) profiles[0]={{0,0,0},2.0,.5,.2};
    add_atom(1,profiles[0].position,1,7,"CA");
    add_atom(2,profiles[1].position,1,7,"CB");
    add_atom(3,profiles[2].position,3,8,"CG");
    if(qualified_donors)
    {
        std::vector<std::array<double,3>> candidates;
        for(double x:{-6.,0.,6.}) for(double y:{-6.,0.,6.}) for(double z:{-6.,0.,6.})
        {
            const std::array<double,3> position{x,y,z};
            if(position==std::array<double,3>{0,0,0} || position==std::array<double,3>{6,0,0}) continue;
            candidates.push_back(position);
        }
        for(std::size_t donor=0;donor<19;++donor)
        {
            profiles.push_back({candidates[donor],2.0,.65,-.2});
            add_atom(static_cast<int>(4+donor),candidates[donor],10+static_cast<int>(donor),
                100+static_cast<int>(donor),"CA");
        }
    }
    if(coupled_targets)
    {
        profiles.push_back({{1,0,0},2.0,.5,.2});
        add_atom(23,profiles.back().position,23,123,"CA");
    }
    ModelObject model(std::move(atoms));
    model.SelectAtoms([qualified_donors](const auto & atom) {
        return atom.GetSerialID()==1 || (qualified_donors && atom.GetSerialID()>=4);
    });

    const std::array<int,3> dimensions=qualified_donors ? std::array<int,3>{61,61,61} :
        std::array<int,3>{45,45,45};
    constexpr std::array<double,3> spacing{.3,.3,.3};
    const std::array<double,3> origin=qualified_donors ? std::array<double,3>{-9.,-9.,-9.} :
        std::array<double,3>{-6.6,-6.6,-6.6};
    auto values=std::make_unique<double[]>(static_cast<std::size_t>(dimensions[0])*
        static_cast<std::size_t>(dimensions[1])*static_cast<std::size_t>(dimensions[2]));
    for (int z=0;z<dimensions[2];++z) for (int y=0;y<dimensions[1];++y) for (int x=0;x<dimensions[0];++x)
    {
        const auto position=core::simulation::GridPosition({x,y,z},spacing,origin);
        double value{};
        for (const auto & profile:profiles)
        {
            const auto square=core::simulation::SupportSquare(position,profile.position);
            if (square>6.25) continue;
            const auto basis=core::joint_component::EvaluateKernel(square,profile.width,2.5);
            value+=profile.amplitude*basis.gaussian+profile.charge*basis.charge;
        }
        const auto index=static_cast<std::size_t>(x)+static_cast<std::size_t>(dimensions[0]) *
            (static_cast<std::size_t>(y)+static_cast<std::size_t>(dimensions[1])*static_cast<std::size_t>(z));
        values[index]=value;
    }
    MapObject map(dimensions,spacing,origin,std::move(values));
    return {std::move(model),std::move(map)};
}

MapFixture MultiTargetHaloWorkflowFixture()
{
    std::vector<std::unique_ptr<AtomObject>> atoms;
    const auto add_atom=[&](int serial,double x,int component,int atom_key,
        const char * atom_id) {
        auto atom=std::make_unique<AtomObject>();
        atom->SetSerialID(serial);
        atom->SetElement(Element::CARBON);
        atom->SetComponentKey(static_cast<ComponentKey>(component));
        atom->SetAtomKey(static_cast<AtomKey>(atom_key));
        atom->SetChainID("A");
        atom->SetComponentID("ALA");
        atom->SetAtomID(atom_id);
        atom->SetPosition(x,0,0);
        atoms.push_back(std::move(atom));
    };
    add_atom(1,0.0,1,7,"CA");
    add_atom(2,3.5,1,7,"C");
    add_atom(3,7.0,1,7,"N");
    add_atom(4,10.0,1,7,"O");
    add_atom(5,-4.5,3,8,"CB");
    ModelObject model(std::move(atoms));
    model.SelectAtoms([](const auto & atom){return atom.GetSerialID()<=3;});

    constexpr std::array<int,3> dimensions{70,25,25};
    constexpr std::array<double,3> spacing{.3,.3,.3};
    constexpr std::array<double,3> origin{-7.5,-3.75,-3.75};
    const std::array<std::array<double,3>,5> positions{{
        {0,0,0},{3.5,0,0},{7.0,0,0},{10.0,0,0},{-4.5,0,0}}};
    const std::array<double,5> amplitudes{{2.0,1.6,1.2,.5,.3}};
    const std::array<double,5> widths{{.65,.75,.85,.8,1.0}};
    const std::array<double,5> charges{{-.2,.15,-.1,.3,-.25}};
    auto values=std::make_unique<double[]>(static_cast<std::size_t>(dimensions[0])*dimensions[1]*dimensions[2]);
    for (int z=0;z<dimensions[2];++z) for (int y=0;y<dimensions[1];++y) for (int x=0;x<dimensions[0];++x)
    {
        const auto position=core::simulation::GridPosition({x,y,z},spacing,origin);
        double value{};
        for (std::size_t atom=0;atom<positions.size();++atom)
        {
            const auto square=core::simulation::SupportSquare(position,positions[atom]);
            if (square>6.25) continue;
            const auto basis=core::joint_component::EvaluateKernel(square,widths[atom],2.5);
            value+=amplitudes[atom]*basis.gaussian+charges[atom]*basis.charge;
        }
        const auto index=static_cast<std::size_t>(x)+static_cast<std::size_t>(dimensions[0]) *
            (static_cast<std::size_t>(y)+static_cast<std::size_t>(dimensions[1])*static_cast<std::size_t>(z));
        values[index]=value;
    }
    MapObject map(dimensions,spacing,origin,std::move(values));
    return {std::move(model),std::move(map)};
}

MapFixture ReconstructionOnlyFixture()
{
    constexpr std::array<double,3> target_position{0.0,0.0,0.0};
    constexpr std::array<double,3> other_target_position{-4.5,-4.5,-4.5};
    constexpr std::array<std::array<double,3>,2> halo_positions{{
        {4.72,1.18,1.18},{4.72,-1.18,1.18}}};
    const GaussianModel3D point{2.0,.65,-.2};
    const GaussianModel3D other_target_point{10.0,1.1,.4};
    std::vector<std::unique_ptr<AtomObject>> atoms;
    const auto add_atom=[&](int serial,const std::array<double,3> & position,
        ComponentKey component,AtomKey atom_key) {
        auto atom=std::make_unique<AtomObject>();
        atom->SetSerialID(serial);
        atom->SetElement(Element::CARBON);
        atom->SetComponentKey(component);
        atom->SetAtomKey(atom_key);
        atom->SetChainID("A");
        atom->SetComponentID("ALA");
        atom->SetAtomID("CA");
        atom->SetPosition(position);
        atoms.push_back(std::move(atom));
    };
    add_atom(1,target_position,static_cast<ComponentKey>(1),static_cast<AtomKey>(7));
    add_atom(2,halo_positions[0],static_cast<ComponentKey>(1),static_cast<AtomKey>(7));
    add_atom(3,halo_positions[1],static_cast<ComponentKey>(1),static_cast<AtomKey>(7));
    add_atom(4,other_target_position,static_cast<ComponentKey>(2),static_cast<AtomKey>(9));
    ModelObject model(std::move(atoms));
    model.SelectAtoms([](const auto & atom){return atom.GetSerialID()==1 || atom.GetSerialID()==4;});

    constexpr std::array<int,3> dimensions{45,45,45};
    constexpr std::array<double,3> spacing{.3,.3,.3};
    constexpr std::array<double,3> origin{-6.6,-6.6,-6.6};
    auto values=std::make_unique<double[]>(static_cast<std::size_t>(dimensions[0])*dimensions[1]*dimensions[2]);
    for (int z=0;z<dimensions[2];++z) for (int y=0;y<dimensions[1];++y) for (int x=0;x<dimensions[0];++x)
    {
        const auto position=core::simulation::GridPosition({x,y,z},spacing,origin);
        double value{};
        const auto add_profile=[&](const auto & center,const GaussianModel3D & profile) {
            const double square=core::simulation::SupportSquare(position,center);
            if (square>6.25) return;
            const auto basis=core::joint_component::EvaluateKernel(square,profile.GetWidth(),2.5);
            value+=profile.GetAmplitude()*basis.gaussian+profile.GetOffset()*basis.charge;
        };
        add_profile(target_position,point);
        for (const auto & halo:halo_positions) add_profile(halo,point);
        add_profile(other_target_position,other_target_point);
        const auto index=static_cast<std::size_t>(x)+static_cast<std::size_t>(dimensions[0]) *
            (static_cast<std::size_t>(y)+static_cast<std::size_t>(dimensions[1])*static_cast<std::size_t>(z));
        values[index]=value;
    }
    MapObject map(dimensions,spacing,origin,std::move(values));
    return {std::move(model),std::move(map)};
}

JointAnalysisResult QualifiedReconstructionSnapshot(const JointProblem & problem,const ModelObject & model)
{
    auto result=PurePeelingSnapshot(problem,model);
    result.halo_median_reconstruction.reset();
    result.initialization=TargetInitialization();
    result.search_completed=true;
    const auto targets=problem.OriginalInput().selection_domain->target_indices;
    result.components.front().atoms={targets.front()};
    result.components.front().state->ac={2.0,-.2};
    result.components.front().state->tau={.65};
    result.components.front().search_completed=true;
    result.components.front().stop_reason="converged";
    result.assembled_state=result.components.front().state;
    result.objective=0.0;
    if (targets.size()>1)
    {
        JointAnalysisComponent component;
        component.id="second-target-component";
        component.atoms={targets[1]};
        component.rows=problem.FittingToOriginalRows();
        component.runtime_convergence=JointCheckStatus::Passed;
        component.target_runtime_convergence=JointCheckStatus::Passed;
        component.search_completed=true;
        component.stop_reason="converged";
        component.state=JointState{{10.0,.4},{1.1},{0.0},{0.0},0.0};
        result.components.push_back(std::move(component));
    }
    result.target_evidence.emplace();
    result.target_evidence->contract=std::string(JointTargetPointValidationContractVersion);
    result.target_evidence->original_rows=problem.OriginalInput().row_ids.size();
    for (const auto target:targets)
    {
        JointAtomPointValidation point;
        point.atom=target;
        point.component_convergence=JointCheckStatus::Passed;
        point.parameter_validity=JointCheckStatus::Passed;
        point.local_conditioning="Adequate";
        point.point_quality="NotAssessed";
        point.reason="test-qualified-point";
        result.target_evidence->point_atoms.push_back(std::move(point));
    }
    auto medians=core::detail::BuildHaloMedianReconstruction(
        model,problem,result,"reconstruction-only-run");
    if (problem.OriginalInput().atom_ids.size()>targets.size())
        result.halo_median_reconstruction=medians.TakeSnapshot();
    return result;
}

;


double ReferenceContribution(double square,const GaussianModel3D & point,bool charge_only=false)
{
    if (square>6.25) return 0.0;
    const double tau=point.GetWidth();
    const double r=std::sqrt(square);
    const double charge=r<1e-5 ? std::sqrt(2.0/std::acos(-1.0))/tau : std::erf(r/(std::sqrt(2.0)*tau))/r;
    const double gaussian=std::pow(2.0*std::acos(-1.0)*tau*tau,-1.5)*std::exp(-square/(2.0*tau*tau));
    return point.GetOffset()*charge+(charge_only ? 0.0 : point.GetAmplitude()*gaussian);
}
}


TEST(PureHaloSupportExclusionTest, ExcludesExactBoundarySupportAndPreservesOriginalDomain)
{
    auto fixture=BoundaryFixture();
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    const auto & fitting=problem.Input();
    const auto & original=problem.OriginalInput();
    ASSERT_TRUE(fitting.selection_domain);
    EXPECT_EQ(fitting.selection_domain->contributor_policy,PureHaloSupportExclusionContributorPolicy);
    EXPECT_EQ(original.row_ids.size(),11u);
    EXPECT_EQ(fitting.row_ids.size(),10u);
    EXPECT_EQ(fitting.row_ids.back(),"9");
    EXPECT_EQ(original.row_ids.back(),"10");
    EXPECT_EQ(problem.FittingToOriginalRows(),(std::vector<std::size_t>{0,1,2,3,4,5,6,7,8,9}));
    EXPECT_EQ(problem.FittingToOriginalContributors(),(std::vector<std::size_t>{0}));
    ASSERT_TRUE(original.selection_domain->pure_halo_support_exclusion);
    const auto & provenance=*original.selection_domain->pure_halo_support_exclusion;
    EXPECT_EQ(provenance.original_rows,11u);
    EXPECT_EQ(provenance.fit_rows,10u);
    EXPECT_EQ(provenance.excluded_rows,1u);
    EXPECT_EQ(provenance.candidate_halos,1u);
    EXPECT_EQ(provenance.retained_halos,1u);
    EXPECT_EQ(provenance.zero_support_halos,0u);
    double norm{};
    for (const double value:original.observations) norm+=value*value;
    EXPECT_DOUBLE_EQ(problem.ObservationScale(),std::max(1.0,std::sqrt(norm)));
    for (const auto & support:fitting.support)
        for (const auto & entry:support) EXPECT_LT(entry.row,fitting.row_ids.size());
}

TEST(PureHaloSupportExclusionTest, BuildsFromSharedGeometryWithOverlappingHaloRowsAndStableMappings)
{
    auto fixture=MultipleHaloBoundaryFixture();
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    const auto & original=problem.OriginalInput();
    const auto & fitting=problem.Input();
    ASSERT_TRUE(original.selection_domain);
    EXPECT_EQ(original.selection_domain->contributor_policy,PureHaloSupportExclusionContributorPolicy);
    EXPECT_EQ(original.atom_ids,(std::vector<std::string>{"10","20","30"}));
    EXPECT_EQ(original.selection_domain->target_indices,(std::vector<std::size_t>{2}));
    EXPECT_EQ(fitting.atom_ids,(std::vector<std::string>{"30"}));
    EXPECT_EQ(fitting.row_ids,(std::vector<std::string>{"0","1","2","3","4","5","6","7","8"}));
    EXPECT_EQ(problem.FittingToOriginalContributors(),(std::vector<std::size_t>{2}));
    EXPECT_EQ(problem.FittingToOriginalRows(),(std::vector<std::size_t>{0,1,2,3,4,5,6,7,8}));
    ASSERT_TRUE(original.selection_domain->pure_halo_support_exclusion);
    const auto & provenance=*original.selection_domain->pure_halo_support_exclusion;
    EXPECT_EQ(provenance.original_rows,11u);
    EXPECT_EQ(provenance.fit_rows,9u);
    EXPECT_EQ(provenance.excluded_rows,2u);
    EXPECT_EQ(provenance.candidate_halos,3u);
    EXPECT_EQ(provenance.retained_halos,2u);
    EXPECT_EQ(provenance.zero_support_halos,1u);
    EXPECT_EQ(original.support[0].size(),2u);
    EXPECT_EQ(original.support[1].size(),1u);
    EXPECT_EQ(original.support[0][1].row,10u);
    EXPECT_EQ(original.support[1][0].row,10u);
    EXPECT_EQ(fitting.support[0].size(),9u);
    EXPECT_DOUBLE_EQ(problem.ObservationScale(),std::sqrt(506.0));
}

TEST(PureHaloSupportExclusionTest, NoHaloKeepsTheWholeSelectedObservationClosure)
{
    auto fixture=BoundaryFixture();
    fixture.model.SelectAllAtoms();
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    ASSERT_TRUE(problem.OriginalInput().selection_domain);
    const auto & provenance=*problem.OriginalInput().selection_domain->pure_halo_support_exclusion;
    EXPECT_EQ(provenance.candidate_halos,0u);
    EXPECT_EQ(provenance.retained_halos,0u);
    EXPECT_EQ(provenance.zero_support_halos,0u);
    EXPECT_EQ(problem.Input().row_ids,problem.OriginalInput().row_ids);
    EXPECT_EQ(problem.FittingToOriginalRows(),(std::vector<std::size_t>{0,1,2,3,4,5,6,7,8,9,10}));
}

TEST(PureHaloSupportExclusionTest, RejectsEmptyFitAndTargetWithoutRetainedSupport)
{
    auto empty=MultipleHaloBoundaryFixture();
    empty.model.FindAtomPtr(10)->SetPosition(0,0,0);
    EXPECT_THROW(BuildPureHaloSupportExclusionProblem(empty.map,empty.model),std::invalid_argument);

    auto target_lost=MultiTargetHaloWorkflowFixture();
    target_lost.model.FindAtomPtr(5)->SetPosition(0,0,0);
    EXPECT_THROW(BuildPureHaloSupportExclusionProblem(target_lost.map,target_lost.model),std::invalid_argument);
}

TEST(PureHaloSupportExclusionTest, ProductionInitializedEntryRejectsUnavailableCalibrationWithoutFallback)
{
    const JointProblem problem(PureInput());
    const auto fit=FitJointComponentsInitialized(problem,TargetInitialization());
    EXPECT_FALSE(fit.assembled_state);
    EXPECT_FALSE(fit.search_completed);
    EXPECT_NE(fit.RuntimeConvergence(),JointCheckStatus::Passed);
    EXPECT_TRUE(fit.available_row_mask==std::vector<bool>(problem.Input().row_ids.size(),false));
    const auto calibration=std::find_if(fit.evidence.begin(),fit.evidence.end(),[](const auto & check) {
        return check.name=="global-width-calibration";
    });
    ASSERT_NE(calibration,fit.evidence.end());
    EXPECT_EQ(calibration->status,JointCheckStatus::Unavailable);
    EXPECT_EQ(calibration->reason,"calibration-unavailable:no-qualified-first-stage-target-widths");
    EXPECT_TRUE(std::none_of(fit.evidence.begin(),fit.evidence.end(),[](const auto & check) {
        return check.name=="active-set-lm-subproblems";
    }));
    const auto saved=core::CaptureJointAnalysisResult(fit);
    const auto decoded=joint_result_io::Decode(joint_result_io::Encode(saved));
    ASSERT_TRUE(decoded.metadata.solver);
    EXPECT_EQ(decoded.metadata.solver->search_method,"gauss-seidel-bcd");
    EXPECT_EQ(decoded.metadata.solver->contract_version,"joint-solver-provenance-v6");
    EXPECT_EQ(decoded.metadata.solver->bcd_policy_version,"gauss-seidel-bcd-production-v1");
    EXPECT_FALSE(IsQualifiedJointAnalysisResult(decoded));
    const auto saved_calibration=std::find_if(decoded.evidence.begin(),decoded.evidence.end(),[](const auto & check) {
        return check.name=="global-width-calibration";
    });
    ASSERT_NE(saved_calibration,decoded.evidence.end());
    EXPECT_EQ(saved_calibration->status,JointCheckStatus::Unavailable);
}

TEST(PureHaloSupportExclusionTest, PeelingUsesOriginalRowsAndReconstructedHalo)
{
    auto fixture=BoundaryFixture();
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    auto snapshot=QualifiedReconstructionSnapshot(problem,fixture.model);
    data_internal::ApplyJointStageEstimates(fixture.model,snapshot,"reconstruction-only-run");
    const SamplingPointList points{{2.5,{2.5,0,0},true}};
    fixture.model.EditAnalysis().SetAtomLocalRawSamplingEntries(*fixture.model.FindAtomPtr(1),
        core::detail::BuildLocalPotentialSampleList(fixture.map,points));
    const auto context=postprocessing_test::MakeJointPostprocessingContext(
        fixture.map,fixture.model,problem,snapshot);
    const auto peelings=postprocessing_test::BuildPostFitPeelingSamples(context, std::vector<std::size_t>{0});
    ASSERT_TRUE(peelings.contains(1));
    ASSERT_FALSE(peelings.at(1).samples.empty());
    EXPECT_TRUE(peelings.at(1).samples.front().response);
    EXPECT_TRUE(peelings.at(1).samples.front().reason.empty());
    const auto prepared=context.BuildPostJointLocalSamples(1,peelings.at(1));
    ASSERT_EQ(prepared.samples.size(),peelings.at(1).samples.size());
    EXPECT_TRUE(prepared.samples.front().fixed_charge_response);
    EXPECT_TRUE(prepared.samples.front().gaussian_response);
}

TEST(PureHaloSupportExclusionTest, GridPeelingAndOwnChargeMatchIndependentDirectReference)
{
    auto fixture=ReconstructionOnlyFixture();
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    ASSERT_EQ(problem.OriginalInput().atom_ids,(std::vector<std::string>{"1","4"}));
    auto snapshot=QualifiedReconstructionSnapshot(problem,fixture.model);
    data_internal::ApplyJointStageEstimates(fixture.model,snapshot,"reconstruction-only-run");
    const std::array<double,3> other_target_position{-4.5,-4.5,-4.5};
    const std::array<std::array<double,3>,2> halo_positions{{
        {4.72,1.18,1.18},{4.72,-1.18,1.18}}};
    const GaussianModel3D point{2.0,.65,-.2};
    const GaussianModel3D other_target_point{10.0,1.1,.4};
    const GaussianModel3D reconstructed_halo{6.0,.875,.1};
    const SamplingPointList points{
        {1.6,{1.6,0.0,0.0},true},
        {1.9,{1.9,0.0,0.0},true}};
    const auto raw=core::detail::BuildLocalPotentialSampleList(fixture.map,points);
    fixture.model.EditAnalysis().SetAtomLocalRawSamplingEntries(*fixture.model.FindAtomPtr(1),raw);

    core::detail::JointPostprocessingMetrics metrics;
    const auto context=postprocessing_test::MakeJointPostprocessingContext(
        fixture.map,fixture.model,problem,snapshot,&metrics);
    const auto peelings=postprocessing_test::BuildPostFitPeelingSamples(context, std::vector<std::size_t>{0});
    const auto & peeling=peelings.at(1);
    ASSERT_FALSE(peeling.samples.empty());
    for (const auto & sample:peeling.samples)
    {
        const auto & observation=sample.observation;
        ASSERT_TRUE(sample.response) << sample.reason;
        double neighbor_prediction=ReferenceContribution(core::simulation::SupportSquare(
            observation.position,other_target_position),other_target_point);
        for(const auto & center:halo_positions)
            neighbor_prediction+=ReferenceContribution(core::simulation::SupportSquare(
                observation.position,center),reconstructed_halo);
        EXPECT_NEAR(*sample.response,observation.raw_response-neighbor_prediction,1e-12);
    }
    EXPECT_EQ(metrics.available_peeling_sample_count,peeling.samples.size());
    const auto local=context.BuildPostJointLocalSamples(1,peeling);
    ASSERT_EQ(local.samples.size(),peeling.samples.size());
    for(std::size_t i=0;i<local.samples.size();++i)
    {
        const auto & sample=local.samples[i];
        ASSERT_TRUE(sample.fixed_charge_response); ASSERT_TRUE(sample.gaussian_response);
        const double charge=ReferenceContribution(sample.squared_distance,point,true);
        EXPECT_NEAR(*sample.fixed_charge_response,charge,1e-12);
        EXPECT_NEAR(*sample.gaussian_response,*peeling.samples[i].response-charge,1e-12);
    }
}

TEST(PureHaloSupportExclusionTest, PeelingIsInvariantToExtraContributorChemicalLabels)
{
    struct Output
    {
        std::vector<PeelingSampleEstimate> peeling;
        std::vector<core::detail::PostJointLocalSample> local_samples;
        std::optional<double> tail_ratio;
    };
    const auto peel_sample=[](MapFixture & fixture) {
        const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
        auto snapshot=QualifiedReconstructionSnapshot(problem,fixture.model);
        data_internal::ApplyJointStageEstimates(fixture.model,snapshot,"reconstruction-only-run");
        const SamplingPointList points{{1.5,{1.5,0.0,0.0},true},{2.8,{2.8,0.0,0.0},true}};
        fixture.model.EditAnalysis().SetAtomLocalRawSamplingEntries(*fixture.model.FindAtomPtr(1),
            core::detail::BuildLocalPotentialSampleList(fixture.map,points));
        const auto context=postprocessing_test::MakeJointPostprocessingContext(
            fixture.map,fixture.model,problem,snapshot);
        auto peeling=postprocessing_test::BuildPostFitPeelingSamples(context, std::vector<std::size_t>{0}).at(1);
        fixture.model.EditAnalysis().ApplyPostFitPeelingsBatch({{1,peeling}});
        const auto tail=AtomLocalPotentialView::For(*fixture.model.FindAtomPtr(1))
            .AssessLocalFittingPeelingRatio(1.0,2.0).ratio;
        const auto local=context.BuildPostJointLocalSamples(1,peeling);
        return Output{std::move(peeling.samples),local.samples,tail};
    };
    auto baseline=ReconstructionOnlyFixture();
    auto relabeled=ReconstructionOnlyFixture();
    for (const int id:{2,3})
    {
        relabeled.model.FindAtomPtr(id)->SetComponentKey(2);
        relabeled.model.FindAtomPtr(id)->SetAtomKey(9);
        relabeled.model.FindAtomPtr(id)->SetElement(Element::OXYGEN);
        relabeled.model.FindAtomPtr(id)->SetComponentID("GLY");
        relabeled.model.FindAtomPtr(id)->SetResidue(Residue::GLY);
        relabeled.model.FindAtomPtr(id)->SetSpot(Spot::O);
    }

    const auto original=peel_sample(baseline);
    const auto changed=peel_sample(relabeled);
    ASSERT_EQ(original.peeling.size(),changed.peeling.size());
    ASSERT_EQ(original.local_samples.size(),changed.local_samples.size());
    for (std::size_t sample=0;sample<original.peeling.size();++sample)
    {
        ASSERT_EQ(original.peeling[sample].response.has_value(),changed.peeling[sample].response.has_value());
        if (original.peeling[sample].response)
            EXPECT_DOUBLE_EQ(*original.peeling[sample].response,*changed.peeling[sample].response);
        EXPECT_EQ(original.peeling[sample].reason,changed.peeling[sample].reason);
        const auto & left=original.local_samples[sample];
        const auto & right=changed.local_samples[sample];
        EXPECT_DOUBLE_EQ(left.raw_response,right.raw_response);
        EXPECT_EQ(left.peeled_response.has_value(),right.peeled_response.has_value());
        if (left.peeled_response) EXPECT_DOUBLE_EQ(*left.peeled_response,*right.peeled_response);
        EXPECT_EQ(left.fixed_charge_response.has_value(),right.fixed_charge_response.has_value());
        if (left.fixed_charge_response)
            EXPECT_DOUBLE_EQ(*left.fixed_charge_response,*right.fixed_charge_response);
        EXPECT_EQ(left.gaussian_response.has_value(),right.gaussian_response.has_value());
        if (left.gaussian_response) EXPECT_DOUBLE_EQ(*left.gaussian_response,*right.gaussian_response);
        EXPECT_EQ(left.reason,right.reason);
    }
    ASSERT_EQ(original.tail_ratio.has_value(),changed.tail_ratio.has_value());
    if (original.tail_ratio) EXPECT_DOUBLE_EQ(*original.tail_ratio,*changed.tail_ratio);
}

TEST(PureHaloSupportExclusionTest, RejectsReconstructionMedianFromDifferentJointEndpoint)
{
    auto fixture=ReconstructionOnlyFixture();
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    auto snapshot=QualifiedReconstructionSnapshot(problem,fixture.model);
    snapshot.halo_median_reconstruction=core::detail::BuildHaloMedianReconstruction(
        fixture.model,problem,snapshot,"reconstruction-only-run").TakeSnapshot();
    snapshot.halo_median_reconstruction->global_median[0]+=1.0;
    data_internal::ApplyJointStageEstimates(fixture.model,snapshot,"reconstruction-only-run");
    const auto context=postprocessing_test::MakeJointPostprocessingContext(
        fixture.map,fixture.model,problem,snapshot);
    EXPECT_THROW(postprocessing_test::BuildPostFitPeelingSamples(context, std::vector<std::size_t>{0}),std::invalid_argument);
}

TEST(PureHaloSupportExclusionTest, HaloSnapshotTransferKeepsDonorWorkspaceAlive)
{
    auto fixture=ReconstructionOnlyFixture();
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    const auto source=QualifiedReconstructionSnapshot(problem,fixture.model);
    auto build=core::detail::BuildHaloMedianReconstruction(
        fixture.model,problem,source,"reconstruction-only-run");
    ASSERT_TRUE(build.snapshot);
    const auto expected=*build.snapshot;
    const auto donor_count=build.target_donors.size();

    const auto transferred=build.TakeSnapshot();

    EXPECT_EQ(transferred,expected);
    EXPECT_FALSE(build.snapshot);
    EXPECT_EQ(build.target_donors.size(),donor_count);
}

TEST(PureHaloSupportExclusionTest, PersistsVersionedDualDomainAndMedianProvenance)
{
    auto fixture=BoundaryFixture();
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    auto snapshot=PurePeelingSnapshot(problem,fixture.model);
    snapshot.initialization=TargetInitialization();
    snapshot.search_completed=true;
    snapshot.objective=0.0;
    snapshot.components.front().search_completed=true;
    snapshot.metadata.solver=JointSolverProvenance{};
    auto & solver=*snapshot.metadata.solver;
    solver.search_method="gauss-seidel-bcd";
    solver.contract_version=std::string(JointSolverProvenanceContractVersion);
    solver.sparse_backend="EIGEN";
    solver.bcd_core_atoms=12;
    solver.bcd_policy_version=std::string(BcdPolicyContractVersion);
    solver.bcd_maximum_sweeps=300;
    solver.bcd_order="forward";
    solver.bcd_local_search=std::string(BcdBoundedLocalSearchName);
    snapshot.evidence.push_back({"global-width-calibration",JointCheckStatus::Passed,
        JointEvidenceScope::AssembledGlobal,.65,2.6,"global-first-stage-p95-k4; qualified-donors=1/1"});
    snapshot.evidence.push_back({"qualified-first-stage-donor-coverage",JointCheckStatus::Passed,
        JointEvidenceScope::AssembledGlobal,1.0,.95,"target-only-qualified-first-stage-widths"});
    const auto encoded=joint_result_io::Encode(snapshot);
    EXPECT_EQ(boost::json::parse(encoded).at("schema_version").as_int64(),24);
    const auto encoded_root=boost::json::parse(encoded).as_object();
    const auto & median_json=encoded_root.at("halo_median_reconstruction").as_object();
    EXPECT_EQ(median_json.at("neighbor_count").as_int64(),4);
    EXPECT_EQ(median_json.at("donor_policy").as_string(),GeometryLocalTargetMedianDonorPolicy);
    const auto & halo_json=median_json.at("halos").as_array().front().as_object();
    EXPECT_FALSE(halo_json.contains("group_key"));
    EXPECT_TRUE(halo_json.contains("donor_atom_ids"));
    const auto decoded=joint_result_io::Decode(encoded);
    EXPECT_EQ(decoded.parameterization_contract,JointPureHaloSupportExclusionParameterizationContract);
    EXPECT_EQ(decoded.fitting_to_original_contributors,problem.FittingToOriginalContributors());
    ASSERT_TRUE(decoded.halo_median_reconstruction);
    EXPECT_EQ(decoded.halo_median_reconstruction->halos,snapshot.halo_median_reconstruction->halos);
}

TEST(PureHaloSupportExclusionTest, RunsFullPureWorkflowWithGeometryLocalHaloSources)
{
    auto fixture=WorkflowFixture(true);
    FitOptions options;
    options.estimator=PotentialEstimator::JOINT_COMPONENTS;
    options.quiet_mode=true;
    options.thread_size=1;
    core::detail::JointPostprocessingMetrics metrics;
    ASSERT_NO_THROW(core::detail::RunJointComponentWorkflow(
        fixture.map,fixture.model,options,&metrics));
    const auto result=fixture.model.GetAnalysisView().GetJointResult();
    ASSERT_TRUE(result);
    EXPECT_EQ(result->parameterization_contract,JointPureHaloSupportExclusionParameterizationContract);
    ASSERT_TRUE(result->selection_domain);
    EXPECT_EQ(result->selection_domain->target_indices.size(),20u);
    EXPECT_TRUE(IsQualifiedJointAnalysisResult(*result));
    const auto calibration=std::find_if(result->evidence.begin(),result->evidence.end(),[](const auto & item) {
        return item.name=="global-width-calibration";
    });
    ASSERT_NE(calibration,result->evidence.end());
    EXPECT_EQ(calibration->status,JointCheckStatus::Passed);
    EXPECT_NE(calibration->reason.find("qualified-donors=20/20"),std::string::npos);
    const auto donor_coverage=std::find_if(result->evidence.begin(),result->evidence.end(),[](const auto & item) {
        return item.name=="qualified-first-stage-donor-coverage";
    });
    ASSERT_NE(donor_coverage,result->evidence.end());
    ASSERT_TRUE(donor_coverage->value && donor_coverage->threshold);
    EXPECT_GE(*donor_coverage->value,*donor_coverage->threshold);
    ASSERT_TRUE(result->halo_median_reconstruction);
    EXPECT_EQ(result->halo_median_reconstruction->halos.size(),2u);
    EXPECT_EQ(metrics.halo_median_valid_target_count,20u);
    EXPECT_EQ(metrics.halo_median_geometry_local_halos,2u);
    EXPECT_EQ(metrics.halo_median_unavailable_halos,0u);
    for (const auto & halo_result:result->halo_median_reconstruction->halos)
    {
        EXPECT_EQ(halo_result.source,HaloMedianSource::GeometryLocalTargetMedian);
        EXPECT_EQ(halo_result.donor_count,4u);
        EXPECT_EQ(halo_result.donor_atom_ids,(std::vector<std::string>{"1","22","14","16"}));
    }
    const auto target=AtomLocalPotentialView::For(*fixture.model.FindAtomPtr(1));
    EXPECT_TRUE(target.GetPostFitPeeling());
    EXPECT_TRUE(target.GetPostJointLocalResult());
    EXPECT_EQ(target.GetStageEstimate(FittingStage::Second).source.role,FittingRole::Target);
    const auto halo=AtomLocalPotentialView::For(*fixture.model.FindAtomPtr(2));
    EXPECT_EQ(halo.GetStageEstimate(FittingStage::Second).source.role,FittingRole::Halo);
    EXPECT_FALSE(halo.GetStageEstimate(FittingStage::Second).point);

    boost::json::array local_evidence;
    for (const auto * atom : fixture.model.GetSelectedAtoms())
    {
        const auto view=AtomLocalPotentialView::For(*atom);
        const auto & local_result=view.GetPostJointLocalResult();
        if (!local_result) continue;
        const auto & r=*local_result;
        boost::json::object row{{"atom",atom->GetSerialID()},{"samples",r.sample_count},
            {"retained",r.effective_sample_count},{"corrected",r.corrected_sample_count},{"positive",r.positive_sample_count},
            {"radius",r.fitting_radius},{"reason",r.reason},{"alpha",r.fit.alpha_r},
            {"ols_A",r.fit.ols.GetModel().GetAmplitude()},{"ols_tau",r.fit.ols.GetModel().GetWidth()},
            {"mdpde_A",r.fit.mdpde.GetModel().GetAmplitude()},{"mdpde_tau",r.fit.mdpde.GetModel().GetWidth()}};
        if (r.information)
        {
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> spectrum(r.information->precision);
            row["information_condition"]=spectrum.eigenvalues().maxCoeff()/spectrum.eigenvalues().minCoeff();
        }
        local_evidence.push_back(std::move(row));
    }
    RecordProperty("local_evidence",boost::json::serialize(local_evidence));
    boost::json::array group_evidence;
    for (const auto key:fixture.model.GetAnalysisView().CollectAtomGroupKeys())
    {
        const auto & summary=fixture.model.GetAnalysisView().GetGroupParameterSummary(key);
        if (!summary) continue;
        boost::json::object row{{"eligible",summary->eligible_count},{"reason",summary->reason}};
        if (summary->inference)
        {
            const auto & g=*summary->inference;
            row["alpha"]=g.alpha_g;
            row["mean_A"]=g.mean.GetAmplitude(); row["mean_tau"]=g.mean.GetWidth();
            row["mdpde_A"]=g.mdpde.GetAmplitude(); row["mdpde_tau"]=g.mdpde.GetWidth();
        }
        group_evidence.push_back(std::move(row));
    }
    RecordProperty("group_evidence",boost::json::serialize(group_evidence));
    const boost::json::object timings{{"context",metrics.context_construction_seconds},
        {"peeling",metrics.peeling_seconds},{"charge",metrics.charge_correction_seconds},
        {"local_alpha",metrics.local_alpha_training_seconds},{"local_fit",metrics.local_mdpde_fitting_seconds},
        {"group_alpha",metrics.group_alpha_training_seconds},{"web",metrics.group_web_seconds},
        {"total",metrics.total_postprocessing_seconds},{"available",metrics.available_peeling_sample_count},
        {"unavailable",metrics.unavailable_peeling_sample_count}};
    RecordProperty("postprocessing",boost::json::serialize(timings));

    EXPECT_GT(metrics.estimated_context_lookup_bytes,0u);
    EXPECT_GT(metrics.estimated_peeling_result_bytes,0u);
    EXPECT_GE(metrics.estimated_peak_postprocessing_payload_bytes,metrics.estimated_peeling_result_bytes);
    RecordProperty("estimated_context_lookup_bytes",
        static_cast<long long>(metrics.estimated_context_lookup_bytes));
    RecordProperty("estimated_reconstruction_workspace_bytes",
        static_cast<long long>(metrics.estimated_reconstruction_workspace_bytes));
    RecordProperty("estimated_peeling_result_bytes",
        static_cast<long long>(metrics.estimated_peeling_result_bytes));
    RecordProperty("estimated_peak_postprocessing_payload_bytes",
        static_cast<long long>(metrics.estimated_peak_postprocessing_payload_bytes));
    if (metrics.rss_before_postprocessing_bytes)
        RecordProperty("rss_before_postprocessing_bytes",
            static_cast<long long>(*metrics.rss_before_postprocessing_bytes));
    if (metrics.rss_after_context_bytes)
        RecordProperty("rss_after_context_bytes",static_cast<long long>(*metrics.rss_after_context_bytes));
    if (metrics.rss_after_peeling_bytes)
        RecordProperty("rss_after_peeling_bytes",static_cast<long long>(*metrics.rss_after_peeling_bytes));
    if (metrics.rss_after_local_bytes)
        RecordProperty("rss_after_local_bytes",static_cast<long long>(*metrics.rss_after_local_bytes));
    if (metrics.rss_after_group_bytes)
        RecordProperty("rss_after_group_bytes",static_cast<long long>(*metrics.rss_after_group_bytes));
}

TEST(PureHaloSupportExclusionTest, RunsPureWorkflowWithoutHalos)
{
    auto fixture=WorkflowFixture();
    fixture.model.SelectAllAtoms();
    FitOptions options;
    options.estimator=PotentialEstimator::JOINT_COMPONENTS;
    options.quiet_mode=true;
    options.thread_size=1;
    core::detail::JointPostprocessingMetrics metrics;
    ASSERT_NO_THROW(core::detail::RunJointComponentWorkflow(
        fixture.map,fixture.model,options,&metrics));
    const auto result=fixture.model.GetAnalysisView().GetJointResult();
    ASSERT_TRUE(result);
    ASSERT_TRUE(result->selection_domain);
    EXPECT_EQ(result->selection_domain->target_indices.size(),result->atom_ids.size());
    EXPECT_FALSE(result->halo_median_reconstruction);
    EXPECT_EQ(metrics.halo_median_geometry_local_halos,0u);
}

TEST(PureHaloSupportExclusionTest, ExtendedSearchQualifiesAssemblyAndFeedsEveryDownstreamStage)
{
    auto fixture=WorkflowFixture(true,true);
    static std::size_t maximum_sweeps;
    maximum_sweeps=0;
    struct ResetHook
    {
        ~ResetHook() {core::detail::JointWorkflowFitHookForTesting()=nullptr;}
    } reset_hook;
    core::detail::JointWorkflowFitHookForTesting()=[](core::JointFitResult & fit) {
        namespace n=core::joint_component;
        ASSERT_TRUE(fit.problem);
        auto initialization=fit.initialization;
        std::fill(initialization.tau.begin(),initialization.tau.end(),1.5);
        const auto calibration=n::CalibrateGlobalWidth(initialization,"extension-workflow-test");
        const auto bounded=n::MakeBoundedInitialization(initialization,calibration);
        ASSERT_TRUE(bounded.available)<<bounded.reason;
        n::BcdExecutionPolicy policy; policy.eta_upper_bound=calibration.eta_max4;
        policy.configure_for_testing=[](auto & local) {
            local.core_atoms=1; local.nominal_sweeps=35; local.maximum_sweeps=100;
        };
        policy.observe_search_for_testing=[](const auto & search) {
            maximum_sweeps=std::max(maximum_sweeps,search.sweep_count);
        };
        auto extended=n::FitBcdComponents(*fit.problem,bounded.tau_values,policy);
        extended.initialization=std::move(initialization);
        for(const auto & check:fit.evidence)
            if(check.name=="global-width-calibration" || check.name=="qualified-first-stage-donor-coverage")
                extended.evidence.push_back(check);
        extended.evidence.push_back({"bounded-initialization-projection",JointCheckStatus::Passed,
            JointEvidenceScope::AssembledGlobal,static_cast<double>(bounded.projection_count),{},"test-start-projection"});
        fit=std::move(extended);
    };
    FitOptions options; options.estimator=PotentialEstimator::JOINT_COMPONENTS;
    options.quiet_mode=true; options.thread_size=1;
    core::detail::JointPostprocessingMetrics metrics;
    ASSERT_NO_THROW(core::detail::RunJointComponentWorkflow(fixture.map,fixture.model,options,&metrics));
    ASSERT_GT(maximum_sweeps,35u);
    ASSERT_LE(maximum_sweeps,100u);
    const auto result=fixture.model.GetAnalysisView().GetJointResult();
    ASSERT_TRUE(result);
    ASSERT_TRUE(IsQualifiedJointAnalysisResult(*result));
    ASSERT_TRUE(result->assembled_state);
    ASSERT_TRUE(result->halo_median_reconstruction);
    EXPECT_EQ(metrics.halo_median_geometry_local_halos,2u);
    for(int serial:{1,23})
    {
        const auto target=AtomLocalPotentialView::For(*fixture.model.FindAtomPtr(serial));
        EXPECT_TRUE(target.GetStageEstimate(FittingStage::Second).point);
        EXPECT_TRUE(target.GetPostFitPeeling());
        EXPECT_TRUE(target.GetPostJointLocalResult());
    }
    EXPECT_GT(metrics.local_mdpde_fitting_seconds,0.0);
    EXPECT_GT(metrics.group_alpha_training_seconds,0.0);
}

TEST(PureHaloSupportExclusionTest, SkipsDownstreamWhenJointQualificationFails)
{
    auto fixture=WorkflowFixture();
    fixture.map.ClearMapValueArray();
    FitOptions options;
    options.estimator=PotentialEstimator::JOINT_COMPONENTS;
    options.quiet_mode=true;
    options.thread_size=1;
    core::detail::JointPostprocessingMetrics metrics;
    ASSERT_NO_THROW(core::detail::RunJointComponentWorkflow(
        fixture.map,fixture.model,options,&metrics));
    const auto result=fixture.model.GetAnalysisView().GetJointResult();
    ASSERT_TRUE(result);
    EXPECT_FALSE(result->initialization.valid);
    ASSERT_FALSE(IsQualifiedJointAnalysisResult(*result));
    EXPECT_FALSE(result->halo_median_reconstruction);
    const auto target=AtomLocalPotentialView::For(*fixture.model.FindAtomPtr(1));
    EXPECT_FALSE(target.GetPostFitPeeling());
    EXPECT_FALSE(target.GetPostJointLocalResult());
}

TEST(PureHaloSupportExclusionTest, ProductionWorkflowUsesBoundedPureExclusionAndStopsWhenCalibrationIsUnavailable)
{
    auto fixture=WorkflowFixture();
    FitOptions options;
    options.estimator=PotentialEstimator::JOINT_COMPONENTS;
    options.quiet_mode=true;
    options.thread_size=1;
    ASSERT_NO_THROW(core::RunPotentialFittingWorkflow(fixture.map,fixture.model,options));
    const auto result=fixture.model.GetAnalysisView().GetJointResult();
    ASSERT_TRUE(result);
    EXPECT_EQ(result->parameterization_contract,JointPureHaloSupportExclusionParameterizationContract);
    ASSERT_TRUE(result->selection_domain);
    EXPECT_EQ(result->selection_domain->contributor_policy,PureHaloSupportExclusionContributorPolicy);
    EXPECT_FALSE(IsQualifiedJointAnalysisResult(*result));
    EXPECT_FALSE(result->halo_median_reconstruction);
    const auto calibration=std::find_if(result->evidence.begin(),result->evidence.end(),[](const auto & check) {
        return check.name=="global-width-calibration";
    });
    ASSERT_NE(calibration,result->evidence.end());
    EXPECT_EQ(calibration->status,JointCheckStatus::Unavailable);
    const auto target=AtomLocalPotentialView::For(*fixture.model.FindAtomPtr(1));
    EXPECT_FALSE(target.GetStageEstimate(FittingStage::Second).point);
    EXPECT_FALSE(target.GetPostFitPeeling());
    EXPECT_FALSE(target.GetPostJointLocalResult());
}

TEST(PureHaloSupportExclusionTest, PublicEstimatorUsesPureExclusion)
{
    auto production_fixture=WorkflowFixture();
    const auto production=EstimateJointComponents(production_fixture.map,production_fixture.model);
    ASSERT_TRUE(production.problem);
    ASSERT_TRUE(production.problem->Input().selection_domain);
    EXPECT_EQ(production.problem->Input().selection_domain->contributor_policy,
        PureHaloSupportExclusionContributorPolicy);
    EXPECT_EQ(production.problem->Input().atom_ids.size(),1u);
    EXPECT_TRUE(std::any_of(production.evidence.begin(),production.evidence.end(),[](const auto & check) {
        return check.name=="global-width-calibration";
    }));

}

TEST(PostJointGridContractTest, SelectsOriginalVoxelsAtInclusiveRadiusWithStableIdentity)
{
    for (const auto spacing : {std::array<double,3>{.5,.5,.5}, std::array<double,3>{.4,.7,1.0}})
    for (const auto center : {std::array<double,3>{0.,0.,0.}, std::array<double,3>{.13,.21,-.17}})
    {
        std::vector<std::unique_ptr<AtomObject>> atoms;
        for (int id=1;id<=2;++id)
        {
            auto atom=std::make_unique<AtomObject>();
            atom->SetSerialID(id); atom->SetElement(Element::CARBON);
            atom->SetPosition(id==1 ? center : std::array<double,3>{2.5,0.,0.});
            atoms.push_back(std::move(atom));
        }
        ModelObject model(std::move(atoms));
        model.SelectAtoms([](const auto & atom){return atom.GetSerialID()==1;});
        auto values=std::make_unique<double[]>(9*9*5);
        for(std::size_t grid=0;grid<9*9*5;++grid) values[grid]=static_cast<double>(grid)+.25;
        MapObject map({9,9,5},spacing,{-2.,-2.,-1.},std::move(values));
        const auto problem=BuildPureHaloSupportExclusionProblem(map,model);
        const auto & input=problem.OriginalInput();
        const auto collect=[&] { return core::detail::CollectPostJointGridObservations(map,input,0,center); };
        const auto samples=collect();
        const auto repeated=collect();
        std::set<std::size_t> expected,actual;
        for(std::size_t grid=0;grid<map.GetMapValueArraySize();++grid)
            if(core::simulation::SupportSquare(map.GetGridPosition(grid),center)<=4.0) expected.insert(grid);
        std::size_t last{}; bool first=true; bool outside_fitting=false;
        for(std::size_t i=0;i<samples.size();++i)
        {
            const auto & sample=samples[i];
            EXPECT_TRUE(first || last<sample.original_row_index);
            first=false; last=sample.original_row_index;
            EXPECT_TRUE(actual.insert(sample.global_grid_index).second);
            EXPECT_EQ(input.row_ids.at(sample.original_row_index),std::to_string(sample.global_grid_index));
            EXPECT_EQ(sample.position,map.GetGridPosition(sample.global_grid_index));
            EXPECT_DOUBLE_EQ(sample.raw_response,map.GetMapValue(sample.global_grid_index));
            EXPECT_DOUBLE_EQ(sample.squared_distance,core::simulation::SupportSquare(sample.position,center));
            EXPECT_EQ(sample.global_grid_index,repeated[i].global_grid_index);
            if(std::find(problem.FittingToOriginalRows().begin(),problem.FittingToOriginalRows().end(),
                    sample.original_row_index)==problem.FittingToOriginalRows().end()) outside_fitting=true;
        }
        EXPECT_EQ(actual,expected); EXPECT_TRUE(outside_fitting);
        auto malformed=input;
        malformed.support[0].push_back(malformed.support[0].front());
        EXPECT_THROW(core::detail::CollectPostJointGridObservations(map,malformed,0,center),std::invalid_argument);
        malformed=input; malformed.row_ids[samples.front().original_row_index]="999999";
        EXPECT_THROW(core::detail::CollectPostJointGridObservations(map,malformed,0,center),std::invalid_argument);
        malformed=input; malformed.observations[samples.front().original_row_index]=std::numeric_limits<double>::quiet_NaN();
        EXPECT_THROW(core::detail::CollectPostJointGridObservations(map,malformed,0,center),std::invalid_argument);
    }
}

TEST(PureHaloSupportExclusionTest, GridSubtractionCoversOverlapsHalosChargesAndSupportBoundaries)
{
    for (bool halos:{false,true}) for (double width:{.01,.4,.9})
    for (const auto center:{std::array<double,3>{0.,0.,0.},std::array<double,3>{.13,.21,-.07}})
    {
        const std::array<std::array<double,3>,3> centers{center,{1.5,0.,0.},{3.,0.,0.}};
        const std::array<GaussianModel3D,3> points{
            GaussianModel3D{2.,width,-.3},GaussianModel3D{3.,.8,.2},
            GaussianModel3D{2.5,(width+.8)/2.,-.05}};
        std::vector<std::unique_ptr<AtomObject>> atoms;
        for(int id=1;id<=(halos ? 3:2);++id)
        {
            auto atom=std::make_unique<AtomObject>();
            atom->SetSerialID(id); atom->SetElement(Element::CARBON);
            atom->SetPosition(centers[static_cast<std::size_t>(id-1)]);
            atoms.push_back(std::move(atom));
        }
        auto hydrogen=std::make_unique<AtomObject>();
        hydrogen->SetSerialID(4); hydrogen->SetElement(Element::HYDROGEN); hydrogen->SetPosition(center);
        atoms.push_back(std::move(hydrogen));
        ModelObject model(std::move(atoms));
        model.SelectAtoms([](const auto & atom){return atom.GetSerialID()<=2;});
        MapObject map({19,13,11},{.5,.4,.6},{-3.,-2.4,-3.});
        auto values=std::make_unique<double[]>(map.GetMapValueArraySize());
        for(std::size_t grid=0;grid<map.GetMapValueArraySize();++grid)
            for(std::size_t atom=0;atom<(halos ? 3u:2u);++atom)
                values[grid]+=ReferenceContribution(core::simulation::SupportSquare(map.GetGridPosition(grid),centers[atom]),points[atom]);
        map.SetMapValueArray(std::move(values));
        const auto problem=BuildPureHaloSupportExclusionProblem(map,model);
        auto snapshot=QualifiedReconstructionSnapshot(problem,model);
        for(std::size_t i=0;i<2;++i)
        {
            snapshot.components[i].state->ac={points[i].GetAmplitude(),points[i].GetOffset()};
            snapshot.components[i].state->tau={points[i].GetWidth()};
        }
        snapshot.halo_median_reconstruction.reset();
        if(halos) snapshot.halo_median_reconstruction=core::detail::BuildHaloMedianReconstruction(
            model,problem,snapshot,"grid-reference").TakeSnapshot();
        data_internal::ApplyJointStageEstimates(model,snapshot,"grid-reference");
        const auto context=postprocessing_test::MakeJointPostprocessingContext(map,model,problem,snapshot);
        const auto batch=context.BuildPostFitPeelingBatch(problem.OriginalInput().selection_domain->target_indices);
        ASSERT_EQ(batch.size(),2u);
        for(const auto & [id,peeling]:batch)
        {
            const auto target=static_cast<std::size_t>(id-1);
            const auto local=context.BuildPostJointLocalSamples(id,peeling);
            ASSERT_EQ(local.samples.size(),peeling.samples.size());
            EXPECT_GT(local.samples.size(),10u);
            for(std::size_t row=0;row<local.samples.size();++row)
            {
                const auto & sample=local.samples[row];
                ASSERT_TRUE(sample.gaussian_response) << sample.reason;
                long double neighbor{}; long double total{};
                for(std::size_t atom=0;atom<(halos ? 3u:2u);++atom)
                {
                    const double value=ReferenceContribution(core::simulation::SupportSquare(sample.position,centers[atom]),points[atom]);
                    total+=value; if(atom!=target) neighbor+=value;
                }
                const double charge=ReferenceContribution(sample.squared_distance,points[target],true);
                const double expected=static_cast<double>(static_cast<long double>(sample.raw_response)-neighbor-charge);
                EXPECT_NEAR(*sample.gaussian_response,expected,1e-10*std::max(1.,std::abs(expected)));
                const double gaussian=ReferenceContribution(sample.squared_distance,points[target])-charge;
                const double alternative=static_cast<double>(sample.raw_response-total+gaussian);
                EXPECT_NEAR(*sample.gaussian_response,alternative,1e-10*std::max(1.,std::abs(expected)));
            }
        }
        if(!halos)
        {
            snapshot.components[1].state.reset();
            data_internal::ApplyJointStageEstimates(model,snapshot,"grid-reference");
            const auto incomplete=postprocessing_test::MakeJointPostprocessingContext(map,model,problem,snapshot)
                .BuildPostFitPeelingBatch(std::vector<std::size_t>{0});
            EXPECT_TRUE(std::any_of(incomplete.front().second.samples.begin(),incomplete.front().second.samples.end(),
                [](const auto & sample){return !sample.response && sample.reason=="missing-contributor-state";}));
        }
    }
}


TEST(PureHaloSupportExclusionTest, GridLocalRangesShareAlphaFitAndInformationDataset)
{
    auto fixture=ReconstructionOnlyFixture();
    auto values=std::make_unique<double[]>(fixture.map.GetMapValueArraySize());
    for(std::size_t grid=0;grid<fixture.map.GetMapValueArraySize();++grid)
        values[grid]=fixture.map.GetMapValue(grid)+.001*std::sin(static_cast<double>(grid));
    fixture.map.SetMapValueArray(std::move(values));
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    const auto snapshot=QualifiedReconstructionSnapshot(problem,fixture.model);
    data_internal::ApplyJointStageEstimates(fixture.model,snapshot,"range-dataset");
    const auto context=postprocessing_test::MakeJointPostprocessingContext(fixture.map,fixture.model,problem,snapshot);
    const auto peeling=context.BuildPostFitPeelingBatch(std::vector<std::size_t>{0}).front().second;
    const auto prepared=context.BuildPostJointLocalSamples(1,peeling);
    const auto samples=prepared.BuildGaussianSamples();
    std::vector<double> responses; for(const auto & sample:samples) responses.push_back(sample.response);
    FitOptions options; options.quiet_mode=true;
    std::size_t previous{};
    for(double radius:{1.,1.5,2.})
    {
        const core::detail::PreparedLocalGaussianDesign design(samples,0.,radius);
        const auto dataset=design.BuildDataset(responses,GaussianModel3D{0.,1.,0.});
        const double alpha=core::detail::TrainLocalGaussianAlpha(dataset,1);
        const auto reference=core::detail::EstimateLocalGaussianDataset(dataset,alpha,1,0.);
        core::detail::JointPostprocessingMetrics metrics;
        const auto result=core::detail::FitPostJointLocalTargetForTesting(context,fixture.model,1,peeling,options,metrics,radius);
        ASSERT_EQ(result.status,EvidenceStatus::Available) << result.reason;
        EXPECT_EQ(result.sample_count,peeling.samples.size());
        EXPECT_EQ(result.corrected_sample_count,samples.size());
        EXPECT_EQ(result.effective_sample_count,static_cast<std::size_t>(dataset.y.size()));
        EXPECT_GT(result.effective_sample_count,previous); previous=result.effective_sample_count;
        EXPECT_DOUBLE_EQ(result.fit.alpha_r,alpha);
        EXPECT_EQ(result.fit.mdpde.GetModel().ToVector(),reference.mdpde.GetModel().ToVector());
        const auto inverse=reference.fit_result->data_covariance.diagonal().cwiseInverse().eval();
        const Eigen::Matrix2d precision=dataset.X.transpose()*inverse.asDiagonal()*dataset.X;
        const Eigen::Vector2d moment=dataset.X.transpose()*(inverse.array()*dataset.y.array()).matrix();
        ASSERT_TRUE(result.information);
        EXPECT_TRUE(result.information->precision.isApprox(precision,1e-12));
        EXPECT_TRUE(result.information->moment.isApprox(moment,1e-12));
    }
}

TEST(PureHaloSupportExclusionTest, GridEvidenceTracksMapRevisionGeometryAndSnapshotLifetime)
{
    for(int mutation=0;mutation<5;++mutation)
    {
        auto fixture=WorkflowFixture(true);
        FitOptions options; options.quiet_mode=true;
        core::detail::RunJointComponentWorkflow(fixture.map,fixture.model,options);
        auto * target=fixture.model.FindAtomPtr(1);
        const auto view=AtomLocalPotentialView::For(*target);
        ASSERT_TRUE(view.GetPostJointLocalResult());
        ASSERT_TRUE(view.GetPostFitPeeling());
        const auto count=view.GetPostFitPeeling()->samples.size();
        auto foreign_metadata=fixture.model.GetAnalysisView().GetJointResult()->metadata;
        foreign_metadata.origin[0]+=1.;
        EXPECT_THROW(fixture.model.EditAnalysis().UpdateJointMetadata(foreign_metadata),std::invalid_argument);
        EXPECT_TRUE(view.GetPostFitPeeling());
        fixture.model.EditAnalysis().SetAtomLocalRawSamplingEntries(*target,{});
        EXPECT_EQ(view.GetPostFitPeeling()->samples.size(),count);
        EXPECT_TRUE(view.GetPostJointLocalResult());
        MapObject independent_map(fixture.map);
        independent_map.ClearMapValueArray();
        EXPECT_TRUE(view.GetPostJointLocalResult());
        ModelObject copy(fixture.model);
        if(mutation==0) fixture.map.ClearMapValueArray();
        else if(mutation==1) fixture.map.MapValueArrayNormalization();
        else if(mutation==2)
        {
            auto values=std::make_unique<double[]>(fixture.map.GetMapValueArraySize());
            fixture.map.SetMapValueArray(std::move(values));
        }
        else if(mutation==3)
        {
            auto position=target->GetPosition(); position[0]+=.1; target->SetPosition(position);
        }
        else fixture.model.EditAnalysis().SetJointResult(*fixture.model.GetAnalysisView().GetJointResult());
        EXPECT_FALSE(fixture.model.GetAnalysisView().GetJointResult().has_value() &&
            AtomLocalPotentialView::For(*target).GetPostFitPeeling().has_value());
        if(mutation<=2)
        {
            EXPECT_FALSE(view.IsAvailable());
            EXPECT_FALSE(copy.GetAnalysisView().GetJointResult());
            EXPECT_FALSE(AtomLocalPotentialView::For(*copy.FindAtomPtr(1)).IsAvailable());
        }
        else EXPECT_FALSE(view.IsAvailable() && view.GetPostJointLocalResult().has_value());
    }
}

TEST(PureHaloSupportExclusionTest, GridPersistenceRejectsForeignIdentitiesAndInvalidPairings)
{
    auto fixture=WorkflowFixture(true);
    for(auto * atom:fixture.model.GetSelectedAtoms())
    {
        atom->SetComponentKey(1); atom->SetAtomKey(7);
    }
    auto values=std::make_unique<double[]>(fixture.map.GetMapValueArraySize());
    for(std::size_t grid=0;grid<fixture.map.GetMapValueArraySize();++grid)
    {
        const auto position=fixture.map.GetGridPosition(grid);
        for(const auto & atom:fixture.model.GetAtomList())
        {
            const int id=atom->GetSerialID();
            const auto point=id==2 ? GaussianModel3D{.4,.8,.1} : id==3 ? GaussianModel3D{.3,.9,-.1} :
                GaussianModel3D{2.+.2*std::sin(static_cast<double>(id)),.65+.05*std::cos(static_cast<double>(id)),-.2};
            values[grid]+=ReferenceContribution(core::simulation::SupportSquare(position,atom->GetPosition()),point);
        }
    }
    fixture.map.SetMapValueArray(std::move(values));
    FitOptions options; options.quiet_mode=true;
    core::detail::RunJointComponentWorkflow(fixture.map,fixture.model,options);
    const auto keys=fixture.model.GetAnalysisView().CollectAtomGroupKeys();
    ASSERT_EQ(keys.size(),1u);
    const auto & group=fixture.model.GetAnalysisView().GetGroupParameterSummary(keys.front());
    ASSERT_TRUE(group && group->inference) << (group ? group->reason:"no group");
    EXPECT_GE(group->eligible_count,18u);
    const auto text=stage_result_io::Encode(fixture.model);
    const auto original=boost::json::parse(text);
    for(int mutation=0;mutation<7;++mutation)
    {
        auto document=original;
        auto & atoms=document.at("atoms").as_array();
        auto found=std::find_if(atoms.begin(),atoms.end(),[](const auto & atom){return !atom.at("peeling").is_null();});
        ASSERT_NE(found,atoms.end());
        auto & sample=found->at("peeling").at("samples").as_array().front();
        if(mutation==0) sample.at("grid")=999999999;
        else if(mutation==1) sample.at("original_row")=999999999;
        else if(mutation==2) sample.at("squared_distance")=4.0001;
        else if(mutation==3) sample.at("position").as_array()[0]=9999.;
        else if(mutation==4) sample.at("reason")="unavailable-with-response";
        else if(mutation==5) found->at("post_joint_local").at("samples")=1;
        else found->at("post_joint_local").at("fixed_charge").as_array()[1]=100.;
        ModelObject rejected(fixture.model);
        EXPECT_THROW(stage_result_io::Decode(rejected,boost::json::serialize(document)),std::invalid_argument);
    }
    auto old=original; old.at("version")=4;
    ModelObject rejected(fixture.model);
    EXPECT_THROW(stage_result_io::Decode(rejected,boost::json::serialize(old)),std::invalid_argument);
    ModelObject restored(fixture.model);
    ASSERT_NO_THROW(stage_result_io::Decode(restored,text));
    const auto round_trip=stage_result_io::Encode(restored);
    auto expected=boost::json::parse(text),actual=boost::json::parse(round_trip);
    const auto sort_groups=[](auto & document) {
        auto & groups=document.at("groups").as_array();
        std::sort(groups.begin(),groups.end(),[](const auto & left,const auto & right) {
            return boost::json::value_to<std::uint64_t>(left.at("key"))<boost::json::value_to<std::uint64_t>(right.at("key"));
        });
    };
    sort_groups(expected); sort_groups(actual);
    EXPECT_TRUE(expected==actual);
    const command_test::ScopedTempDir directory{"grid_group_roundtrip"};
    DataRepository repository(directory.path()/"model.sqlite");
    ASSERT_NO_THROW(repository.SaveModel(fixture.model,"grid-group"));
    const auto loaded=repository.LoadModel("grid-group");
    const auto & loaded_group=loaded->GetAnalysisView().GetGroupParameterSummary(keys.front());
    ASSERT_TRUE(loaded_group && loaded_group->inference);
    EXPECT_EQ(loaded_group->member_ids,group->member_ids);
    for(std::size_t i=0;i<group->inference->member_results.size();++i)
    {
        const auto & reference=group->inference->member_results[i];
        const auto & saved=loaded_group->inference->member_results[i];
        EXPECT_EQ(saved.posterior.GetModel().ToVector(),reference.posterior.GetModel().ToVector());
        EXPECT_EQ(*saved.parameter_covariance,*reference.parameter_covariance);
        ASSERT_TRUE(saved.fixed_charge); ASSERT_TRUE(reference.fixed_charge);
        EXPECT_DOUBLE_EQ(saved.fixed_charge->coefficient,reference.fixed_charge->coefficient);
        EXPECT_DOUBLE_EQ(saved.fixed_charge->width,reference.fixed_charge->width);
    }
}

TEST(PureHaloSupportExclusionTest, NonpositiveGridCorrectionsAreNotFailuresOrPhantomSolverRows)
{
    auto fixture=BoundaryFixture();
    auto values=std::make_unique<double[]>(fixture.map.GetMapValueArraySize());
    std::fill_n(values.get(),fixture.map.GetMapValueArraySize(),-1.);
    fixture.map.SetMapValueArray(std::move(values));
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    const auto snapshot=QualifiedReconstructionSnapshot(problem,fixture.model);
    data_internal::ApplyJointStageEstimates(fixture.model,snapshot,"reconstruction-only-run");
    const auto context=postprocessing_test::MakeJointPostprocessingContext(fixture.map,fixture.model,problem,snapshot);
    const auto peeling=context.BuildPostFitPeelingBatch(std::vector<std::size_t>{0}).front().second;
    FitOptions options; options.quiet_mode=true;
    core::detail::JointPostprocessingMetrics metrics;
    const auto result=core::detail::FitPostJointLocalTargetForTesting(context,fixture.model,1,peeling,options,metrics);
    EXPECT_EQ(result.status,EvidenceStatus::Unavailable);
    EXPECT_EQ(result.reason,"insufficient-positive-responses");
    EXPECT_EQ(result.corrected_sample_count,peeling.samples.size());
    EXPECT_EQ(result.positive_sample_count,0u);
    EXPECT_EQ(result.effective_sample_count,0u);
    EXPECT_TRUE(result.sample_failures.empty());
    EXPECT_FALSE(result.information);
    fixture.model.EditAnalysis().SetAtomPostFitPeeling(*fixture.model.FindAtomPtr(1),peeling);
    EXPECT_NO_THROW(fixture.model.EditAnalysis().SetAtomPostJointLocalResult(*fixture.model.FindAtomPtr(1),result));
}
