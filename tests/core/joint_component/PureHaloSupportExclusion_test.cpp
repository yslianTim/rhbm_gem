#include <gtest/gtest.h>
#include <boost/json.hpp>

#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/MapSampler.hpp"
#include "core/detail/joint_component/HaloMedianReconstruction.hpp"
#include "core/detail/GroupPotentialFitting.hpp"
#include "core/detail/joint_component/PostFitPeeling.hpp"
#include "core/detail/PotentialFittingWorkflow.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/command/detail/SimulationGeometry.hpp"
#include "data/detail/AtomClassifier.hpp"
#include "data/detail/JointStageAdapter.hpp"
#include "data/io/detail/JointResultJson.hpp"
#include "data/io/detail/StageResultJson.hpp"
#include <rhbm_gem/core/JointComponentEstimator.hpp>
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
    initialization.b={.65};
    JointInitializationAtom atom;
    atom.id="1";
    atom.target=true;
    atom.reason="valid-width";
    atom.used_b=.65;
    atom.original_b=.65;
    initialization.atoms.push_back(std::move(atom));
    initialization.reason="valid-widths";
    return initialization;
}

void FitWithLowLevelUnboundedOracle(JointFitResult & fit)
{
    if (!fit.problem) throw std::runtime_error("Joint test hook lacks a problem snapshot.");
    const auto problem=*fit.problem;
    auto initialization=fit.initialization;
    fit=core::FitJointComponents(problem,initialization.b);
    fit.initialization=std::move(initialization);
}

struct ScopedLowLevelUnboundedOracle
{
    ScopedLowLevelUnboundedOracle()
    {
        core::detail::JointWorkflowFitHookForTesting()=&FitWithLowLevelUnboundedOracle;
    }
    ~ScopedLowLevelUnboundedOracle()
    {
        core::detail::JointWorkflowFitHookForTesting()=nullptr;
    }
};

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

MapFixture WorkflowFixture()
{
    std::vector<std::unique_ptr<AtomObject>> atoms;
    const auto add_atom=[&](int serial,double x,int component,int atom_key) {
        auto atom=std::make_unique<AtomObject>();
        atom->SetSerialID(serial);
        atom->SetElement(Element::CARBON);
        atom->SetComponentKey(static_cast<ComponentKey>(component));
        atom->SetAtomKey(static_cast<AtomKey>(atom_key));
        atom->SetChainID("A");
        atom->SetComponentID("ALA");
        atom->SetAtomID(serial==1 ? "CA" : serial==2 ? "CB" : "CG");
        atom->SetPosition(x,0,0);
        atoms.push_back(std::move(atom));
    };
    add_atom(1,0.0,1,7);
    add_atom(2,4.5,1,7);
    add_atom(3,4.8,3,8);
    ModelObject model(std::move(atoms));
    model.SelectAtoms([](const auto & atom){return atom.GetSerialID()==1;});

    constexpr std::array<int,3> dimensions{45,45,45};
    constexpr std::array<double,3> spacing{.3,.3,.3};
    constexpr std::array<double,3> origin{-6.6,-6.6,-6.6};
    const std::array<std::array<double,3>,3> positions{{{0,0,0},{4.5,0,0},{4.8,0,0}}};
    const std::array<double,3> amplitudes{{2.0,.4,.3}};
    const std::array<double,3> widths{{.65,.8,.9}};
    const std::array<double,3> charges{{-.2,.1,-.1}};
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
    result.components.front().state->b={.65};
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
        point.parameter_stability="NotAssessed";
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


double ReferenceCubic(double p0, double p1, double p2, double p3, double t)
{
    const double a0 = p1;
    const double a1 = .5 * (p2 - p0);
    const double a2 = .5 * (-p3 + 4.0 * p2 - 5.0 * p1 + 2.0 * p0);
    const double a3 = .5 * (p3 - 3.0 * p2 + 3.0 * p1 - p0);
    return a3 * t * t * t + a2 * t * t + a1 * t + a0;
}

double ReferenceTricubicCharge(
    const MapObject & map, const std::vector<double> & original_charge,
    const std::array<double, 3> & position)
{
    const auto dimensions = map.GetGridSize();
    const auto spacing = map.GetGridSpacing();
    const auto origin = map.GetOrigin();
    std::array<int, 3> base{};
    std::array<double, 3> local{};
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
        const double coordinate = (position[axis] - origin[axis]) / spacing[axis];
        base[axis] = static_cast<int>(std::floor(coordinate));
        local[axis] = coordinate - static_cast<double>(base[axis]);
    }
    const auto read = [&](int x, int y, int z) {
        x = std::clamp(x, 0, dimensions[0] - 1);
        y = std::clamp(y, 0, dimensions[1] - 1);
        z = std::clamp(z, 0, dimensions[2] - 1);
        const auto index = static_cast<std::size_t>(x) + static_cast<std::size_t>(dimensions[0]) *
            (static_cast<std::size_t>(y) + static_cast<std::size_t>(dimensions[1]) *
                static_cast<std::size_t>(z));
        return original_charge[index];
    };
    std::array<std::array<std::array<double, 4>, 4>, 4> values{};
    for (int i = -1; i < 3; ++i) for (int j = -1; j < 3; ++j) for (int k = -1; k < 3; ++k)
        values[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(j + 1)][static_cast<std::size_t>(k + 1)] =
            read(base[0] + i, base[1] + j, base[2] + k);
    std::array<std::array<double, 4>, 4> along_x{};
    for (std::size_t j = 0; j < 4; ++j) for (std::size_t k = 0; k < 4; ++k)
        along_x[j][k] = ReferenceCubic(
            values[0][j][k], values[1][j][k], values[2][j][k], values[3][j][k], local[0]);
    std::array<double, 4> along_y{};
    for (std::size_t k = 0; k < 4; ++k)
        along_y[k] = ReferenceCubic(
            along_x[0][k], along_x[1][k], along_x[2][k], along_x[3][k], local[1]);
    return ReferenceCubic(along_y[0], along_y[1], along_y[2], along_y[3], local[2]);
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
    auto snapshot=PurePeelingSnapshot(problem,fixture.model);
    data_internal::ApplyJointStageEstimates(fixture.model,snapshot,"pure-peeling-run");
    const SamplingPointList points{{2.5,{2.5,0,0},true}};
    fixture.model.EditAnalysis().SetAtomLocalRawSamplingEntries(*fixture.model.FindAtomPtr(1),
        core::detail::BuildLocalPotentialSampleList(fixture.map,points));
    const core::detail::JointPostprocessingContext context(
        fixture.map,fixture.model,problem,&snapshot);
    const auto peelings=context.BuildPostFitPeelingSamples(std::vector<std::size_t>{0});
    ASSERT_TRUE(peelings.contains(1));
    ASSERT_EQ(peelings.at(1).samples.size(),1u);
    EXPECT_TRUE(peelings.at(1).samples.front().response);
    EXPECT_TRUE(peelings.at(1).samples.front().reason.empty());
    const auto prepared=context.BuildPostJointLocalSamples(1,peelings.at(1));
    ASSERT_EQ(prepared.samples.size(),1u);
    EXPECT_TRUE(prepared.samples.front().fixed_charge_response);
    EXPECT_TRUE(prepared.samples.front().gaussian_response);
}

TEST(PureHaloSupportExclusionTest, ReconstructsMissingNodesFromAllQualifiedHaloContributors)
{
    auto fixture=ReconstructionOnlyFixture();
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    ASSERT_EQ(problem.OriginalInput().atom_ids,(std::vector<std::string>{"1","4"}));
    auto snapshot=QualifiedReconstructionSnapshot(problem,fixture.model);
    data_internal::ApplyJointStageEstimates(fixture.model,snapshot,"reconstruction-only-run");
    const std::array<double,3> target_position{0.0,0.0,0.0};
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
    const core::detail::JointPostprocessingContext context(
        fixture.map,fixture.model,problem,&snapshot,&metrics);
    const auto peelings=context.BuildPostFitPeelingSamples(std::vector<std::size_t>{0});
    const auto & peeling=peelings.at(1);
    ASSERT_EQ(peeling.samples.size(),2u);

    const auto map_size=fixture.map.GetMapValueArraySize();
    std::vector<double> target_field(map_size),other_target_field(map_size),halo_field(map_size),charge_field(map_size);
    for (std::size_t grid=0;grid<map_size;++grid)
    {
        const auto position=fixture.map.GetGridPosition(grid);
        const auto contribution=[&](const auto & center,const GaussianModel3D & profile) {
            const double square=core::simulation::SupportSquare(position,center);
            if (square>6.25) return 0.0;
            const auto basis=core::joint_component::EvaluateKernel(square,profile.GetWidth(),2.5);
            return profile.GetAmplitude()*basis.gaussian+profile.GetOffset()*basis.charge;
        };
        target_field[grid]=contribution(target_position,point);
        other_target_field[grid]=contribution(other_target_position,other_target_point);
        for (const auto & halo:halo_positions) halo_field[grid]+=contribution(halo,reconstructed_halo);
        const double target_square=core::simulation::SupportSquare(position,target_position);
        if (target_square<=6.25)
            charge_field[grid]=point.GetOffset()*core::joint_component::EvaluateKernel(
                target_square,point.GetWidth(),2.5).charge;
    }
    for (std::size_t sample=0;sample<points.size();++sample)
    {
        ASSERT_TRUE(peeling.samples[sample].response) << peeling.samples[sample].reason;
        const double neighbor_prediction=ReferenceTricubicCharge(
            fixture.map,halo_field,points[sample].position)+ReferenceTricubicCharge(
                fixture.map,other_target_field,points[sample].position);
        const double expected=raw[sample].response-neighbor_prediction;
        EXPECT_NEAR(*peeling.samples[sample].response,expected,1e-12);
    }
    EXPECT_NEAR(*peeling.samples[0].response,raw[0].response,1e-12);
    EXPECT_EQ(metrics.available_peeling_sample_count,2u);
    EXPECT_GE(metrics.original_support_boundary_samples,1u);
    EXPECT_EQ(metrics.reconstruction_lookup_atom_count,4u);
    EXPECT_GT(metrics.reconstruction_only_node_evaluation_count,0u);
    EXPECT_LE(metrics.reconstruction_node_cache_peak_entries,4096u);
    EXPECT_GT(metrics.estimated_reconstruction_workspace_bytes,0u);

    const auto local=context.BuildPostJointLocalSamples(1,peeling);
    ASSERT_EQ(local.samples.size(),2u);
    ASSERT_TRUE(local.samples[1].fixed_charge_response);
    ASSERT_TRUE(local.samples[1].gaussian_response);
    const double expected_charge=ReferenceTricubicCharge(
        fixture.map,charge_field,points[1].position);
    EXPECT_NEAR(*local.samples[1].fixed_charge_response,expected_charge,1e-12);
    EXPECT_NEAR(*local.samples[1].gaussian_response,
        *peeling.samples[1].response-expected_charge,1e-12);
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
        const core::detail::JointPostprocessingContext context(
            fixture.map,fixture.model,problem,&snapshot);
        auto peeling=context.BuildPostFitPeelingSamples(std::vector<std::size_t>{0}).at(1);
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
    const core::detail::JointPostprocessingContext context(
        fixture.map,fixture.model,problem,&snapshot);
    EXPECT_THROW(context.BuildPostFitPeelingSamples(std::vector<std::size_t>{0}),std::invalid_argument);
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
    solver.search_method="fixed-neighbor";
    solver.contract_version=std::string(JointSolverProvenanceContractVersion);
    solver.sparse_backend="EIGEN";
    solver.fixed_neighbor_core_atoms=12;
    solver.fixed_neighbor_policy_version=std::string(FixedNeighborPolicyContractVersion);
    solver.fixed_neighbor_maximum_sweeps=200;
    solver.fixed_neighbor_order="forward";
    solver.fixed_neighbor_local_search=std::string(FixedNeighborBoundedLocalSearchName);
    snapshot.evidence.push_back({"global-width-calibration",JointCheckStatus::Passed,
        JointEvidenceScope::AssembledGlobal,.65,2.6,"global-first-stage-p95-k4; qualified-donors=1/1"});
    snapshot.evidence.push_back({"qualified-first-stage-donor-coverage",JointCheckStatus::Passed,
        JointEvidenceScope::AssembledGlobal,1.0,.95,"target-only-qualified-first-stage-widths"});
    const auto encoded=joint_result_io::Encode(snapshot);
    EXPECT_EQ(boost::json::parse(encoded).at("schema_version").as_int64(),19);
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
    auto fixture=WorkflowFixture();
    FitOptions options;
    options.estimator=PotentialEstimator::JOINT_COMPONENTS;
    options.quiet_mode=true;
    options.thread_size=1;
    core::detail::JointPostprocessingMetrics metrics;
    ScopedLowLevelUnboundedOracle unbounded_test_solver;
    ASSERT_NO_THROW(core::detail::RunJointComponentWorkflow(
        fixture.map,fixture.model,options,&metrics));
    const auto result=fixture.model.GetAnalysisView().GetJointResult();
    ASSERT_TRUE(result);
    EXPECT_EQ(result->parameterization_contract,JointPureHaloSupportExclusionParameterizationContract);
    ASSERT_TRUE(result->halo_median_reconstruction);
    EXPECT_EQ(result->halo_median_reconstruction->halos.size(),2u);
    EXPECT_EQ(metrics.halo_median_valid_target_count,1u);
    EXPECT_EQ(metrics.halo_median_geometry_local_halos,2u);
    EXPECT_EQ(metrics.halo_median_unavailable_halos,0u);
    for (const auto & halo_result:result->halo_median_reconstruction->halos)
    {
        EXPECT_EQ(halo_result.source,HaloMedianSource::GeometryLocalTargetMedian);
        EXPECT_EQ(halo_result.donor_count,1u);
        EXPECT_EQ(halo_result.donor_atom_ids,(std::vector<std::string>{"1"}));
    }
    const auto target=AtomLocalPotentialView::For(*fixture.model.FindAtomPtr(1));
    EXPECT_TRUE(target.GetPostFitPeeling());
    EXPECT_TRUE(target.GetPostJointLocalResult());
    EXPECT_EQ(target.GetStageEstimate(FittingStage::Second).source.role,FittingRole::Target);
    const auto halo=AtomLocalPotentialView::For(*fixture.model.FindAtomPtr(2));
    EXPECT_EQ(halo.GetStageEstimate(FittingStage::Second).source.role,FittingRole::Halo);
    EXPECT_FALSE(halo.GetStageEstimate(FittingStage::Second).point);

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
