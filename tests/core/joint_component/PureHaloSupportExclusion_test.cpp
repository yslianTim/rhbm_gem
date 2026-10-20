#include <gtest/gtest.h>
#include <boost/json.hpp>

#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/MapSampler.hpp"
#include "core/detail/HaloMedianReconstruction.hpp"
#include "core/detail/GroupPotentialFitting.hpp"
#include "core/detail/PostFitPeeling.hpp"
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

void FitWithInternalUnboundedOracle(JointFitResult & fit)
{
    if (!fit.problem) throw std::runtime_error("Joint test hook lacks a problem snapshot.");
    fit=core::detail::FitJointComponentsInitializedWithProgress(
        *fit.problem,fit.initialization,{},false);
}

struct ScopedInternalUnboundedOracle
{
    ScopedInternalUnboundedOracle()
    {
        core::detail::JointWorkflowFitHookForTesting()=&FitWithInternalUnboundedOracle;
    }
    ~ScopedInternalUnboundedOracle()
    {
        core::detail::JointWorkflowFitHookForTesting()=nullptr;
    }
};

JointAnalysisResult PurePeelingSnapshot(const JointProblem & problem,const ModelObject & model)
{
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
    component.state=JointState{{1.0,0.0},{1.0},{0.0},{0.0},0.0,{}};
    result.components.push_back(component);
    result.assembled_state=component.state;
    JointHaloMedianSnapshot medians;
    medians.joint_run_id="pure-peeling-run";
    medians.valid_target_count=1;
    medians.global_median={.5,.8,0.0};
    const auto * halo=model.FindAtomPtr(2);
    if (problem.OriginalInput().atom_ids.size()>
        problem.OriginalInput().selection_domain->target_indices.size())
    {
        medians.halos.push_back({1,"2",data_internal::GetGroupKey(halo),.5,.8,0.0,
            HaloMedianSource::GlobalTargetMedian,1,"pure-peeling-run",
            std::string(HaloMedianReconstructionContractVersion)});
        result.halo_median_reconstruction=std::move(medians);
    }
    return result;
}

JointAnalysisResult QualifiedPurePeelingSnapshot(const JointProblem & problem,const ModelObject & model)
{
    auto result=PurePeelingSnapshot(problem,model);
    result.initialization=TargetInitialization();
    result.search_completed=true;
    result.components.front().search_completed=true;
    result.components.front().stop_reason="converged";
    result.objective=0.0;
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

MapFixture WideBoundaryFixture()
{
    auto fixture=BoundaryFixture();
    constexpr std::array<int,3> dimensions{21,9,9};
    auto values=std::make_unique<double[]>(static_cast<std::size_t>(dimensions[0])*dimensions[1]*dimensions[2]);
    for (std::size_t i=0;i<static_cast<std::size_t>(dimensions[0])*dimensions[1]*dimensions[2];++i)
        values[i]=static_cast<double>(i);
    MapObject map(dimensions,{.5,.5,.5},{-5,-2,-2},std::move(values));
    return {std::move(fixture.model),std::move(map)};
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
        component.state=JointState{{10.0,.4},{1.1},{0.0},{0.0},0.0,{}};
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
    const auto medians=core::detail::BuildHaloMedianReconstruction(
        model,problem,result,"reconstruction-only-run");
    if (problem.OriginalInput().atom_ids.size()>targets.size())
        result.halo_median_reconstruction=medians.snapshot;
    return result;
}

struct IndexedChargeCorrectionFixture
{
    ModelObject model;
    MapObject map;
    JointProblem problem;
};

IndexedChargeCorrectionFixture MakeIndexedChargeCorrectionFixture()
{
    std::vector<std::unique_ptr<AtomObject>> atoms;
    for (int serial = 1; serial <= 3; ++serial)
    {
        auto atom = std::make_unique<AtomObject>();
        atom->SetSerialID(serial);
        atom->SetElement(Element::CARBON);
        atom->SetPosition(static_cast<double>(serial + 3), 0, 0);
        atoms.push_back(std::move(atom));
    }
    ModelObject model(std::move(atoms));
    model.SelectAtoms([](const auto & atom) { return atom.GetSerialID() != 1; });

    auto values = std::make_unique<double[]>(10);
    std::fill_n(values.get(), 10, 0.0);
    MapObject map({10, 1, 1}, {1, 1, 1}, {0, 0, 0}, std::move(values));

    JointProblemInput original;
    original.atom_ids = {"1", "2", "3"};
    original.row_ids = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"};
    original.observations.assign(original.row_ids.size(), 0.0);
    original.support = {
        {{1, .3}, {4, 1.1}, {7, 2.2}},
        {{0, .2}, {2, 1.3}, {5, 2.1}, {8, .7}},
        {{3, .4}, {6, 1.7}, {9, 2.8}}
    };
    original.selection_domain.emplace();
    original.selection_domain->contributor_policy = std::string(PureHaloSupportExclusionContributorPolicy);
    original.selection_domain->target_indices = {1, 2};
    original.selection_domain->pure_halo_support_exclusion =
        JointPureHaloSupportExclusionProvenance{10, 7, 3, 1, 1, 0};

    JointProblemInput fitting;
    fitting.atom_ids = {"2", "3"};
    fitting.row_ids = {"0", "2", "3", "5", "6", "8", "9"};
    fitting.observations.assign(fitting.row_ids.size(), 0.0);
    fitting.support = {
        {{0, .2}, {1, 1.3}, {3, 2.1}, {5, .7}},
        {{2, .4}, {4, 1.7}, {6, 2.8}}
    };
    fitting.selection_domain.emplace();
    fitting.selection_domain->contributor_policy = std::string(PureHaloSupportExclusionContributorPolicy);
    fitting.selection_domain->target_indices = {0, 1};
    fitting.selection_domain->pure_halo_support_exclusion =
        JointPureHaloSupportExclusionProvenance{7, 7, 0, 0, 0, 0};

    const auto problem = JointProblem(std::move(fitting)).WithReconstructionDomain(
        std::move(original), {1, 2}, {0, 2, 3, 5, 6, 8, 9});
    return {std::move(model), std::move(map), std::move(problem)};
}

JointAnalysisResult IndexedChargeCorrectionSnapshot(
    const JointProblem & problem, const ModelObject & model, double first_charge,
    double second_charge, double halo_charge = .2)
{
    const auto & original = problem.OriginalInput();
    JointAnalysisResult result;
    result.parameterization_contract = std::string(JointPureHaloSupportExclusionParameterizationContract);
    result.atom_ids = original.atom_ids;
    result.row_ids = original.row_ids;
    result.fitting_row_ids = problem.Input().row_ids;
    result.fitting_to_original_contributors = problem.FittingToOriginalContributors();
    result.fitting_to_original_rows = problem.FittingToOriginalRows();
    result.selection_domain = original.selection_domain;
    result.available_row_mask.assign(result.row_ids.size(), false);
    for (const auto row : result.fitting_to_original_rows) result.available_row_mask[row] = true;
    result.runtime_convergence = JointCheckStatus::Passed;
    result.target_runtime_convergence = JointCheckStatus::Passed;

    JointAnalysisComponent component;
    component.id = "indexed-target-component";
    component.atoms = {1, 2};
    component.rows = result.fitting_to_original_rows;
    component.runtime_convergence = JointCheckStatus::Passed;
    component.target_runtime_convergence = JointCheckStatus::Passed;
    component.state = JointState{
        {2.0, first_charge, 1.4, second_charge},
        {.65, 1.15},
        {std::log(.65), std::log(1.15)},
        {0.0, 0.0},
        0.0,
        {}};
    result.components.push_back(component);
    result.assembled_state = component.state;

    JointHaloMedianSnapshot median;
    median.joint_run_id = "indexed-charge-run";
    median.valid_target_count = 2;
    median.global_median = {.5, .8, halo_charge};
    const auto * halo = model.FindAtomPtr(1);
    median.halos.push_back({0, "1", data_internal::GetGroupKey(halo), .5, .8, halo_charge,
        HaloMedianSource::GlobalTargetMedian, 2, "indexed-charge-run",
        std::string(HaloMedianReconstructionContractVersion)});
    result.halo_median_reconstruction = std::move(median);
    return result;
}

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
    EXPECT_FALSE(original.selection_domain->fixed_halo_selection);
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

TEST(PureHaloSupportExclusionTest, RetiredFixedHaloSelectionCannotCreateRuntimeProblem)
{
    auto input=PureInput();
    input.atom_ids.push_back("2");
    input.support.push_back({{0,.25}});
    input.selection_domain->contributor_policy=std::string(FixedHaloContributorPolicy);
    input.selection_domain->pure_halo_support_exclusion.reset();
    input.selection_domain->fixed_halo_selection=JointFixedHaloSelectionProvenance{1,1,1,0};
    EXPECT_THROW(JointProblem(std::move(input)),std::invalid_argument);
}

TEST(PureHaloSupportExclusionTest, InternalUnboundedOracleFitsTargetOnly)
{
    const JointProblem problem(PureInput());
    const auto fit=core::detail::FitJointComponentsInitializedWithProgress(
        problem,TargetInitialization(),{},false);
    ASSERT_TRUE(fit.assembled_state);
    EXPECT_EQ(fit.RuntimeConvergence(),JointCheckStatus::Passed);
    EXPECT_EQ(fit.TargetRuntimeConvergence(),JointCheckStatus::Passed);
    EXPECT_TRUE(fit.available_row_mask==std::vector<bool>(problem.Input().row_ids.size(),true));
    EXPECT_TRUE(fit.evidence.end()!=std::find_if(fit.evidence.begin(),fit.evidence.end(),[](const auto & check) {
        return check.name=="independent-replay" && check.status==JointCheckStatus::Passed;
    }));
    EXPECT_FALSE(std::any_of(fit.evidence.begin(),fit.evidence.end(),[](const auto & check) {
        return check.name=="fixed-halo-replay";
    }));
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

TEST(PureHaloSupportExclusionTest, CapturesOriginalAndFittingIdentitiesSeparately)
{
    auto fixture=BoundaryFixture();
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    const auto snapshot=CaptureJointAnalysisResult(FitJointComponents(problem,{.65}));
    EXPECT_EQ(snapshot.parameterization_contract,JointPureHaloSupportExclusionParameterizationContract);
    EXPECT_EQ(snapshot.atom_ids,problem.OriginalInput().atom_ids);
    EXPECT_EQ(snapshot.row_ids,problem.OriginalInput().row_ids);
    EXPECT_EQ(snapshot.fitting_row_ids,problem.Input().row_ids);
    EXPECT_EQ(snapshot.fitting_to_original_contributors,problem.FittingToOriginalContributors());
    EXPECT_EQ(snapshot.fitting_to_original_rows,problem.FittingToOriginalRows());
    ASSERT_EQ(snapshot.available_row_mask.size(),snapshot.row_ids.size());
    EXPECT_FALSE(snapshot.available_row_mask.back());
    for (const auto & component:snapshot.components)
    {
        for (const auto atom:component.atoms) EXPECT_LT(atom,snapshot.atom_ids.size());
        for (const auto row:component.rows) EXPECT_LT(row,snapshot.row_ids.size());
    }
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

TEST(PureHaloSupportExclusionTest, ReproducesPostJointLocalChargeCorrectionDomainMismatch)
{
    auto fixture = MakeIndexedChargeCorrectionFixture();
    const auto & original = fixture.problem.OriginalInput();
    const auto & fitting = fixture.problem.Input();
    ASSERT_EQ(fixture.problem.FittingToOriginalContributors(), (std::vector<std::size_t>{1, 2}));
    ASSERT_EQ(fixture.problem.FittingToOriginalRows(),
        (std::vector<std::size_t>{0, 2, 3, 5, 6, 8, 9}));
    ASSERT_NE(fixture.problem.FittingToOriginalContributors()[0], 0u);
    ASSERT_NE(fixture.problem.FittingToOriginalRows()[1], 1u);
    ASSERT_EQ(original.support[1][2].row, 5u);
    ASSERT_EQ(fitting.support[0][2].row, 3u);

    const auto result = IndexedChargeCorrectionSnapshot(fixture.problem, fixture.model, .35, -.45);
    data_internal::ApplyJointStageEstimates(fixture.model, result, "indexed-charge-run");
    fixture.model.EditAnalysis().SetAtomLocalRawSamplingEntries(*fixture.model.FindAtomPtr(2),
        core::detail::BuildLocalPotentialSampleList(fixture.map,
            SamplingPointList{{0.0, {5.25, 0.0, 0.0}, true}}));
    fixture.model.EditAnalysis().SetAtomLocalRawSamplingEntries(*fixture.model.FindAtomPtr(3),
        core::detail::BuildLocalPotentialSampleList(fixture.map,
            SamplingPointList{{0.0, {6.25, 0.0, 0.0}, true}}));

    core::detail::JointPostprocessingContext context(
        fixture.map, fixture.model, fixture.problem, &result);
    const auto peelings = context.BuildPostFitPeelingSamples(std::vector<std::size_t>{1, 2});
    ASSERT_TRUE(peelings.contains(2));
    ASSERT_TRUE(peelings.contains(3));

    const auto target_one = context.BuildPostJointLocalSamples(2, peelings.at(2));
    ASSERT_EQ(target_one.samples.size(), 1u);
    ASSERT_TRUE(target_one.samples.front().peeled_response);
    ASSERT_TRUE(target_one.samples.front().fixed_charge_response);
    ASSERT_TRUE(target_one.samples.front().gaussian_response);
    std::vector<double> target_one_charge(original.row_ids.size(), 0.0);
    for (const auto & support : original.support[1])
        target_one_charge[support.row] = .35 * core::joint_component::EvaluateKernel(
            support.squared_distance, .65, 2.5).charge;
    const double expected_target_one_charge = ReferenceTricubicCharge(
        fixture.map, target_one_charge, {5.25, 0.0, 0.0});
    EXPECT_NEAR(*target_one.samples.front().fixed_charge_response,
        expected_target_one_charge, 1e-12);
    EXPECT_NEAR(*target_one.samples.front().gaussian_response,
        *target_one.samples.front().peeled_response - expected_target_one_charge, 1e-12);
    EXPECT_GT(std::abs(expected_target_one_charge), 1e-8);

    const auto target_two = context.BuildPostJointLocalSamples(3, peelings.at(3));
    ASSERT_EQ(target_two.samples.size(), 1u);
    ASSERT_TRUE(target_two.samples.front().peeled_response);
    ASSERT_TRUE(target_two.samples.front().fixed_charge_response);
    ASSERT_TRUE(target_two.samples.front().gaussian_response);
    std::vector<double> target_two_charge(original.row_ids.size(), 0.0);
    for (const auto & support : original.support[2])
        target_two_charge[support.row] = -.45 * core::joint_component::EvaluateKernel(
            support.squared_distance, 1.15, 2.5).charge;
    const double expected_target_two_charge = ReferenceTricubicCharge(
        fixture.map, target_two_charge, {6.25, 0.0, 0.0});
    EXPECT_NEAR(*target_two.samples.front().fixed_charge_response,
        expected_target_two_charge, 1e-12);
    EXPECT_NEAR(*target_two.samples.front().gaussian_response,
        *target_two.samples.front().peeled_response - expected_target_two_charge, 1e-12);
    EXPECT_LT(expected_target_two_charge, -1e-8);

    auto zero_charge_result = result;
    zero_charge_result.components.front().state->ac[1] = 0.0;
    data_internal::ApplyJointStageEstimates(fixture.model, zero_charge_result, "indexed-charge-run");
    const core::detail::JointPostprocessingContext zero_context(
        fixture.map, fixture.model, fixture.problem, &zero_charge_result);
    const auto zero_peelings = zero_context.BuildPostFitPeelingSamples(std::vector<std::size_t>{1});
    const auto zero_target = zero_context.BuildPostJointLocalSamples(2, zero_peelings.at(2));
    ASSERT_TRUE(zero_target.samples.front().fixed_charge_response);
    EXPECT_DOUBLE_EQ(*zero_target.samples.front().fixed_charge_response, 0.0);

    auto changed_halo_result = result;
    changed_halo_result.halo_median_reconstruction->global_median[2] = 1.7;
    changed_halo_result.halo_median_reconstruction->halos.front().charge = 1.7;
    data_internal::ApplyJointStageEstimates(fixture.model, changed_halo_result, "indexed-charge-run");
    const core::detail::JointPostprocessingContext changed_halo_context(
        fixture.map, fixture.model, fixture.problem, &changed_halo_result);
    const auto changed_halo_peelings =
        changed_halo_context.BuildPostFitPeelingSamples(std::vector<std::size_t>{1});
    const auto changed_halo_target = changed_halo_context.BuildPostJointLocalSamples(
        2, changed_halo_peelings.at(2));
    ASSERT_TRUE(changed_halo_target.samples.front().fixed_charge_response);
    EXPECT_DOUBLE_EQ(*changed_halo_target.samples.front().fixed_charge_response,
        *target_one.samples.front().fixed_charge_response);
}

TEST(PureHaloSupportExclusionTest, ClassifiesPeelingGeometrySupportAndMaskIntersections)
{
    auto fixture=WideBoundaryFixture();
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    auto snapshot=PurePeelingSnapshot(problem,fixture.model);
    data_internal::ApplyJointStageEstimates(fixture.model,snapshot,"peeling-domain-audit");
    const SamplingPointList points{
        {1.75,{1.75,0,0},true},
        {1.9,{1.9,0,0},true},
        {10.0,{10.0,0,0},true}};
    fixture.model.EditAnalysis().SetAtomLocalRawSamplingEntries(*fixture.model.FindAtomPtr(1),
        core::detail::BuildLocalPotentialSampleList(fixture.map,points));
    core::detail::JointPostprocessingMetrics metrics;
    const core::detail::JointPostprocessingContext context(
        fixture.map,fixture.model,problem,&snapshot,&metrics);
    const auto peelings=context.BuildPostFitPeelingSamples(std::vector<std::size_t>{0});
    const auto & samples=peelings.at(1).samples;
    ASSERT_EQ(samples.size(),3u);
    EXPECT_TRUE(samples[0].response);
    EXPECT_TRUE(samples[1].response);
    EXPECT_TRUE(samples[2].response);
    EXPECT_EQ(metrics.raw_peeling_sample_count,3u);
    EXPECT_EQ(metrics.available_peeling_sample_count,3u);
    EXPECT_EQ(metrics.unavailable_peeling_sample_count,0u);
    EXPECT_EQ(metrics.map_geometry_boundary_samples,1u);
    EXPECT_EQ(metrics.mapping_mismatch_samples,0u);
    EXPECT_GE(metrics.fitting_mask_only_intersection_samples,1u);
    EXPECT_GT(metrics.nonzero_stencil_node_count,0u);
    EXPECT_GT(metrics.map_geometry_boundary_nodes,0u);
    EXPECT_GT(metrics.original_support_boundary_nodes,0u);
    EXPECT_GT(metrics.fitting_mask_only_intersection_nodes,0u);

    fixture.model.EditAnalysis().SetAtomPostFitPeeling(*fixture.model.FindAtomPtr(1),peelings.at(1));
    const auto assessment=AtomLocalPotentialView::For(*fixture.model.FindAtomPtr(1))
        .AssessLocalFittingPeelingRatio(1.0,2.0);
    EXPECT_TRUE(assessment.ratio);
    EXPECT_EQ(assessment.total_samples,2u);
    EXPECT_EQ(assessment.valid_samples,2u);
    EXPECT_TRUE(assessment.reason.empty());
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
        for (const auto & halo:halo_positions) halo_field[grid]+=contribution(halo,point);
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
        const double direct_targets=ReferenceTricubicCharge(
            fixture.map,target_field,points[sample].position)+ReferenceTricubicCharge(
                fixture.map,other_target_field,points[sample].position);
        EXPECT_NEAR(*peeling.samples[sample].response,direct_targets,1e-12);
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

TEST(PureHaloSupportExclusionTest, CurrentReconstructionOnlyPolicyChangesPredictionWhenExtraContributorGroupChanges)
{
    const auto peel_sample=[](MapFixture & fixture) {
        const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
        auto snapshot=QualifiedReconstructionSnapshot(problem,fixture.model);
        data_internal::ApplyJointStageEstimates(fixture.model,snapshot,"reconstruction-only-run");
        const SamplingPointList points{{2.8,{2.8,0.0,0.0},true}};
        fixture.model.EditAnalysis().SetAtomLocalRawSamplingEntries(*fixture.model.FindAtomPtr(1),
            core::detail::BuildLocalPotentialSampleList(fixture.map,points));
        const core::detail::JointPostprocessingContext context(
            fixture.map,fixture.model,problem,&snapshot);
        return context.BuildPostFitPeelingSamples(std::vector<std::size_t>{0})
            .at(1).samples.front().response;
    };
    auto baseline=ReconstructionOnlyFixture();
    auto relabeled=ReconstructionOnlyFixture();
    for (const int id:{2,3})
    {
        relabeled.model.FindAtomPtr(id)->SetComponentKey(2);
        relabeled.model.FindAtomPtr(id)->SetAtomKey(9);
    }

    const auto original_response=peel_sample(baseline);
    const auto relabeled_response=peel_sample(relabeled);
    ASSERT_TRUE(original_response);
    ASSERT_TRUE(relabeled_response);
    EXPECT_GT(std::abs(*original_response-*relabeled_response),1e-6);
}

TEST(PureHaloSupportExclusionTest, MissingHaloMedianKeepsOffDomainContributorUnavailable)
{
    auto fixture=ReconstructionOnlyFixture();
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    auto snapshot=PurePeelingSnapshot(problem,fixture.model);
    data_internal::ApplyJointStageEstimates(fixture.model,snapshot,"no-median-run");
    const SamplingPointList points{{1.9,{1.9,0.0,0.0},true}};
    fixture.model.EditAnalysis().SetAtomLocalRawSamplingEntries(*fixture.model.FindAtomPtr(1),
        core::detail::BuildLocalPotentialSampleList(fixture.map,points));
    core::detail::JointPostprocessingMetrics metrics;
    const core::detail::JointPostprocessingContext context(
        fixture.map,fixture.model,problem,&snapshot,&metrics);
    const auto peelings=context.BuildPostFitPeelingSamples(std::vector<std::size_t>{0});
    ASSERT_EQ(peelings.at(1).samples.size(),1u);
    EXPECT_FALSE(peelings.at(1).samples.front().response);
    EXPECT_EQ(peelings.at(1).samples.front().reason,"missing-contributor-state");
    EXPECT_EQ(metrics.missing_contributor_samples,1u);
}

TEST(PureHaloSupportExclusionTest, RejectsReconstructionMedianFromDifferentJointEndpoint)
{
    auto fixture=ReconstructionOnlyFixture();
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    auto snapshot=QualifiedReconstructionSnapshot(problem,fixture.model);
    snapshot.halo_median_reconstruction=core::detail::BuildHaloMedianReconstruction(
        fixture.model,problem,snapshot,"reconstruction-only-run").snapshot;
    snapshot.halo_median_reconstruction->global_median[0]+=1.0;
    data_internal::ApplyJointStageEstimates(fixture.model,snapshot,"reconstruction-only-run");
    const core::detail::JointPostprocessingContext context(
        fixture.map,fixture.model,problem,&snapshot);
    EXPECT_THROW(context.BuildPostFitPeelingSamples(std::vector<std::size_t>{0}),std::invalid_argument);
}

TEST(PureHaloSupportExclusionTest, BoundsReconstructionNodeCacheForLargeSamplingSet)
{
    auto fixture=ReconstructionOnlyFixture();
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    auto snapshot=QualifiedReconstructionSnapshot(problem,fixture.model);
    data_internal::ApplyJointStageEstimates(fixture.model,snapshot,"reconstruction-only-run");

    SamplingPointList points;
    constexpr std::size_t sample_axis_count=25;
    points.reserve(sample_axis_count*sample_axis_count*sample_axis_count);
    const auto dimensions=fixture.map.GetGridSize();
    for (std::size_t z=0;z<sample_axis_count;++z)
        for (std::size_t y=0;y<sample_axis_count;++y)
            for (std::size_t x=0;x<sample_axis_count;++x)
    {
        const std::array<int,3> index{
            static_cast<int>(x*static_cast<std::size_t>(dimensions[0]-1)/(sample_axis_count-1)),
            static_cast<int>(y*static_cast<std::size_t>(dimensions[1]-1)/(sample_axis_count-1)),
            static_cast<int>(z*static_cast<std::size_t>(dimensions[2]-1)/(sample_axis_count-1))};
        const auto position=core::simulation::GridPosition(
            index,fixture.map.GetGridSpacing(),fixture.map.GetOrigin());
        points.push_back({std::sqrt(position[0]*position[0]+position[1]*position[1]+position[2]*position[2]),
            position,true});
    }
    fixture.model.EditAnalysis().SetAtomLocalRawSamplingEntries(*fixture.model.FindAtomPtr(1),
        core::detail::BuildLocalPotentialSampleList(fixture.map,points));
    core::detail::JointPostprocessingMetrics metrics;
    const core::detail::JointPostprocessingContext context(
        fixture.map,fixture.model,problem,&snapshot,&metrics);
    std::map<int,PostFitPeelingResult> peelings;
    ASSERT_NO_THROW(peelings=context.BuildPostFitPeelingSamples(std::vector<std::size_t>{0}));
    ASSERT_TRUE(peelings.contains(1));
    ASSERT_EQ(peelings.at(1).samples.size(),points.size());
    EXPECT_EQ(metrics.raw_peeling_sample_count,points.size());
    EXPECT_EQ(metrics.reconstruction_node_cache_peak_entries,4096u);
    EXPECT_GT(metrics.reconstruction_only_node_evaluation_count,4096u);
    EXPECT_GT(metrics.estimated_reconstruction_workspace_bytes,0u);
    EXPECT_GE(metrics.estimated_peak_postprocessing_payload_bytes,
        metrics.estimated_reconstruction_workspace_bytes);
}

TEST(PureHaloSupportExclusionTest, PeelingRejectsFitReconstructionIdentityMismatch)
{
    auto fixture=BoundaryFixture();
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    auto snapshot=PurePeelingSnapshot(problem,fixture.model);
    snapshot.fitting_row_ids.front()="stale";
    EXPECT_THROW((core::detail::JointPostprocessingContext(fixture.map,fixture.model,problem,&snapshot)),
        std::invalid_argument);
}

TEST(PureHaloSupportExclusionTest, PersistsVersionedDualDomainAndMedianProvenance)
{
    auto fixture=BoundaryFixture();
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    const auto snapshot=PurePeelingSnapshot(problem,fixture.model);
    const auto encoded=joint_result_io::Encode(snapshot);
    EXPECT_EQ(boost::json::parse(encoded).at("schema_version").as_int64(),17);
    const auto decoded=joint_result_io::Decode(encoded);
    EXPECT_EQ(decoded.parameterization_contract,JointPureHaloSupportExclusionParameterizationContract);
    EXPECT_EQ(decoded.fitting_to_original_contributors,problem.FittingToOriginalContributors());
    ASSERT_TRUE(decoded.halo_median_reconstruction);
    EXPECT_EQ(decoded.halo_median_reconstruction->halos,snapshot.halo_median_reconstruction->halos);
}

TEST(PureHaloSupportExclusionTest, PersistsUnqualifiedFailureDiagnosticsWithoutFormalEstimate)
{
    auto fixture=BoundaryFixture();
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    auto result=PurePeelingSnapshot(problem,fixture.model);
    result.initialization.valid=false;
    result.initialization.reason="target-width-initialization-failed";
    result.initialization.b.clear();
    result.initialization.atoms.clear();
    result.search_completed=false;
    result.runtime_convergence=JointCheckStatus::Failed;
    result.target_runtime_convergence=JointCheckStatus::NotRun;
    result.assembled_state.reset();
    result.objective.reset();
    result.halo_median_reconstruction.reset();
    result.components.front().state.reset();
    result.components.front().search_completed=false;
    result.components.front().runtime_convergence=JointCheckStatus::Failed;
    result.components.front().stop_reason="block-recovery-exhausted";

    const auto decoded=joint_result_io::Decode(joint_result_io::Encode(result));
    EXPECT_FALSE(IsQualifiedJointAnalysisResult(decoded));
    EXPECT_EQ(decoded.initialization.reason,"target-width-initialization-failed");
    EXPECT_EQ(decoded.components.front().stop_reason,"block-recovery-exhausted");
    EXPECT_FALSE(decoded.components.front().state);
    EXPECT_FALSE(decoded.halo_median_reconstruction);
}

TEST(PureHaloSupportExclusionTest, PersistsSearchFailureWithValidInitializationWithoutEstimate)
{
    auto fixture=BoundaryFixture();
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    auto result=QualifiedPurePeelingSnapshot(problem,fixture.model);
    result.search_completed=false;
    result.runtime_convergence=JointCheckStatus::Failed;
    result.target_runtime_convergence=JointCheckStatus::NotRun;
    result.assembled_state.reset();
    result.objective.reset();
    result.halo_median_reconstruction.reset();
    result.components.front().state.reset();
    result.components.front().search_completed=false;
    result.components.front().runtime_convergence=JointCheckStatus::Failed;
    result.components.front().target_runtime_convergence=JointCheckStatus::NotRun;
    result.components.front().stop_reason="block-recovery-exhausted";

    const auto decoded=joint_result_io::Decode(joint_result_io::Encode(result));
    EXPECT_TRUE(decoded.initialization.valid);
    EXPECT_FALSE(IsQualifiedJointAnalysisResult(decoded));
    EXPECT_EQ(decoded.components.front().stop_reason,"block-recovery-exhausted");
    EXPECT_FALSE(decoded.components.front().state);
    EXPECT_FALSE(decoded.halo_median_reconstruction);
}

TEST(PureHaloSupportExclusionTest, PersistsQualifiedNoDonorAsExplicitDownstreamUnavailable)
{
    auto fixture=BoundaryFixture();
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    auto result=QualifiedPurePeelingSnapshot(problem,fixture.model);
    result.halo_median_reconstruction.reset();
    result.evidence.push_back({std::string(JointHaloMedianReconstructionEvidenceName),
        JointCheckStatus::Unavailable,JointEvidenceScope::AssembledGlobal,{}, {},
        "halo-median-no-valid-targets"});

    ASSERT_TRUE(IsQualifiedJointAnalysisResult(result));
    ASSERT_TRUE(HasUnavailableJointHaloMedianReconstruction(result));
    const auto decoded=joint_result_io::Decode(joint_result_io::Encode(result));
    EXPECT_TRUE(IsQualifiedJointAnalysisResult(decoded));
    EXPECT_TRUE(HasUnavailableJointHaloMedianReconstruction(decoded));
    EXPECT_FALSE(decoded.halo_median_reconstruction);
}

TEST(PureHaloSupportExclusionTest, RejectsQualifiedHaloResultWithoutMedianOrUnavailableEvidence)
{
    auto fixture=BoundaryFixture();
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    auto result=QualifiedPurePeelingSnapshot(problem,fixture.model);
    result.halo_median_reconstruction.reset();
    EXPECT_THROW(joint_result_io::Encode(result),std::invalid_argument);
}

TEST(PureHaloSupportExclusionTest, RunsFullPureWorkflowWithGroupAndGlobalHaloSources)
{
    auto fixture=WorkflowFixture();
    FitOptions options;
    options.estimator=PotentialEstimator::JOINT_COMPONENTS;
    options.quiet_mode=true;
    options.thread_size=1;
    core::detail::JointPostprocessingMetrics metrics;
    ScopedInternalUnboundedOracle unbounded_test_solver;
    ASSERT_NO_THROW(core::detail::RunJointComponentWorkflow(
        fixture.map,fixture.model,options,&metrics));
    const auto result=fixture.model.GetAnalysisView().GetJointResult();
    ASSERT_TRUE(result);
    EXPECT_EQ(result->parameterization_contract,JointPureHaloSupportExclusionParameterizationContract);
    ASSERT_TRUE(result->halo_median_reconstruction);
    EXPECT_EQ(result->halo_median_reconstruction->halos.size(),2u);
    EXPECT_EQ(metrics.halo_median_valid_target_count,1u);
    EXPECT_EQ(metrics.halo_median_same_group_halos,1u);
    EXPECT_EQ(metrics.halo_median_fallback_halos,1u);
    EXPECT_GT(metrics.halo_median_fallback_support_rows,0u);
    const auto target=AtomLocalPotentialView::For(*fixture.model.FindAtomPtr(1));
    EXPECT_TRUE(target.GetPostFitPeeling());
    EXPECT_TRUE(target.GetPostJointLocalResult());
    EXPECT_EQ(target.GetStageEstimate(FittingStage::Second).source.role,FittingRole::Target);
    const auto halo=AtomLocalPotentialView::For(*fixture.model.FindAtomPtr(2));
    EXPECT_EQ(halo.GetStageEstimate(FittingStage::Second).source.role,FittingRole::Halo);
    EXPECT_FALSE(halo.GetStageEstimate(FittingStage::Second).point);
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
    EXPECT_EQ(metrics.halo_median_fallback_halos,0u);
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

TEST(PureHaloSupportExclusionTest, RejectsUnqualifiedJointSnapshotWithCurrentJointDownstreamState)
{
    auto fixture=WorkflowFixture();
    FitOptions options;
    options.estimator=PotentialEstimator::JOINT_COMPONENTS;
    options.quiet_mode=true;
    options.thread_size=1;
    ScopedInternalUnboundedOracle unbounded_test_solver;
    ASSERT_NO_THROW(core::detail::RunJointComponentWorkflow(
        fixture.map,fixture.model,options,nullptr));
    ASSERT_TRUE(fixture.model.GetAnalysisView().GetJointResult());

    auto stale=*fixture.model.GetAnalysisView().GetJointResult();
    stale.search_completed=false;
    stale.runtime_convergence=JointCheckStatus::Failed;
    stale.target_runtime_convergence=JointCheckStatus::NotRun;
    stale.halo_median_reconstruction.reset();
    for (auto & component:stale.components)
    {
        component.search_completed=false;
        component.runtime_convergence=JointCheckStatus::Failed;
        component.target_runtime_convergence=JointCheckStatus::NotRun;
    }
    fixture.model.EditAnalysis().SetJointResult(std::move(stale));

    EXPECT_THROW(stage_result_io::Encode(fixture.model),std::invalid_argument);
}

TEST(PureHaloSupportExclusionTest, RejectsMixedCurrentJointRunIdentities)
{
    auto fixture=WorkflowFixture();
    FitOptions options;
    options.estimator=PotentialEstimator::JOINT_COMPONENTS;
    options.quiet_mode=true;
    options.thread_size=1;
    ScopedInternalUnboundedOracle unbounded_test_solver;
    ASSERT_NO_THROW(core::detail::RunJointComponentWorkflow(
        fixture.map,fixture.model,options,nullptr));
    ASSERT_TRUE(fixture.model.GetAnalysisView().GetJointResult());
    auto snapshot=*fixture.model.GetAnalysisView().GetJointResult();
    snapshot.halo_median_reconstruction.reset();
    const auto * atom=fixture.model.FindAtomPtr(1);
    auto changed=AtomLocalPotentialView::For(*atom).GetStageEstimate(FittingStage::Second);
    changed.source.run_id="different-joint-run";
    fixture.model.EditAnalysis().SetAtomStageEstimate(FittingStage::Second,*atom,changed);
    fixture.model.EditAnalysis().SetJointResult(std::move(snapshot));

    EXPECT_THROW(stage_result_io::Encode(fixture.model),std::invalid_argument);
}

TEST(PureHaloSupportExclusionTest, FullWorkflowRetainsMultipleMemberGroupContract)
{
    auto fixture=WorkflowFixture();
    fixture.model.SelectAllAtoms();
    fixture.model.FindAtomPtr(3)->SetComponentKey(static_cast<ComponentKey>(1));
    fixture.model.FindAtomPtr(3)->SetAtomKey(static_cast<AtomKey>(7));
    FitOptions options;
    options.estimator=PotentialEstimator::JOINT_COMPONENTS;
    options.quiet_mode=true;
    options.thread_size=1;
    core::detail::JointPostprocessingMetrics metrics;
    ScopedInternalUnboundedOracle unbounded_test_solver;
    ASSERT_NO_THROW(core::detail::RunJointComponentWorkflow(
        fixture.map,fixture.model,options,&metrics));
    const auto result=fixture.model.GetAnalysisView().GetJointResult();
    ASSERT_TRUE(result);
    ASSERT_TRUE(result->selection_domain);
    EXPECT_EQ(result->selection_domain->target_indices.size(),3u);
    EXPECT_FALSE(result->halo_median_reconstruction);
    const auto keys=fixture.model.GetAnalysisView().CollectAtomGroupKeys();
    bool found_multi_member=false;
    for (const auto key:keys)
    {
        const auto summary=fixture.model.GetAnalysisView().GetGroupParameterSummary(key);
        if (!summary) continue;
        const auto & members=fixture.model.GetAnalysisView().GetAtomObjectList(key);
        if (members.size()!=3) continue;
        found_multi_member=true;
        EXPECT_EQ(summary->member_ids,(std::vector<int>{1,2,3}));
        EXPECT_EQ(summary->point_count,summary->eligible_count);
        EXPECT_EQ(summary->excluded_count,0u);
        EXPECT_EQ(summary->status,EvidenceStatus::Available) << summary->reason;
        ASSERT_TRUE(summary->inference);
        EXPECT_TRUE(std::isfinite(summary->inference->mdpde.GetAmplitude()));
        EXPECT_TRUE(std::isfinite(summary->inference->mdpde.GetWidth()));
    }
    EXPECT_TRUE(found_multi_member);
    EXPECT_EQ(metrics.halo_median_same_group_halos,0u);
    EXPECT_EQ(metrics.halo_median_fallback_halos,0u);
}

TEST(PureHaloSupportExclusionTest, RunsFullPureWorkflowWithMultipleTargetsAndHaloGroupSources)
{
    auto fixture=MultiTargetHaloWorkflowFixture();
    FitOptions options;
    options.estimator=PotentialEstimator::JOINT_COMPONENTS;
    options.quiet_mode=true;
    options.thread_size=1;
    core::detail::JointPostprocessingMetrics metrics;
    ScopedInternalUnboundedOracle unbounded_test_solver;
    ASSERT_NO_THROW(core::detail::RunJointComponentWorkflow(
        fixture.map,fixture.model,options,&metrics));

    const auto result=fixture.model.GetAnalysisView().GetJointResult();
    ASSERT_TRUE(result);
    ASSERT_TRUE(result->selection_domain);
    EXPECT_EQ(result->selection_domain->target_indices,(std::vector<std::size_t>{0,1,2}));
    ASSERT_TRUE(result->halo_median_reconstruction);
    ASSERT_EQ(result->halo_median_reconstruction->halos.size(),2u);
    EXPECT_EQ(result->halo_median_reconstruction->valid_target_count,3u);
    EXPECT_EQ(metrics.halo_median_valid_target_count,3u);
    EXPECT_EQ(metrics.halo_median_same_group_halos,1u);
    EXPECT_EQ(metrics.halo_median_fallback_halos,1u);

    const auto same_group_halo=std::find_if(result->halo_median_reconstruction->halos.begin(),
        result->halo_median_reconstruction->halos.end(),
        [](const auto & halo){return halo.atom_id=="4";});
    const auto fallback_halo=std::find_if(result->halo_median_reconstruction->halos.begin(),
        result->halo_median_reconstruction->halos.end(),
        [](const auto & halo){return halo.atom_id=="5";});
    ASSERT_NE(same_group_halo,result->halo_median_reconstruction->halos.end());
    ASSERT_NE(fallback_halo,result->halo_median_reconstruction->halos.end());
    EXPECT_EQ(same_group_halo->source,HaloMedianSource::GroupKeyMedian);
    EXPECT_EQ(same_group_halo->donor_count,3u);
    EXPECT_EQ(fallback_halo->source,HaloMedianSource::GlobalTargetMedian);
    EXPECT_EQ(fallback_halo->donor_count,3u);
    EXPECT_TRUE(std::isfinite(same_group_halo->amplitude));
    EXPECT_TRUE(std::isfinite(fallback_halo->amplitude));

    double local_amplitude_sum{};
    double local_width_sum{};
    const std::array<double,3> expected_amplitudes{{2.0,1.6,1.2}};
    const std::array<double,3> expected_widths{{.65,.75,.85}};
    for (std::size_t target=0;target<3;++target)
    {
        const int atom_id=static_cast<int>(target+1);
        const auto local=AtomLocalPotentialView::For(*fixture.model.FindAtomPtr(atom_id)).GetPostJointLocalResult();
        ASSERT_TRUE(local);
        EXPECT_EQ(local->status,EvidenceStatus::Available) << local->reason;
        EXPECT_GT(local->sample_count,0u);
        EXPECT_TRUE(std::isfinite(local->fit.mdpde.GetModel().GetAmplitude()));
        EXPECT_TRUE(std::isfinite(local->fit.mdpde.GetModel().GetWidth()));
        EXPECT_TRUE(std::isfinite(local->fixed_charge.coefficient));
        EXPECT_NEAR(local->fit.mdpde.GetModel().GetAmplitude(),expected_amplitudes[target],.25);
        EXPECT_NEAR(local->fit.mdpde.GetModel().GetWidth(),expected_widths[target],.08);
        local_amplitude_sum+=local->fit.mdpde.GetModel().GetAmplitude();
        local_width_sum+=local->fit.mdpde.GetModel().GetWidth();
    }

    const auto group_key=data_internal::GetGroupKey(fixture.model.FindAtomPtr(1));
    const auto summary=fixture.model.GetAnalysisView().GetGroupParameterSummary(group_key);
    ASSERT_TRUE(summary);
    EXPECT_EQ(summary->member_ids,(std::vector<int>{1,2,3}));
    EXPECT_EQ(summary->point_count,3u);
    EXPECT_EQ(summary->eligible_count,3u);
    EXPECT_EQ(summary->excluded_count,0u);
    EXPECT_EQ(summary->status,EvidenceStatus::Available) << summary->reason;
    ASSERT_TRUE(summary->descriptive_mean);
    EXPECT_NEAR(summary->descriptive_mean->GetAmplitude(),local_amplitude_sum/3.0,1e-12);
    EXPECT_NEAR(summary->descriptive_mean->GetWidth(),local_width_sum/3.0,1e-12);
    ASSERT_TRUE(summary->inference);
    EXPECT_TRUE(std::isfinite(summary->inference->mdpde.GetAmplitude()));
    EXPECT_TRUE(std::isfinite(summary->inference->mdpde.GetWidth()));
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

namespace {
struct Phase3Target
{
    int id{};
    GaussianModel3D point;
};

struct Phase3Halo
{
    std::size_t contributor{};
    int id{};
    GaussianModel3D truth;
};

struct Phase3Row
{
    std::string policy;
    std::vector<GaussianModel3D> halos;
    std::vector<std::size_t> donor_counts;
};

double ComponentMedian(std::vector<double> values)
{
    std::sort(values.begin(),values.end());
    const auto middle=values.size()/2;
    return values.size()%2 ? values[middle] : (values[middle-1]+values[middle])/2.0;
}

GaussianModel3D MedianPoint(const std::vector<Phase3Target> & donors)
{
    std::vector<double> amplitudes,widths,charges;
    for (const auto & donor:donors)
    {
        amplitudes.push_back(donor.point.GetAmplitude());
        widths.push_back(donor.point.GetWidth());
        charges.push_back(donor.point.GetOffset());
    }
    return {ComponentMedian(std::move(amplitudes)),ComponentMedian(std::move(widths)),
        ComponentMedian(std::move(charges))};
}

std::vector<Phase3Target> Phase3JointTargets(
    const JointAnalysisResult & result)
{
    std::vector<Phase3Target> targets;
    if (!result.selection_domain) return targets;
    for (const auto index:result.selection_domain->target_indices)
    {
        std::optional<GaussianModel3D> point;
        for (const auto & component:result.components)
        {
            if (!component.state) continue;
            const auto & atoms=component.layout ? component.layout->full_atoms : component.atoms;
            const auto found=std::find(atoms.begin(),atoms.end(),index);
            if (found==atoms.end()) continue;
            const auto offset=static_cast<std::size_t>(std::distance(atoms.begin(),found));
            point.emplace(component.state->ac[2*offset],component.state->b[offset],
                component.state->ac[2*offset+1]);
            break;
        }
        if (!point) continue;
        const auto id=std::stoi(result.atom_ids[index]);
        targets.push_back({id,*point});
    }
    return targets;
}

std::vector<GaussianModel3D> Phase3KnnProfiles(
    const ModelObject & model,const std::vector<Phase3Target> & targets,
    const std::vector<JointHaloMedianReconstruction> & halos,std::size_t k)
{
    std::vector<GaussianModel3D> profiles;
    for (const auto & halo:halos)
    {
        const auto halo_position=model.FindAtomPtr(std::stoi(halo.atom_id))->GetPosition();
        std::vector<std::pair<double,Phase3Target>> ordered;
        for (const auto & target:targets)
        {
            const auto target_position=model.FindAtomPtr(target.id)->GetPosition();
            double square{};
            for (std::size_t axis=0;axis<3;++axis)
            {
                const auto delta=target_position[axis]-halo_position[axis];
                square+=delta*delta;
            }
            ordered.push_back({square,target});
        }
        std::sort(ordered.begin(),ordered.end(),[](const auto & left,const auto & right) {
            if (left.first!=right.first) return left.first<right.first;
            return left.second.id<right.second.id;
        });
        std::vector<Phase3Target> donors;
        for (std::size_t index=0;index<std::min(k,ordered.size());++index)
            donors.push_back(ordered[index].second);
        profiles.push_back(MedianPoint(donors));
    }
    return profiles;
}

std::optional<std::vector<GaussianModel3D>> Phase3ObservationAC(
    const JointProblem & problem,
    const std::vector<Phase3Target> & targets,
    const std::vector<JointHaloMedianReconstruction> & halos,
    const std::vector<GaussianModel3D> & width_profiles)
{
    const auto & input=problem.OriginalInput();
    std::map<int,GaussianModel3D> target_by_id;
    for (const auto & target:targets) target_by_id.emplace(target.id,target.point);
    std::vector<double> target_prediction(input.observations.size());
    for (std::size_t atom=0;atom<input.atom_ids.size();++atom)
    {
        const auto found=target_by_id.find(std::stoi(input.atom_ids[atom]));
        if (found==target_by_id.end()) continue;
        for (const auto & support:input.support[atom])
        {
            const auto basis=core::joint_component::EvaluateKernel(support.squared_distance,
                found->second.GetWidth(),2.5);
            target_prediction[support.row]+=found->second.GetAmplitude()*basis.gaussian+
                found->second.GetOffset()*basis.charge;
        }
    }
    std::vector<std::size_t> halo_indices;
    for (std::size_t atom=0;atom<input.atom_ids.size();++atom)
        if (!target_by_id.contains(std::stoi(input.atom_ids[atom]))) halo_indices.push_back(atom);
    if (halo_indices.size()!=halos.size() || halos.size()!=width_profiles.size() || halos.empty())
        return std::nullopt;

    Eigen::MatrixXd design=Eigen::MatrixXd::Zero(
        static_cast<Eigen::Index>(input.observations.size()),static_cast<Eigen::Index>(2*halos.size()));
    for (std::size_t h=0;h<halo_indices.size();++h)
        for (const auto & support:input.support[halo_indices[h]])
        {
            const auto basis=core::joint_component::EvaluateKernel(support.squared_distance,
                width_profiles[h].GetWidth(),2.5);
            design(static_cast<Eigen::Index>(support.row),static_cast<Eigen::Index>(2*h))=basis.gaussian;
            design(static_cast<Eigen::Index>(support.row),static_cast<Eigen::Index>(2*h+1))=basis.charge;
        }
    Eigen::VectorXd response(static_cast<Eigen::Index>(input.observations.size()));
    for (std::size_t row=0;row<input.observations.size();++row)
        response[static_cast<Eigen::Index>(row)]=input.observations[row]-target_prediction[row];

    double best_error=std::numeric_limits<double>::infinity();
    std::optional<Eigen::VectorXd> best;
    const auto active_sets=std::size_t{1}<<halos.size();
    for (std::size_t active=0;active<active_sets;++active)
    {
        std::vector<std::pair<std::size_t,bool>> columns;
        for (std::size_t h=0;h<halos.size();++h)
            if (active&(std::size_t{1}<<h)) columns.push_back({h,false});
        for (std::size_t h=0;h<halos.size();++h) columns.push_back({h,true});
        Eigen::MatrixXd reduced(design.rows(),static_cast<Eigen::Index>(columns.size()));
        for (std::size_t column=0;column<columns.size();++column)
            reduced.col(static_cast<Eigen::Index>(column))=design.col(static_cast<Eigen::Index>(
                2*columns[column].first+(columns[column].second ? 1 : 0)));
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solver(reduced);
        solver.setThreshold(1.0e-10);
        if (solver.rank()!=reduced.cols()) continue;
        const auto estimate=solver.solve(response);
        bool feasible=true;
        for (std::size_t column=0;column<columns.size();++column)
            if (!columns[column].second && estimate[static_cast<Eigen::Index>(column)]<-1.0e-10)
                feasible=false;
        if (!feasible) continue;
        const double error=(reduced*estimate-response).squaredNorm();
        if (error<best_error)
        {
            best_error=error;
            Eigen::VectorXd full=Eigen::VectorXd::Zero(static_cast<Eigen::Index>(2*halos.size()));
            for (std::size_t column=0;column<columns.size();++column)
                full[static_cast<Eigen::Index>(2*columns[column].first+(columns[column].second ? 1 : 0))]=
                    estimate[static_cast<Eigen::Index>(column)];
            best=std::move(full);
        }
    }
    if (!best) return std::nullopt;
    std::vector<GaussianModel3D> profiles;
    for (std::size_t h=0;h<halos.size();++h)
        profiles.emplace_back(std::max(0.0,(*best)[static_cast<Eigen::Index>(2*h)]),
            width_profiles[h].GetWidth(),(*best)[static_cast<Eigen::Index>(2*h+1)]);
    return profiles;
}

std::vector<double> Phase3TargetTruthField(const MapObject & map)
{
    const std::array<int,3> ids{{1,2,3}};
    const std::array<GaussianModel3D,3> truth{{{2.0,.65,-.2},{1.6,.75,.15},{1.2,.85,-.1}}};
    std::vector<double> field(map.GetMapValueArraySize());
    for (std::size_t grid=0;grid<field.size();++grid)
    {
        const auto position=map.GetGridPosition(grid);
        for (std::size_t atom=0;atom<ids.size();++atom)
        {
            const std::array<std::array<double,3>,3> centers{{{0,0,0},{3.5,0,0},{7,0,0}}};
            double square{};
            for (std::size_t axis=0;axis<3;++axis)
            {
                const auto delta=position[axis]-centers[atom][axis];
                square+=delta*delta;
            }
            const auto basis=core::joint_component::EvaluateKernel(square,truth[atom].GetWidth(),2.5);
            field[grid]+=truth[atom].GetAmplitude()*basis.gaussian+truth[atom].GetOffset()*basis.charge;
        }
    }
    return field;
}

template<typename T>
void Phase3CsvValue(const std::optional<T> & value)
{
    if (value) std::cout<<*value;
    else std::cout<<"NA";
}
}

TEST(PureHaloSupportExclusionTest, ComparesLabelFreePoliciesThroughPeelingAndLocalGroupFitting)
{
    auto geometry_fixture=MultiTargetHaloWorkflowFixture();
    auto noisy_values=std::make_unique<double[]>(geometry_fixture.map.GetMapValueArraySize());
    for (std::size_t grid=0;grid<geometry_fixture.map.GetMapValueArraySize();++grid)
        noisy_values[grid]=geometry_fixture.map.GetMapValue(grid)+.01*std::sin(.127*static_cast<double>(grid+1))+
            .005*std::cos(.071*static_cast<double>(grid+3));
    MapObject noisy_map(geometry_fixture.map.GetGridSize(),geometry_fixture.map.GetGridSpacing(),
        geometry_fixture.map.GetOrigin(),std::move(noisy_values));
    MapFixture fixture{std::move(geometry_fixture.model),std::move(noisy_map)};

    FitOptions options;
    options.estimator=PotentialEstimator::JOINT_COMPONENTS;
    options.quiet_mode=true;
    options.thread_size=1;
    ScopedInternalUnboundedOracle unbounded_test_solver;
    ASSERT_NO_THROW(core::detail::RunJointComponentWorkflow(fixture.map,fixture.model,options));
    auto result=fixture.model.GetAnalysisView().GetJointResult();
    ASSERT_TRUE(result);
    ASSERT_TRUE(result->selection_domain);
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    const auto targets=Phase3JointTargets(*result);
    ASSERT_EQ(targets.size(),3u);
    ASSERT_TRUE(result->halo_median_reconstruction);
    const auto baseline=*result->halo_median_reconstruction;
    ASSERT_EQ(baseline.halos.size(),2u);
    const auto run_id=AtomLocalPotentialView::For(*fixture.model.FindAtomPtr(1))
        .GetStageEstimate(FittingStage::Second).source.run_id;
    const auto baseline_build=core::detail::BuildHaloMedianReconstruction(
        fixture.model,problem,*result,run_id);
    const auto global=GaussianModel3D{baseline.global_median[0],baseline.global_median[1],baseline.global_median[2]};
    const auto knn1=Phase3KnnProfiles(fixture.model,targets,baseline.halos,1);
    const auto knn2=Phase3KnnProfiles(fixture.model,targets,baseline.halos,2);
    const auto observation_ac=Phase3ObservationAC(problem,targets,baseline.halos,knn2);
    const std::array<int,2> halo_ids{{4,5}};
    const std::array<GaussianModel3D,2> halo_truth{{{.5,.8,.3},{.3,1.0,-.25}}};
    const auto target_field=Phase3TargetTruthField(fixture.map);
    std::vector<Phase3Row> rows;
    const auto add_row=[&](std::string name,std::vector<GaussianModel3D> profiles,
        std::vector<std::size_t> donor_counts) {
        rows.push_back({std::move(name),std::move(profiles),std::move(donor_counts)});
    };
    std::vector<GaussianModel3D> baseline_profiles;
    for (const auto & halo:baseline.halos)
        baseline_profiles.emplace_back(halo.amplitude,halo.width,halo.charge);
    add_row("group-key-global-fallback",baseline_profiles,{3,3});
    add_row("global-target-median",{global,global},{3,3});
    add_row("geometry-knn-1",knn1,{1,1});
    add_row("geometry-knn-2",knn2,{2,2});
    if (observation_ac) add_row("observation-ac-knn-2",*observation_ac,{2,2});

    std::cout<<"PHASE3_CSV_BEGIN\n";
    std::cout<<"policy,donors_per_halo,halo_a_rmse,halo_b_rmse,halo_c_rmse,halo_signal_rmse,"
        <<"halo_signal_abs_median,halo_signal_abs_max,peeling_rmse,peeling_coverage,tail_ratio_mae,"
        <<"local_a_rmse,local_b_rmse,local_available,mean_effective_samples,mean_alpha_r,"
        <<"group_alpha,group_mdpde_a,group_mdpde_b,group_truth_rmse,group_prior_sd_a,group_prior_sd_b,"
        <<"group_outliers,reconstruction_only_nodes,cache_peak,estimated_workspace_bytes,"
        <<"observation_design_rows,observation_design_columns,observation_dense_design_bytes\n";

    for (const auto & row:rows)
    {
        double halo_a_sq{},halo_b_sq{},halo_c_sq{},halo_signal_sq{};
        std::size_t halo_signal_count{};
        std::vector<double> halo_abs_errors;
        for (std::size_t h=0;h<row.halos.size();++h)
        {
            const auto & point=row.halos[h];
            halo_a_sq+=std::pow(point.GetAmplitude()-halo_truth[h].GetAmplitude(),2);
            halo_b_sq+=std::pow(point.GetWidth()-halo_truth[h].GetWidth(),2);
            halo_c_sq+=std::pow(point.GetOffset()-halo_truth[h].GetOffset(),2);
            const auto center=fixture.model.FindAtomPtr(halo_ids[h])->GetPosition();
            for (std::size_t grid=0;grid<fixture.map.GetMapValueArraySize();++grid)
            {
                const auto position=fixture.map.GetGridPosition(grid);
                double square{};
                for (std::size_t axis=0;axis<3;++axis)
                {
                    const auto delta=position[axis]-center[axis];
                    square+=delta*delta;
                }
                if (square>6.25) continue;
                const auto estimated=core::joint_component::EvaluateKernel(square,point.GetWidth(),2.5);
                const auto truth=core::joint_component::EvaluateKernel(square,halo_truth[h].GetWidth(),2.5);
                const double error=point.GetAmplitude()*estimated.gaussian+point.GetOffset()*estimated.charge-
                    halo_truth[h].GetAmplitude()*truth.gaussian-halo_truth[h].GetOffset()*truth.charge;
                halo_signal_sq+=error*error;
                ++halo_signal_count;
                halo_abs_errors.push_back(std::abs(error));
            }
        }
        std::sort(halo_abs_errors.begin(),halo_abs_errors.end());
        const std::optional<double> halo_signal_median=halo_abs_errors.empty() ? std::nullopt :
            std::optional<double>{halo_abs_errors[halo_abs_errors.size()/2]};
        const std::optional<double> halo_signal_max=halo_abs_errors.empty() ? std::nullopt :
            std::optional<double>{halo_abs_errors.back()};

        auto candidate=*result;
        auto build=baseline_build;
        for (std::size_t h=0;h<row.halos.size();++h)
        {
            auto & saved=candidate.halo_median_reconstruction->halos[h];
            auto & rebuilt=build.snapshot.halos[h];
            saved.amplitude=rebuilt.amplitude=row.halos[h].GetAmplitude();
            saved.width=rebuilt.width=row.halos[h].GetWidth();
            saved.charge=rebuilt.charge=row.halos[h].GetOffset();
            const auto count=row.donor_counts[h];
            saved.donor_count=rebuilt.donor_count=count;
            saved.source=rebuilt.source=count==baseline.valid_target_count ?
                HaloMedianSource::GlobalTargetMedian : HaloMedianSource::GroupKeyMedian;
        }
        core::detail::JointPostprocessingMetrics metrics;
        core::detail::JointPostprocessingContext context(fixture.map,fixture.model,problem,
            &candidate,&metrics,&build);
        auto peelings=context.BuildPostFitPeelingSamples(
            std::span<const std::size_t>(candidate.selection_domain->target_indices));
        double peeling_sq{};
        std::size_t peeling_count{},peeling_total{};
        for (const auto & target:targets)
        {
            const auto peeling=peelings.find(target.id);
            ASSERT_NE(peeling,peelings.end());
            const auto raw=AtomLocalPotentialView::For(*fixture.model.FindAtomPtr(target.id)).GetRawSamplingEntriesView();
            for (std::size_t sample=0;sample<raw.size();++sample)
            {
                ++peeling_total;
                if (!peeling->second.samples[sample].response) continue;
                const auto truth=ReferenceTricubicCharge(fixture.map,target_field,raw[sample].point.position);
                const auto error=*peeling->second.samples[sample].response-truth;
                peeling_sq+=error*error;
                ++peeling_count;
            }
        }
        const std::optional<double> peeling_rmse=peeling_count ?
            std::optional<double>{std::sqrt(peeling_sq/static_cast<double>(peeling_count))} : std::nullopt;
        const double peeling_coverage=peeling_total ? static_cast<double>(peeling_count)/
            static_cast<double>(peeling_total) : 0.0;
        std::vector<std::pair<int,PostFitPeelingResult>> peeling_updates;
        for (const auto & [id,peeling]:peelings) peeling_updates.emplace_back(id,peeling);
        fixture.model.EditAnalysis().ApplyPostFitPeelingsBatch(std::move(peeling_updates));

        double tail_ratio_abs_error{};
        std::size_t tail_count{};
        for (const auto & target:targets)
        {
            const auto view=AtomLocalPotentialView::For(*fixture.model.FindAtomPtr(target.id));
            const auto assessment=view.AssessLocalFittingPeelingRatio(1.0,2.0);
            if (!assessment.ratio) continue;
            const auto raw=view.GetRawSamplingEntriesView();
            double raw_sum{},true_target_sum{};
            for (const auto & sample:raw)
                if (sample.point.distance>=1.0 && sample.point.distance<=2.0)
                {
                    raw_sum+=sample.response;
                    true_target_sum+=ReferenceTricubicCharge(fixture.map,target_field,sample.point.position);
                }
            if (raw_sum==0.0) continue;
            const double truth_ratio=1.0-true_target_sum/raw_sum;
            tail_ratio_abs_error+=std::abs(*assessment.ratio-truth_ratio);
            ++tail_count;
        }

        std::vector<std::pair<int,PostJointLocalResult>> local_updates;
        for (const auto & target:targets)
        {
            const auto found=std::find_if(peelings.begin(),peelings.end(),[&](const auto & item) {
                return item.first==target.id;
            });
            if (found==peelings.end()) continue;
            local_updates.emplace_back(target.id,core::detail::FitPostJointLocalTargetForResearch(
                context,fixture.model,target.id,found->second,options,metrics));
        }
        fixture.model.EditAnalysis().ApplyPostJointLocalResultsBatch(std::move(local_updates));
        core::detail::RunJointGroupPotentialFitting(fixture.model,options,&metrics);

        double local_a_sq{},local_b_sq{},effective_samples{},alpha_r{};
        std::size_t local_available{};
        const std::array<double,3> truth_a{{2.0,1.6,1.2}},truth_b{{.65,.75,.85}};
        for (std::size_t index=0;index<targets.size();++index)
        {
            const auto local=AtomLocalPotentialView::For(*fixture.model.FindAtomPtr(targets[index].id))
                .GetPostJointLocalResult();
            if (!local || local->status!=EvidenceStatus::Available) continue;
            const auto & point=local->fit.mdpde.GetModel();
            local_a_sq+=std::pow(point.GetAmplitude()-truth_a[index],2);
            local_b_sq+=std::pow(point.GetWidth()-truth_b[index],2);
            effective_samples+=static_cast<double>(local->effective_sample_count);
            alpha_r+=local->fit.alpha_r;
            ++local_available;
        }
        const auto key=data_internal::GetGroupKey(fixture.model.FindAtomPtr(1));
        const auto summary=fixture.model.GetAnalysisView().GetGroupParameterSummary(key);
        std::optional<double> group_alpha,group_a,group_b,group_truth_rmse,prior_sd_a,prior_sd_b;
        std::size_t outliers{};
        if (summary)
        {
            group_alpha=fixture.model.GetAnalysisView().GetAtomAlphaG(key);
            if (summary->inference)
            {
                const auto & inference=*summary->inference;
                group_a=inference.mdpde.GetAmplitude();
                group_b=inference.mdpde.GetWidth();
                group_truth_rmse=std::sqrt((std::pow(*group_a-1.6,2)+std::pow(*group_b-.75,2))/2.0);
                prior_sd_a=inference.prior.GetStandardDeviationModel().GetAmplitude();
                prior_sd_b=inference.prior.GetStandardDeviationModel().GetWidth();
                outliers=static_cast<std::size_t>(std::count_if(inference.member_results.begin(),
                    inference.member_results.end(),[](const auto & member){return member.is_outlier;}));
            }
        }

        std::cout<<row.policy<<','<<row.donor_counts.front()<<','
            <<std::sqrt(halo_a_sq/static_cast<double>(row.halos.size()))<<','
            <<std::sqrt(halo_b_sq/static_cast<double>(row.halos.size()))<<','
            <<std::sqrt(halo_c_sq/static_cast<double>(row.halos.size()))<<','
            <<std::sqrt(halo_signal_sq/static_cast<double>(halo_signal_count))<<',';
        Phase3CsvValue(halo_signal_median); std::cout<<','; Phase3CsvValue(halo_signal_max); std::cout<<',';
        Phase3CsvValue(peeling_rmse); std::cout<<','<<peeling_coverage<<',';
        const std::optional<double> tail_ratio_mae=tail_count ?
            std::optional<double>{tail_ratio_abs_error/static_cast<double>(tail_count)} : std::nullopt;
        Phase3CsvValue(tail_ratio_mae); std::cout<<',';
        Phase3CsvValue(local_available ? std::optional<double>{std::sqrt(local_a_sq/
            static_cast<double>(local_available))} : std::nullopt);
        std::cout<<',';
        Phase3CsvValue(local_available ? std::optional<double>{std::sqrt(local_b_sq/
            static_cast<double>(local_available))} : std::nullopt);
        std::cout<<','<<local_available<<',';
        Phase3CsvValue(local_available ? std::optional<double>{effective_samples/
            static_cast<double>(local_available)} : std::nullopt);
        std::cout<<',';
        Phase3CsvValue(local_available ? std::optional<double>{alpha_r/
            static_cast<double>(local_available)} : std::nullopt);
        std::cout<<','; Phase3CsvValue(group_alpha); std::cout<<','; Phase3CsvValue(group_a); std::cout<<',';
        Phase3CsvValue(group_b); std::cout<<','; Phase3CsvValue(group_truth_rmse); std::cout<<',';
        Phase3CsvValue(prior_sd_a); std::cout<<','; Phase3CsvValue(prior_sd_b); std::cout<<','<<outliers<<','
            <<metrics.reconstruction_only_node_evaluation_count<<','
            <<metrics.reconstruction_node_cache_peak_entries<<','
            <<metrics.estimated_reconstruction_workspace_bytes<<','
            <<problem.OriginalInput().observations.size()<<','<<2*baseline.halos.size()<<','
            <<problem.OriginalInput().observations.size()*2*baseline.halos.size()*sizeof(double)<<'\n';
    }
    std::cout<<"PHASE3_OBSERVATION_AC_AVAILABLE,"<<(observation_ac ? "true" : "false")<<'\n';
}
