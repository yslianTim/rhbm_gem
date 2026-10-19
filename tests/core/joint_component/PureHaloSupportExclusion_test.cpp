#include <gtest/gtest.h>
#include <boost/json.hpp>

#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/MapSampler.hpp"
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

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
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
    medians.halos.push_back({1,"2",data_internal::GetGroupKey(halo),.5,.8,0.0,
        HaloMedianSource::GlobalTargetMedian,1,"pure-peeling-run",
        std::string(HaloMedianReconstructionContractVersion)});
    result.halo_median_reconstruction=std::move(medians);
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

TEST(PureHaloSupportExclusionTest, FitsTargetOnlyWithoutFixedHaloSnapshot)
{
    const JointProblem problem(PureInput());
    EXPECT_TRUE(problem.FixedHalos().empty());
    const auto fit=FitJointComponentsInitialized(problem,TargetInitialization());
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
        {2.25,{2.25,0,0},true},
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
    EXPECT_FALSE(samples[1].response);
    EXPECT_FALSE(samples[2].response);
    EXPECT_EQ(samples[1].reason,"outside-joint-domain");
    EXPECT_EQ(samples[2].reason,"outside-joint-domain");
    EXPECT_EQ(metrics.raw_peeling_sample_count,3u);
    EXPECT_EQ(metrics.available_peeling_sample_count,1u);
    EXPECT_EQ(metrics.unavailable_peeling_sample_count,2u);
    EXPECT_EQ(metrics.map_geometry_boundary_samples,1u);
    EXPECT_EQ(metrics.original_support_boundary_samples,1u);
    EXPECT_EQ(metrics.mapping_mismatch_samples,0u);
    EXPECT_GE(metrics.fitting_mask_only_intersection_samples,1u);
    EXPECT_GT(metrics.nonzero_stencil_node_count,0u);
    EXPECT_GT(metrics.map_geometry_boundary_nodes,0u);
    EXPECT_GT(metrics.original_support_boundary_nodes,0u);
    EXPECT_GT(metrics.fitting_mask_only_intersection_nodes,0u);
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
    ASSERT_NO_THROW(core::detail::RunJointComponentWorkflow(
        fixture.map,fixture.model,options,&metrics,true));
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
        fixture.map,fixture.model,options,&metrics,true));
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
        fixture.map,fixture.model,options,&metrics,true));
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
    ASSERT_NO_THROW(core::detail::RunJointComponentWorkflow(
        fixture.map,fixture.model,options,nullptr,true));
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
    ASSERT_NO_THROW(core::detail::RunJointComponentWorkflow(
        fixture.map,fixture.model,options,nullptr,true));
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
    ASSERT_NO_THROW(core::detail::RunJointComponentWorkflow(
        fixture.map,fixture.model,options,&metrics,true));
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
    ASSERT_NO_THROW(core::detail::RunJointComponentWorkflow(
        fixture.map,fixture.model,options,&metrics,true));

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

TEST(PureHaloSupportExclusionTest, ProductionWorkflowUsesPureExclusionDefault)
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
    ASSERT_TRUE(result->halo_median_reconstruction);
    EXPECT_EQ(result->halo_median_reconstruction->halos.size(),2u);
}

TEST(PureHaloSupportExclusionTest, PublicEstimatorEntrypointsDeclareTheirPolicies)
{
    auto production_fixture=WorkflowFixture();
    const auto production=EstimateJointComponents(production_fixture.map,production_fixture.model);
    ASSERT_TRUE(production.problem);
    ASSERT_TRUE(production.problem->Input().selection_domain);
    EXPECT_EQ(production.problem->Input().selection_domain->contributor_policy,
        PureHaloSupportExclusionContributorPolicy);
    EXPECT_EQ(production.problem->Input().atom_ids.size(),1u);

    auto legacy_fixture=WorkflowFixture();
    const auto legacy=EstimateJointComponentsFixedHalo(legacy_fixture.map,legacy_fixture.model);
    ASSERT_TRUE(legacy.problem);
    ASSERT_TRUE(legacy.problem->Input().selection_domain);
    EXPECT_EQ(legacy.problem->Input().selection_domain->contributor_policy,FixedHaloContributorPolicy);
    EXPECT_EQ(legacy.problem->Input().atom_ids.size(),3u);
}
