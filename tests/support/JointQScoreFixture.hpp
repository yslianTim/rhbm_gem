#pragma once

#include "core/estimator/joint_component/PostJointQScore.hpp"
#include "support/JointComponentChecks.hpp"
#include <rhbm_gem/core/QScoreHelper.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>

namespace qscore_test {
using namespace rhbm_gem;
using namespace rhbm_gem::core;
using namespace rhbm_gem::core::detail;
struct Atom
{
    std::array<double,3> position{};
    GaussianModel3D point{2,.6,0};
    bool target{true}, hydrogen{};
};
struct Fixture
{
    std::unique_ptr<ModelObject> model;
    std::unique_ptr<MapObject> map;
    std::unique_ptr<JointProblem> problem;
    JointAnalysisResult snapshot;
    HaloMedianBuildResult halos;
    JointPostprocessingMetrics metrics;
    std::vector<std::pair<int,PostJointQScoreResult>> results;
    void Calculate()
    {
        metrics={};
        results=BuildPostJointQScores(*map,*model,*problem,snapshot,&halos,metrics);
    }
};
inline Fixture Make(std::vector<Atom> specs={Atom{}},bool coarse_x=false,bool boundary=false)
{
    Fixture f;
    std::vector<std::unique_ptr<AtomObject>> atoms;
    for (std::size_t i=0;i<specs.size();++i)
    {
        auto atom=std::make_unique<AtomObject>();
        atom->SetSerialID(static_cast<int>(i+1)); atom->SetPosition(specs[i].position);
        atom->SetElement(specs[i].hydrogen ? Element::HYDROGEN : Element::CARBON);
        atom->SetChainID("A"); atom->SetSequenceID(static_cast<int>(i+1));
        atoms.push_back(std::move(atom));
    }
    f.model=std::make_unique<ModelObject>(std::move(atoms));
    f.model->SelectAtoms([&](const auto & atom){return specs[static_cast<std::size_t>(atom.GetSerialID()-1)].target;});
    const std::array<int,3> dimensions=coarse_x ? std::array{5,25,25}:std::array{45,37,37};
    const std::array<double,3> spacing=coarse_x ? std::array{3.,.4,.4}:std::array{.25,.25,.25};
    std::array<double,3> origin=coarse_x ? std::array{-6.,-4.8,-4.8}:std::array{-5.5,-4.5,-4.5};
    if (boundary) origin[0]=0;
    const auto size=static_cast<std::size_t>(dimensions[0]*dimensions[1]*dimensions[2]);
    auto values=std::make_unique<double[]>(size);
    for (int z=0;z<dimensions[2];++z) for (int y=0;y<dimensions[1];++y) for (int x=0;x<dimensions[0];++x)
    {
        const std::array position{origin[0]+spacing[0]*x,origin[1]+spacing[1]*y,origin[2]+spacing[2]*z};
        double value{};
        for (const auto & atom:specs) if (!atom.hydrogen)
        {
            double square{};
            for (std::size_t axis=0;axis<3;++axis) {const double d=position[axis]-atom.position[axis];square+=d*d;}
            const auto kernel=joint_component::EvaluateKernel(square,atom.point.GetWidth(),2.5);
            value+=atom.point.GetAmplitude()*kernel.gaussian+atom.point.GetOffset()*kernel.charge;
        }
        values[static_cast<std::size_t>(x+dimensions[0]*(y+dimensions[1]*z))]=value;
    }
    f.map=std::make_unique<MapObject>(dimensions,spacing,origin,std::move(values));
    const auto [height,offset]=GetReferenceGaussianParameters(*f.map);
    f.model->SetReferenceHeight(height); f.model->SetReferenceOffset(offset);
    f.problem=std::make_unique<JointProblem>(BuildPureHaloSupportExclusionProblem(*f.map,*f.model));
    std::vector<double> widths;
    for (const auto & atom:specs) if (atom.target && !atom.hydrogen) widths.push_back(atom.point.GetWidth());
    f.snapshot=CaptureJointAnalysisResult(joint_component_test::matched::joint_component::FitBoundedJointComponentsForTesting(*f.problem,widths,2.0));
    if (!IsQualifiedJointAnalysisResult(f.snapshot)) throw std::runtime_error("Synthetic QScore endpoint not qualified.");
    // The bounded test adapter bypasses production donor-count calibration;
    // persist its explicit 2 A bound using the ordinary provenance contract.
    f.snapshot.evidence.push_back({"global-width-calibration",JointCheckStatus::Passed,
        JointEvidenceScope::AssembledGlobal,.5,2.,"global-initial-stage-p95-k4; test-adapter-bound"});
    f.snapshot.evidence.push_back({"qualified-initial-stage-donor-coverage",JointCheckStatus::Passed,
        JointEvidenceScope::AssembledGlobal,1.,.95,"target-only-qualified-initial-stage-widths"});
    f.snapshot.metadata.grid_size=f.map->GetGridSize(); f.snapshot.metadata.grid_spacing=f.map->GetGridSpacing();
    f.snapshot.metadata.origin=f.map->GetOrigin();
    f.halos=BuildHaloMedianReconstruction(*f.model,*f.problem,f.snapshot,"qscore-test-run");
    if (!f.halos.snapshot->halos.empty()) f.snapshot.halo_median_reconstruction=*f.halos.snapshot;
    f.model->EditAnalysis().ApplyJointResult(f.snapshot,"qscore-test-run");
    f.model->EditAnalysis().BindJointMapObservations(*f.map);
    for (const auto * atom:f.model->GetSelectedAtoms())
        if (!boundary && atom->GetElement()!=Element::HYDROGEN)
            const_cast<AtomObject *>(atom)->SetStandardQScore(CalculateQScoreForAtom(*atom,*f.map,*f.model,height,offset));
    f.Calculate();
    return f;
}
}
