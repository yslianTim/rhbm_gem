#include "core/detail/joint_component/Problem.hpp"
#include "support/JointReproducibilityFingerprint.hpp"

#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <boost/json.hpp>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {
namespace core=rhbm_gem::core;
namespace json=boost::json;

void Write(const std::filesystem::path & path,const json::value & value)
{
    std::ofstream output(path);
    output.exceptions(std::ios::failbit|std::ios::badbit);
    output<<json::serialize(value)<<'\n';
}

const char * Status(rhbm_gem::JointCheckStatus status)
{
    using Status=rhbm_gem::JointCheckStatus;
    switch(status)
    {
    case Status::Passed:return "Passed";
    case Status::Failed:return "Failed";
    case Status::Unavailable:return "Unavailable";
    case Status::NotRun:return "NotRun";
    }
    return "NotRun";
}

void Run(const std::filesystem::path & model_path,const std::filesystem::path & map_path,
    const std::filesystem::path & output)
{
    if(std::filesystem::exists(output)) throw std::runtime_error("Use a fresh reproducibility output directory.");
    std::filesystem::create_directories(output);
    auto model=rhbm_gem::ReadModel(model_path); auto map=rhbm_gem::ReadMap(map_path);
    const auto normalization=map->MapValueArrayNormalization();
    if(!normalization) throw std::runtime_error("Map normalization unavailable.");
    model->SelectAtoms([](const rhbm_gem::AtomObject & atom) {
        return atom.GetElement()!=Element::HYDROGEN && atom.GetChainID()=="A";
    });
    model->EditAnalysis().InitializeFromSelection();
    const auto started=std::chrono::steady_clock::now();
    auto fit=core::EstimateJointComponents(*map,*model);
    if(!fit.problem) throw std::runtime_error("Production fit did not retain its problem snapshot.");
    const auto & data=core::JointProblemAccess::Get(*fit.problem);
    const auto fingerprints=rhbm_gem::joint_reproducibility::WriteFingerprint(output,
        model_path,map_path,*normalization,*map,*model,*fit.problem,fit.initialization.b,
        {"production-no-instrumentation",false});
    json::array components;
    for(const auto & component:fit.components)
        components.push_back(json::object{
            {"component_id",component.id},
            {"search_sweeps",component.search_sweeps},
            {"profile_evaluations",component.profile_evaluations},
            {"accepted_updates",component.accepted_updates},
            {"search_completed",component.search_completed},
            {"stop_reason",component.stop_reason},
            {"runtime_convergence",Status(component.RuntimeConvergence())},
            {"candidate_objective",component.candidate_diagnostic && component.candidate_diagnostic->candidate_objective ?
                json::value(*component.candidate_diagnostic->candidate_objective) : json::value(nullptr)},
            {"candidate_ac_kkt",component.candidate_diagnostic && component.candidate_diagnostic->candidate_ac_kkt ?
                json::value(*component.candidate_diagnostic->candidate_ac_kkt) : json::value(nullptr)},
            {"candidate_width_gradient_inf_norm",component.candidate_diagnostic && component.candidate_diagnostic->candidate_width_gradient_inf_norm ?
                json::value(*component.candidate_diagnostic->candidate_width_gradient_inf_norm) : json::value(nullptr)}});
    Write(output/"summary.json",json::object{
        {"execution_role","production-no-instrumentation"},
        {"instrumentation_enabled",false},
        {"target_count",data.input->selection_domain->target_indices.size()},
        {"observation_count",data.input->observations.size()},
        {"formal_runtime_convergence",Status(fit.RuntimeConvergence())},
        {"search_completed",fit.search_completed},
        {"objective",fit.objective ? json::value(*fit.objective) : json::value(nullptr)},
        {"initialization_valid",fit.initialization.valid},
        {"initialization_reason",fit.initialization.reason},
        {"components",std::move(components)},
        {"fingerprints",fingerprints},
        {"wall_seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()}});
}
}

int main(int argc,char ** argv)
{
    try
    {
        if(argc!=4) throw std::invalid_argument("Usage: joint_width_repro_no_instrumentation MODEL MAP OUTPUT_DIR");
        Run(argv[1],argv[2],argv[3]);
        return 0;
    }
    catch(const std::exception & error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
