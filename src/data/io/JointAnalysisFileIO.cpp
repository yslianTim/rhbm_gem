#include <rhbm_gem/data/io/JointAnalysisFileIO.hpp>
#include "detail/JointResultJson.hpp"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace rhbm_gem {
namespace {
std::string CsvText(const std::string & text)
{
    std::string result="\"";
    for(char c:text) {if(c=='"') result+='"'; result+=c;}
    return result+'"';
}
}
void WriteJointAnalysisResult(const JointAnalysisResult & result,
    const std::filesystem::path & json_path, const std::filesystem::path & csv_path)
{
    const auto json=joint_result_io::Encode(result);
    try
    {
        std::ofstream json_file,csv_file;
        json_file.exceptions(std::ios::failbit|std::ios::badbit);
        csv_file.exceptions(std::ios::failbit|std::ios::badbit);
        json_file.open(json_path); csv_file.open(csv_path);
        json_file << json << '\n';
        csv_file << "AtomID,ComponentID,A,B,C,StateAvailable,SearchCompleted,StopReason,RuntimeConvergence,RegularCertificate,SelectionRole,Parameterization,ContributionGroupID,TargetRuntimeConvergence,PointQuality\n"
            << std::setprecision(std::numeric_limits<double>::max_digits10);

        std::vector<std::pair<const JointAnalysisComponent *,std::size_t>> mapping(result.atom_ids.size());
        for(const auto & component:result.components)
            for(std::size_t i=0;i<component.atoms.size();++i) mapping[component.atoms[i]]={&component,i};
        std::vector<const JointHaloMedianReconstruction *> median_halos(result.atom_ids.size());
        if(result.halo_median_reconstruction)
            for(const auto & halo:result.halo_median_reconstruction->halos)
                median_halos.at(halo.halo_contributor_index)=&halo;

        for(std::size_t atom=0;atom<result.atom_ids.size();++atom)
        {
            const auto [component,component_index]=mapping[atom];
            const auto * median_halo=median_halos[atom];
            const bool target=std::binary_search(result.selection_domain->target_indices.begin(),
                result.selection_domain->target_indices.end(),atom);
            std::size_t state_index=component_index;
            bool state_mapping_valid=component!=nullptr;
            if(component && component->layout)
            {
                const auto found=std::find(component->layout->full_atoms.begin(),component->layout->full_atoms.end(),atom);
                state_mapping_valid=found!=component->layout->full_atoms.end();
                if(state_mapping_valid)
                    state_index=static_cast<std::size_t>(found-component->layout->full_atoms.begin());
            }
            std::string point_quality="NotAssessed";
            if(component && component->target_evidence)
            {
                const auto & points=component->target_evidence->point_atoms;
                const auto point=std::find_if(points.begin(),points.end(),[&](const auto & value) {
                    return value.atom==atom;
                });
                if(point!=points.end())
                    point_quality=point->point_quality;
            }
            const bool export_point=target && component && component->state && state_mapping_valid;
            const bool state_available=median_halo || export_point;
            const char * parameterization=median_halo ? "HaloMedianReconstructed" : target ?
                state_available ? "JointEstimated" : "TargetOnlyUnestimated" :
                "HaloReconstructionUnavailable";
            const std::string stop_reason=median_halo ? "geometry-local-target-median" : !target ?
                "halo-reconstruction-unavailable" : component ? component->stop_reason : result.initialization.reason;
            const auto runtime_convergence=median_halo ? JointCheckStatus::NotRun :
                component && export_point ? component->runtime_convergence : JointCheckStatus::Unavailable;
            const auto regular_certificate=median_halo ? JointCheckStatus::NotRun :
                component ? component->regular_certificate : result.regular_certificate;
            const auto target_convergence=median_halo ? JointCheckStatus::NotRun :
                component && target ? component->target_runtime_convergence : JointCheckStatus::NotRun;

            csv_file << CsvText(result.atom_ids[atom]) << ',' << CsvText(component ? component->id : "") << ',';
            if(median_halo)
                csv_file << median_halo->amplitude << ',' << median_halo->width << ',' << median_halo->charge;
            else if(export_point)
            {
                const auto & state=*component->state;
                csv_file << state.ac[2*state_index] << ',' << state.b[state_index] << ',' << state.ac[2*state_index+1];
            }
            else csv_file << ",,";
            csv_file << ',' << state_available << ',' << (component && component->search_completed)
                << ',' << CsvText(stop_reason)
                << ',' << joint_result_io::StatusText(runtime_convergence)
                << ',' << joint_result_io::StatusText(regular_certificate)
                << ',' << (target ? "target" : "halo") << ',' << parameterization
                << ',' << CsvText("") << ',' << joint_result_io::StatusText(target_convergence)
                << ',' << (target ? point_quality : "NotApplicable") << '\n';
        }
        json_file.close(); csv_file.close();
    }
    catch(const std::ios_base::failure & error)
    {
        throw std::runtime_error("Failed to export joint result to '"+json_path.string()+"' / '"+csv_path.string()+"': "+error.what());
    }
}
}
