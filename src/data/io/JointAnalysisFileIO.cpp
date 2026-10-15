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
        std::ofstream json_file, csv_file, contribution_file;
        json_file.exceptions(std::ios::failbit|std::ios::badbit);
        csv_file.exceptions(std::ios::failbit|std::ios::badbit);
        json_file.open(json_path); csv_file.open(csv_path);
        json_file << json << '\n';
        auto contribution_path=csv_path; contribution_path.replace_extension(".contributions.csv");
        const bool has_legacy_groups=result.layout && !result.layout->groups.empty();
        if(has_legacy_groups)
        {
            contribution_file.exceptions(std::ios::failbit|std::ios::badbit);
            contribution_file.open(contribution_path);
            contribution_file << "GroupID,RowID,AtomIDs,Lambda,StateAvailable\n" << std::setprecision(std::numeric_limits<double>::max_digits10);
        }
        std::vector<std::string> group_ids(result.atom_ids.size());
        if(result.layout) for(const auto & group:result.layout->groups)
        {
            const auto id="row:"+result.row_ids.at(group.row);
            std::string members;
            for(auto atom:group.atoms) {group_ids.at(atom)=id; if(!members.empty()) members+=';'; members+=result.atom_ids.at(atom);}
            std::optional<double> amplitude;
            for(const auto & component:result.components) if(component.layout && component.state)
                for(std::size_t k=0;k<component.layout->groups.size();++k)
                    if(component.layout->groups[k].row==group.row) amplitude=component.state->nuisance_amplitudes.at(k);
            contribution_file << CsvText(id) << ',' << CsvText(result.row_ids.at(group.row)) << ',' << CsvText(members) << ',';
            if(amplitude) contribution_file << *amplitude;
            contribution_file << ',' << amplitude.has_value() << '\n';
        }
        csv_file << "AtomID,ComponentID,A,B,C,StateAvailable,SearchCompleted,StopReason,RuntimeConvergence,RegularCertificate,SelectionRole,Parameterization,ContributionGroupID,TargetRuntimeConvergence,PointQuality\n"
            << std::setprecision(std::numeric_limits<double>::max_digits10);
        std::vector<std::pair<const JointAnalysisComponent *,std::size_t>> mapping(result.atom_ids.size());
        for(const auto & component:result.components)
            for(std::size_t i=0;i<component.atoms.size();++i) mapping[component.atoms[i]]={&component,i};
        const bool target_only_fixed_halo=
            result.parameterization_contract==JointTargetOnlyFixedHaloParameterizationContract;
        std::vector<const JointFixedHaloParameter *> fixed_halos(result.atom_ids.size());
        if(target_only_fixed_halo)
            for(const auto & halo:result.fixed_halos) fixed_halos.at(halo.contributor_index)=&halo;
        for(std::size_t atom=0;atom<result.atom_ids.size();++atom)
        {
            const auto [component,original_index]=mapping[atom];
            const auto * fixed_halo=fixed_halos[atom];
            auto index=original_index;
            const bool full=group_ids[atom].empty();
            if(component && component->layout && full)
                index=static_cast<std::size_t>(std::find(component->layout->full_atoms.begin(),component->layout->full_atoms.end(),atom)-component->layout->full_atoms.begin());
            const bool target=result.selection_domain && std::binary_search(result.selection_domain->target_indices.begin(),result.selection_domain->target_indices.end(),atom);
            JointCheckStatus identified=JointCheckStatus::NotRun;
            std::string point_quality="NotAssessed";
            bool point_valid=false;
            if(component && component->target_evidence)
            {
                if(component->target_evidence->contract==JointTargetPointValidationContractVersion)
                {
                    const auto & points=component->target_evidence->point_atoms;
                    const auto point=std::find_if(points.begin(),points.end(),[&](const auto & a){return a.atom==atom;});
                    if(point!=points.end())
                    {
                        point_quality=point->point_quality;
                        point_valid=point->parameter_validity==JointCheckStatus::Passed &&
                            point->component_convergence==JointCheckStatus::Passed;
                    }
                }
                else
                    for(const auto & a:component->target_evidence->atoms) if(a.atom==atom)
                    {identified=a.status; point_valid=target || a.status==JointCheckStatus::Passed;}
            }
            const bool export_point=full && (target || !component || !component->target_evidence || point_valid);
            const bool state_available=fixed_halo || (component && component->state && export_point);
            const char * parameterization=full ? "FullABC" : "ObservableContributionOnly";
            if(target_only_fixed_halo)
                parameterization=fixed_halo ? "FixedFromFirstStage" : target ?
                    state_available ? "JointEstimated" : "TargetOnlyUnestimated" : "FixedHaloUnavailable";
            const std::string stop_reason=fixed_halo ? "fixed-from-first-stage" :
                target_only_fixed_halo && !target ? "fixed-halo-unavailable" :
                component ? component->stop_reason : result.initialization.reason;
            const auto runtime_convergence=fixed_halo ? JointCheckStatus::NotRun :
                component && full ? component->runtime_convergence : JointCheckStatus::Unavailable;
            const auto regular_certificate=fixed_halo ? JointCheckStatus::NotRun :
                component ? component->regular_certificate : result.regular_certificate;
            const auto target_convergence=fixed_halo ? JointCheckStatus::NotRun :
                component ? component->target_runtime_convergence : JointCheckStatus::NotRun;
            const std::string point_quality_output=target_only_fixed_halo && !target ? "NotApplicable" :
                component && component->target_evidence &&
                component->target_evidence->contract==JointTargetPointValidationContractVersion ?
                point_quality : std::string(joint_result_io::StatusText(identified));
            csv_file << CsvText(result.atom_ids[atom]) << ',' << CsvText(component ? component->id : "") << ',';
            if(fixed_halo)
                csv_file << fixed_halo->amplitude << ',' << fixed_halo->width << ',' << fixed_halo->charge;
            else if(component && component->state && export_point)
            {
                const auto & state=*component->state;
                csv_file << state.ac[2*index] << ',' << state.b[index] << ',' << state.ac[2*index+1];
            }
            else csv_file << ",,";
            csv_file << ',' << state_available << ',' << (component && component->search_completed)
                << ',' << CsvText(stop_reason)
                << ',' << joint_result_io::StatusText(runtime_convergence)
                << ',' << joint_result_io::StatusText(regular_certificate)
                << ',' << (!result.selection_domain ? "not-recorded" :
                    target ? "target" : "halo") << ',' <<
                    parameterization
                << ',' << CsvText(group_ids[atom]) << ','
                << joint_result_io::StatusText(target_convergence) << ','
                << point_quality_output << '\n';
        }
        json_file.close(); csv_file.close();
        if(contribution_file.is_open()) contribution_file.close();
    }
    catch(const std::ios_base::failure & error)
    {
        throw std::runtime_error("Failed to export joint result to '"+json_path.string()+"' / '"+csv_path.string()+"': "+error.what());
    }
}
}
