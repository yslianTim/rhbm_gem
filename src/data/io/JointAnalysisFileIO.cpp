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
        contribution_file.exceptions(std::ios::failbit|std::ios::badbit);
        json_file.open(json_path); csv_file.open(csv_path);
        json_file << json << '\n';
        auto contribution_path=csv_path; contribution_path.replace_extension(".contributions.csv");
        contribution_file.open(contribution_path);
        contribution_file << "GroupID,RowID,AtomIDs,Lambda,StateAvailable\n" << std::setprecision(std::numeric_limits<double>::max_digits10);
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
        csv_file << "AtomID,ComponentID,A,B,C,StateAvailable,SearchCompleted,StopReason,RuntimeConvergence,RegularCertificate,SelectionRole,Parameterization,ContributionGroupID\n"
            << std::setprecision(std::numeric_limits<double>::max_digits10);
        std::vector<std::pair<const JointAnalysisComponent *,std::size_t>> mapping(result.atom_ids.size());
        for(const auto & component:result.components)
            for(std::size_t i=0;i<component.atoms.size();++i) mapping[component.atoms[i]]={&component,i};
        for(std::size_t atom=0;atom<result.atom_ids.size();++atom)
        {
            const auto [component,original_index]=mapping[atom];
            auto index=original_index;
            const bool full=group_ids[atom].empty();
            if(component && component->layout && full)
                index=static_cast<std::size_t>(std::find(component->layout->full_atoms.begin(),component->layout->full_atoms.end(),atom)-component->layout->full_atoms.begin());
            csv_file << CsvText(result.atom_ids[atom]) << ',' << CsvText(component ? component->id : "") << ',';
            if(component && component->state && full)
            {
                const auto & state=*component->state;
                csv_file << state.ac[2*index] << ',' << state.b[index] << ',' << state.ac[2*index+1];
            }
            else csv_file << ",,";
            csv_file << ',' << (component && component->state.has_value() && full) << ',' << (component && component->search_completed)
                << ',' << CsvText(component ? component->stop_reason : result.initialization.reason)
                << ',' << joint_result_io::StatusText(component && full ? component->runtime_convergence : JointCheckStatus::Unavailable)
                << ',' << joint_result_io::StatusText(component ? component->regular_certificate : result.regular_certificate)
                << ',' << (!result.selection_domain ? "not-recorded" :
                    std::binary_search(result.selection_domain->target_indices.begin(),result.selection_domain->target_indices.end(),atom) ? "target" : "halo") << ',' << (full ? "FullABC" : "ObservableContributionOnly")
                << ',' << CsvText(group_ids[atom]) << '\n';
        }
        json_file.close(); csv_file.close(); contribution_file.close();
    }
    catch(const std::ios_base::failure & error)
    {
        throw std::runtime_error("Failed to export joint result to '"+json_path.string()+"' / '"+csv_path.string()+"': "+error.what());
    }
}
}
