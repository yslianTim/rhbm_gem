#include <rhbm_gem/data/io/JointAnalysisFileIO.hpp>
#include "detail/JointResultJson.hpp"
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
        std::ofstream json_file, csv_file;
        json_file.exceptions(std::ios::failbit|std::ios::badbit);
        csv_file.exceptions(std::ios::failbit|std::ios::badbit);
        json_file.open(json_path); csv_file.open(csv_path);
        json_file << json << '\n';
        csv_file << "AtomID,ComponentID,A,B,C,StateAvailable,SearchCompleted,StopReason,RuntimeConvergence,RegularCertificate\n"
            << std::setprecision(std::numeric_limits<double>::max_digits10);
        std::vector<std::pair<const JointAnalysisComponent *,std::size_t>> mapping(result.atom_ids.size());
        for(const auto & component:result.components)
            for(std::size_t i=0;i<component.atoms.size();++i) mapping[component.atoms[i]]={&component,i};
        for(std::size_t atom=0;atom<result.atom_ids.size();++atom)
        {
            const auto [component,index]=mapping[atom];
            csv_file << CsvText(result.atom_ids[atom]) << ',' << CsvText(component ? component->id : "") << ',';
            if(component && component->state)
            {
                const auto & state=*component->state;
                csv_file << state.ac[2*index] << ',' << state.b[index] << ',' << state.ac[2*index+1];
            }
            else csv_file << ",,";
            csv_file << ',' << (component && component->state.has_value()) << ',' << (component && component->search_completed)
                << ',' << CsvText(component ? component->stop_reason : result.initialization.reason)
                << ',' << joint_result_io::StatusText(component ? component->runtime_convergence : JointCheckStatus::Unavailable)
                << ',' << joint_result_io::StatusText(component ? component->regular_certificate : result.regular_certificate) << '\n';
        }
        json_file.close(); csv_file.close();
    }
    catch(const std::ios_base::failure & error)
    {
        throw std::runtime_error("Failed to export joint result to '"+json_path.string()+"' / '"+csv_path.string()+"': "+error.what());
    }
}
}
