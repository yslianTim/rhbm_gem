#include "JointReproducibilityFingerprint.hpp"
#include "core/detail/joint_component/StructuralPartition.hpp"

#include "SimulationBuildInfo.hpp"
#include "utils/domain/FileFingerprint.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <boost/json.hpp>

#include <bit>
#include <array>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace rhbm_gem::joint_reproducibility {
namespace {
namespace json=boost::json;
namespace joint=core::joint_component;

void Token(std::ostream & output,const std::string & value)
{
    output<<value.size()<<':'<<value<<'\n';
}

void Bits(std::ostream & output,double value)
{
    output<<std::hex<<std::setw(16)<<std::setfill('0')
        <<std::bit_cast<std::uint64_t>(value)<<std::dec<<std::setfill(' ')<<'\n';
}

void Indices(std::ostream & output,const std::vector<Eigen::Index> & values)
{
    output<<values.size()<<'\n';
    for(const auto value:values) output<<value<<'\n';
}

void Strings(std::ostream & output,const std::vector<std::string> & values)
{
    output<<values.size()<<'\n';
    for(const auto & value:values) Token(output,value);
}

json::array Strings(const std::vector<std::string> & values)
{
    json::array out;
    for(const auto & value:values) out.push_back(json::value(value));
    return out;
}

std::string WriteMaterial(const std::filesystem::path & path,
    const std::filesystem::path & model_path,const std::filesystem::path & map_path,
    double normalization_divisor,const MapObject & map,const ModelObject & model,
    const core::JointProblem & problem,const std::vector<double> & initial_widths,
    const FingerprintOptions & options)
{
    const auto & data=core::JointProblemAccess::Get(problem);
    if(initial_widths.size()!=data.input->atom_ids.size())
        throw std::invalid_argument("Fingerprint initialization dimensions do not match contributors.");
    std::ofstream input(path/"fingerprint-input.material");
    std::ofstream problem_output(path/"fingerprint-problem.material");
    std::ofstream initialization(path/"fingerprint-initialization.material");
    std::ofstream partition(path/"fingerprint-partition.material");
    std::ofstream policy(path/"fingerprint-policy.material");
    for(auto * stream:{&input,&problem_output,&initialization,&partition,&policy})
        stream->exceptions(std::ios::failbit|std::ios::badbit);

    Token(input,FileSha256(model_path));
    Token(input,FileSha256(map_path));
    Token(input,"MapValueArrayNormalization:standard-deviation-divisor");
    Bits(input,normalization_divisor);
    input<<map.GetMapValueArraySize()<<'\n';
    for(const auto value:map.GetGridSize()) input<<value<<'\n';
    for(const auto value:map.GetGridSpacing()) Bits(input,value);
    for(const auto value:map.GetOrigin()) Bits(input,value);
    std::vector<std::string> selected_ids;
    for(const auto * atom:model.GetSelectedAtoms()) selected_ids.push_back(std::to_string(atom->GetSerialID()));
    Strings(input,selected_ids);

    const auto & joint_input=*data.input;
    Strings(problem_output,joint_input.atom_ids);
    Strings(problem_output,joint_input.row_ids);
    problem_output<<joint_input.observations.size()<<'\n';
    for(const auto value:joint_input.observations) Bits(problem_output,value);
    problem_output<<joint_input.support.size()<<'\n';
    for(const auto & supports:joint_input.support)
    {
        problem_output<<supports.size()<<'\n';
        for(const auto & support:supports)
        {
            problem_output<<support.row<<'\n';
            Bits(problem_output,support.squared_distance);
        }
    }
    if(!joint_input.selection_domain) throw std::invalid_argument("Joint problem has no selection domain.");
    std::vector<Eigen::Index> target_indices;
    target_indices.reserve(joint_input.selection_domain->target_indices.size());
    for(const auto index:joint_input.selection_domain->target_indices)
        target_indices.push_back(static_cast<Eigen::Index>(index));
    Indices(problem_output,target_indices);

    initialization<<initial_widths.size()<<'\n';
    for(std::size_t index=0;index<initial_widths.size();++index)
    {
        Token(initialization,joint_input.atom_ids.at(index));
        Bits(initialization,initial_widths.at(index));
    }

    partition<<data.partition.components.size()<<'\n';
    for(const auto & component:data.partition.components)
    {
        Token(partition,component.id);
        Indices(partition,component.atoms);
        Indices(partition,component.rows);
        const auto prepared=joint::MakeComponentProblemView(data.layout,component);
        const auto local_ids=data.context.atom_ids.Select(prepared.local_to_parent_atoms);
        joint::StructuralTopology topology(prepared.domain,local_ids,prepared.local_layout);
        const auto blocks=topology.BuildCores(12,false);
        partition<<blocks.cores.size()<<'\n';
        for(const auto & block:blocks.cores)
        {
            Token(partition,block.id);
            Indices(partition,block.atoms);
            Indices(partition,block.affected_rows);
        }
    }

    Token(policy,options.execution_role);
    policy<<options.instrumentation_enabled<<'\n';
    policy<<12<<'\n'<<200<<'\n'<<"forward\n"<<"fixed-neighbor\n"<<"profile-lm\n";
    Bits(policy,joint::PracticalJointConvergenceCriteria.global_ac_kkt);
    Bits(policy,joint::PracticalJointConvergenceCriteria.global_width_gradient);
    Bits(policy,joint::PracticalJointConvergenceCriteria.confirmed_log_width_change);
    policy<<data.context.profile_budget<<'\n'<<data.context.update_budget<<'\n';
    Bits(policy,data.context.scale);
    Bits(policy,data.context.linear.rank_relative);
    Bits(policy,data.context.linear.release_factor);
    Bits(policy,data.context.linear.release_response_norm);
    Token(policy,RHBM_GEM_SIMULATION_VERSION);
    Token(policy,RHBM_GEM_SIMULATION_SOURCE_SHA256);
    Token(policy,RHBM_GEM_SIMULATION_CONFIG_SHA256);
    Token(policy,RHBM_GEM_SIMULATION_BUILD_SHA256);
#ifdef __clang_version__
    Token(policy,__clang_version__);
#elif defined(__VERSION__)
    Token(policy,__VERSION__);
#else
    Token(policy,"unknown-compiler");
#endif
    input.close(); problem_output.close(); initialization.close(); partition.close(); policy.close();
    return {};
}

} // namespace

boost::json::object WriteFingerprint(
    const std::filesystem::path & output,const std::filesystem::path & model_path,
    const std::filesystem::path & map_path,double normalization_divisor,const MapObject & map,
    const ModelObject & model,const core::JointProblem & problem,const std::vector<double> & initial_widths,
    const FingerprintOptions & options)
{
    std::filesystem::create_directories(output);
    WriteMaterial(output,model_path,map_path,normalization_divisor,map,model,problem,initial_widths,options);
    const auto & data=core::JointProblemAccess::Get(problem);
    const auto target_count=data.input->selection_domain ? data.input->selection_domain->target_indices.size() : 0;
    return {
        {"schema","joint-reproducibility-fingerprint-v1"},
        {"input_sha256",json::object{
            {"model",FileSha256(model_path)},
            {"map",FileSha256(map_path)},
            {"input_material",FileSha256(output/"fingerprint-input.material")}}},
        {"map_normalization",json::object{
            {"operation","MapValueArrayNormalization"},
            {"divisor",normalization_divisor},
            {"normalized_sd",map.GetMapValueSD()}}},
        {"selection",json::object{
            {"selected_atom_ids",[&] {
                std::vector<std::string> ids;
                for(const auto * atom:model.GetSelectedAtoms()) ids.push_back(std::to_string(atom->GetSerialID()));
                return Strings(ids);
            }()},
            {"target_ids",[&] {
                std::vector<std::string> ids;
                for(const auto index:data.input->selection_domain->target_indices)
                    ids.push_back(data.input->atom_ids.at(index));
                return Strings(ids);
            }()}}},
        {"problem",json::object{
            {"rows",data.input->observations.size()},
            {"contributors",data.input->atom_ids.size()},
            {"targets",target_count},
            {"material_sha256",FileSha256(output/"fingerprint-problem.material")}}},
        {"initialization",json::object{
            {"width_count",initial_widths.size()},
            {"material_sha256",FileSha256(output/"fingerprint-initialization.material")}}},
        {"partition",json::object{
            {"component_count",data.partition.components.size()},
            {"core_atoms",12},
            {"block_order","forward"},
            {"material_sha256",FileSha256(output/"fingerprint-partition.material")}}},
        {"policy",json::object{
            {"execution_role",options.execution_role},
            {"instrumentation_enabled",options.instrumentation_enabled},
            {"material_sha256",FileSha256(output/"fingerprint-policy.material")}}},
        {"build",json::object{
            {"version",RHBM_GEM_SIMULATION_VERSION},
            {"source_sha256",RHBM_GEM_SIMULATION_SOURCE_SHA256},
            {"configuration_sha256",RHBM_GEM_SIMULATION_CONFIG_SHA256},
            {"build_sha256",RHBM_GEM_SIMULATION_BUILD_SHA256}}}
    };
}

} // namespace rhbm_gem::joint_reproducibility
