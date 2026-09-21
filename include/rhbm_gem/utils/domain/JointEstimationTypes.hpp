#pragma once
#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace rhbm_gem {
struct JointSelectionDomain
{
    std::string contract{"fixed-selected-voxel-closure-v1"};
    std::vector<std::size_t> target_indices;
    double observation_radius{2.5},support_radius{2.5};
    std::string contributor_policy{"all-non-hydrogen"};
    bool IsValid(std::size_t atoms) const
    {
        if(contract!="fixed-selected-voxel-closure-v1" || contributor_policy!="all-non-hydrogen" ||
            observation_radius!=2.5 || support_radius!=2.5 || target_indices.empty()) return false;
        for(std::size_t i=0;i<target_indices.size();++i)
            if(target_indices[i]>=atoms || (i && target_indices[i-1]>=target_indices[i])) return false;
        return true;
    }
};
enum class JointCheckStatus {Passed,Failed,Unavailable,NotRun};
enum class JointEvidenceScope {ComponentLocal,AssembledGlobal};
struct JointCheck
{
    std::string name;
    JointCheckStatus status{JointCheckStatus::NotRun};
    JointEvidenceScope scope{JointEvidenceScope::ComponentLocal};
    std::optional<double> value,threshold;
    std::string reason;
};
struct JointRankEvidence
{
    std::string name;
    JointEvidenceScope scope{JointEvidenceScope::ComponentLocal};
    std::size_t rank{};
    double threshold{};
    std::vector<double> singular_values;
};
struct JointInitializationAtom
{
    std::string id;
    std::array<double,3> ols{},mdpde{};
    double alpha{};
    std::size_t sample_count{};
    std::optional<int> native_status;
    std::string reason;
};
struct JointInitialization
{
    bool valid{};
    std::string reason;
    std::vector<double> b;
    std::vector<JointInitializationAtom> atoms;
    std::string data_scope{"caller-provided-widths"};
};
struct JointState
{
    std::vector<double> ac,b,log_b,width_gradient;
    // 0.5 * squared residual norm / parent ObservationScale() squared.
    double objective{};
};
struct JointComponentData
{
    std::string id,stop_reason;
    std::vector<std::size_t> atoms,rows;
    bool search_completed{};
    int profile_evaluations{},reference_evaluations{},accepted_updates{},native_status{};
    std::optional<JointState> state;
    std::vector<JointCheck> evidence;
    std::vector<JointRankEvidence> ranks;
    JointCheckStatus regular_certificate{JointCheckStatus::NotRun};
};
struct JointFitCosts
{
    double initialization_seconds{},search_seconds{},search_reference_seconds{},assessment_seconds{},assembly_seconds{};
    double construction_seconds{};
};
}
