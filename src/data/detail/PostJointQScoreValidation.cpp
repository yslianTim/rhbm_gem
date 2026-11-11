#include "PostJointQScoreValidation.hpp"
#include "JointStageAdapter.hpp"
#include "ModelAnalysisData.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <algorithm>
#include <cmath>

namespace rhbm_gem::data_internal {
namespace {
bool SameSource(const EstimateSource & a,const EstimateSource & b)
{
    return a.method==b.method && a.role==b.role && a.atom_id==b.atom_id &&
        a.component_id==b.component_id && a.run_id==b.run_id &&
        a.fixed_from_initial_stage==b.fixed_from_initial_stage;
}
}
PostJointQScoreValidation::PostJointQScoreValidation(const ModelObject & model):m_model(model)
{
    const auto & result=ModelAnalysisData::Of(model).joint_result;
    m_snapshot=result && IsQualifiedJointAnalysisResult(*result) ? &*result : nullptr;
    if (!m_snapshot) return;
    m_expected=BuildJointStageEstimates(*m_snapshot,"");
    m_contributors.insert(m_snapshot->atom_ids.begin(),m_snapshot->atom_ids.end());
    if (m_snapshot->target_evidence)
        for (const auto & point:m_snapshot->target_evidence->point_atoms)
            if (point.component_convergence==JointCheckStatus::Passed &&
                point.parameter_validity==JointCheckStatus::Passed && point.local_conditioning=="Adequate")
                m_donors.insert(m_snapshot->atom_ids.at(point.atom));
    for (const auto & id:m_donors)
    {
        const auto & endpoint=m_expected.at(std::stoi(id));
        if (!endpoint.point) continue;
        const auto & point=*endpoint.point;
        m_donor_records.push_back({0,id,model.FindAtomPtr(std::stoi(id))->GetPositionRef(),
            {point.GetAmplitude(),point.GetWidth(),point.GetOffset()}});
    }
    m_donor_index.emplace(m_donor_records);
}
std::string PostJointQScoreValidation::Validate(const PostJointQScoreResult & q,
    const AtomObject & atom,const LocalStageEstimate & stage) const
{
    const auto found=m_expected.find(atom.GetSerialID());
    if (!m_snapshot || found==m_expected.end() || atom.GetElement()==Element::HYDROGEN ||
        q.contract!=PostJointQScoreContract || !SameSource(q.source,stage.source) ||
        q.source.method!=EstimateMethod::JointComponents || q.source.role!=FittingRole::Target ||
        q.source.atom_id!=std::to_string(atom.GetSerialID()) || q.source.run_id.empty() ||
        q.source.component_id.empty() || q.source.fixed_from_initial_stage)
        return "QScore source/contract mismatch.";
    const auto & expected=found->second;
    if (stage.source.component_id!=expected.source.component_id || stage.source.role!=expected.source.role ||
        stage.convergence!=expected.convergence || stage.point.has_value()!=expected.point.has_value() ||
        (stage.point && stage.point->ToVector()!=expected.point->ToVector()))
        return "QScore endpoint differs from Joint snapshot.";
    if (!std::isfinite(q.reference_height) || !std::isfinite(q.reference_offset) ||
        q.reference_height!=m_model.GetReferenceHeight() || q.reference_offset!=m_model.GetReferenceOffset())
        return "QScore standard reference mismatch.";
    if ((q.status!=EvidenceStatus::Available && q.status!=EvidenceStatus::Unavailable) ||
        (q.status==EvidenceStatus::Available)!=q.scores.has_value() || q.scores.has_value()!=q.reason.empty())
        return "QScore availability/reason mismatch.";
    if (q.sample_count>2128 || (q.sample_count!=0 && q.sample_count!=8 && q.sample_count<16) ||
        (q.scores && q.sample_count<=8) || (q.sample_count==0 && q.contributor_count!=0) || q.contributor_count>=m_model.GetNumberOfAtom() || q.halo_count>q.contributor_count ||
        q.extra_halos.size()>q.halo_count || q.missing_contributor_count>q.contributor_count ||
        q.evaluated_contributions>q.sample_count*q.contributor_count ||
        q.evaluated_contributions<q.contributor_count-q.missing_contributor_count ||
        (q.scores && q.missing_contributor_count))
        return "QScore sample/contributor coverage mismatch.";
    if (q.scores)
    {
        if (!stage.point || stage.convergence!=JointCheckStatus::Passed)
            return "Available QScore requires qualified endpoint.";
        for (const double value:{q.scores->standard,q.scores->joint})
            if (!std::isfinite(value) || value < -1.0-1.0e-12 || value > 1.0+1.0e-12)
                return "Invalid QScore correlation.";
    }
    if (q.halo_count>q.extra_halos.size())
    {
        const auto & halo=m_snapshot->halo_median_reconstruction;
        if (!halo)
        {
            if (q.scores || !q.missing_contributor_count)
                return "QScore missing saved Halo provenance.";
        }
        else if (halo->joint_run_id!=q.source.run_id || halo->contract!=HaloMedianReconstructionContractVersion ||
            halo->donor_policy!=GeometryLocalTargetMedianDonorPolicy)
            return "QScore saved Halo provenance mismatch.";
    }
    int previous{};
    for (const auto & halo:q.extra_halos)
    {
        if (halo.atom_id<=previous || m_contributors.contains(std::to_string(halo.atom_id)) ||
            halo.donor_atom_ids.size()!=std::min(GeometryLocalTargetMedianNeighborCount,m_donors.size()) ||
            halo.donor_atom_ids.empty()) return "Invalid QScore extra Halo coverage.";
        previous=halo.atom_id;
        const AtomObject * object=nullptr;
        try { object=m_model.FindAtomPtr(halo.atom_id); }
        catch (const std::out_of_range &) { return "Unknown QScore extra Halo."; }
        if (!object || object->GetElement()==Element::HYDROGEN) return "Invalid QScore extra Halo.";
        auto donor_ids_found=m_extra_donor_ids.find(halo.atom_id);
        if (donor_ids_found==m_extra_donor_ids.end())
        {
            std::vector<std::string> nearest;
            for (const auto * donor:m_donor_index->Query(object->GetPositionRef())) nearest.push_back(donor->atom_id);
            donor_ids_found=m_extra_donor_ids.emplace(halo.atom_id,std::move(nearest)).first;
        }
        if (halo.donor_atom_ids!=donor_ids_found->second)
            return "QScore Halo donors differ from qualified geometry-local selection.";
    }
    return {};
}
}
