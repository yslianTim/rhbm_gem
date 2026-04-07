#pragma once

#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "data/detail/GroupPotentialEntry.hpp"
#include "data/detail/LocalPotentialEntry.hpp"
#include "data/detail/ModelAnalysisData.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/KeyPacker.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>

namespace rhbm_gem {

class ModelPotentialView
{
    const ModelObject & m_model_object;

public:
    explicit ModelPotentialView(const ModelObject & model_object) :
        m_model_object(model_object) {}

    double GetAtomGausEstimateMinimum(int par_id, Element element) const
    {
        std::vector<double> gaus_estimate_list;
        gaus_estimate_list.reserve(m_model_object.GetSelectedAtomCount());
        for (const auto * atom : m_model_object.GetSelectedAtoms())
        {
            if (atom->GetElement() != element) continue;
            auto * entry{ ModelAnalysisData::FindLocalEntry(*atom) };
            gaus_estimate_list.emplace_back(entry->GetEstimateMDPDE().GetParameter(par_id));
        }
        return ArrayStats<double>::ComputeMin(gaus_estimate_list.data(), gaus_estimate_list.size());
    }

    double GetBondGausEstimateMinimum(int par_id) const
    {
        std::vector<double> gaus_estimate_list;
        gaus_estimate_list.reserve(m_model_object.GetSelectedBondCount());
        for (const auto * bond : m_model_object.GetSelectedBonds())
        {
            auto * entry{ ModelAnalysisData::FindLocalEntry(*bond) };
            gaus_estimate_list.emplace_back(entry->GetEstimateMDPDE().GetParameter(par_id));
        }
        return ArrayStats<double>::ComputeMin(gaus_estimate_list.data(), gaus_estimate_list.size());
    }

    const GroupPotentialBucket * TryGetAtomGroup(
        GroupKey group_key,
        const std::string & class_key,
        bool verbose = false) const
    {
        auto * entry{ FindAtomGroupEntry(class_key) };
        if (entry == nullptr)
        {
            if (verbose)
            {
                Logger::Log(LogLevel::Error, "Atom class key not found: " + class_key);
            }
            return nullptr;
        }
        const auto * group{ entry->FindGroup(group_key) };
        if (group == nullptr)
        {
            if (verbose)
            {
                LogMissingAtomGroup(group_key, class_key);
            }
            return nullptr;
        }
        return group;
    }

    const GroupPotentialBucket * TryGetBondGroup(
        GroupKey group_key,
        const std::string & class_key,
        bool verbose = false) const
    {
        auto * entry{ FindBondGroupEntry(class_key) };
        if (entry == nullptr)
        {
            if (verbose)
            {
                Logger::Log(LogLevel::Error, "Bond class key not found: " + class_key);
            }
            return nullptr;
        }
        const auto * group{ entry->FindGroup(group_key) };
        if (group == nullptr)
        {
            if (verbose)
            {
                LogMissingBondGroup(group_key, class_key);
            }
            return nullptr;
        }
        return group;
    }

    bool HasAtomGroup(GroupKey group_key, const std::string & class_key, bool verbose=false) const
    {
        return TryGetAtomGroup(group_key, class_key, verbose) != nullptr;
    }

    bool HasBondGroup(GroupKey group_key, const std::string & class_key, bool verbose=false) const
    {
        return TryGetBondGroup(group_key, class_key, verbose) != nullptr;
    }

    const GroupPotentialBucket & GetAtomGroup(
        GroupKey group_key,
        const std::string & class_key) const
    {
        return *RequireAtomGroup(group_key, class_key);
    }

    const GroupPotentialBucket & GetBondGroup(
        GroupKey group_key,
        const std::string & class_key) const
    {
        return *RequireBondGroup(group_key, class_key);
    }

    double GetAtomGausEstimatePrior(GroupKey group_key, const std::string & class_key, int par_id) const
    {
        return RequireAtomGroup(group_key, class_key)->prior.GetParameter(par_id);
    }

    double GetBondGausEstimatePrior(GroupKey group_key, const std::string & class_key, int par_id) const
    {
        return RequireBondGroup(group_key, class_key)->prior.GetParameter(par_id);
    }

    double GetAtomGausVariancePrior(GroupKey group_key, const std::string & class_key, int par_id) const
    {
        return BuildPriorPosterior(*RequireAtomGroup(group_key, class_key)).GetVariance(par_id);
    }

    double GetBondGausVariancePrior(GroupKey group_key, const std::string & class_key, int par_id) const
    {
        return BuildPriorPosterior(*RequireBondGroup(group_key, class_key)).GetVariance(par_id);
    }

    const std::vector<AtomObject *> & GetAtomObjectList(GroupKey group_key, const std::string & class_key) const
    {
        return RequireAtomGroup(group_key, class_key)->atom_members;
    }

    const std::vector<BondObject *> & GetBondObjectList(GroupKey group_key, const std::string & class_key) const
    {
        return RequireBondGroup(group_key, class_key)->bond_members;
    }

    std::vector<AtomObject *> GetOutlierAtomObjectList(GroupKey group_key, const std::string & class_key) const
    {
        std::vector<AtomObject *> outlier_atom_list;
        const auto & atom_list{ GetAtomObjectList(group_key, class_key) };
        outlier_atom_list.reserve(atom_list.size());
        for (auto * atom : atom_list)
        {
            const auto * annotation{
                ModelAnalysisData::RequireLocalEntry(*atom).FindAnnotation(class_key) };
            if (annotation != nullptr && annotation->is_outlier)
            {
                outlier_atom_list.emplace_back(atom);
            }
        }
        return outlier_atom_list;
    }

    std::unordered_map<int, AtomObject *> GetAtomObjectMap(GroupKey group_key, const std::string & class_key) const
    {
        std::unordered_map<int, AtomObject *> atom_object_map;
        for (auto * atom_object : GetAtomObjectList(group_key, class_key))
        {
            atom_object_map[atom_object->GetSerialID()] = atom_object;
        }
        return atom_object_map;
    }

    double GetAtomAlphaR(GroupKey group_key, const std::string & class_key) const
    {
        const auto & atom_list{ GetAtomObjectList(group_key, class_key) };
        if (atom_list.empty())
        {
            throw std::runtime_error("Atom group has no members.");
        }
        return ModelAnalysisData::RequireLocalEntry(*atom_list.front()).GetAlphaR();
    }

    double GetAtomAlphaG(GroupKey group_key, const std::string & class_key) const
    {
        return RequireAtomGroup(group_key, class_key)->alpha_g;
    }

    Residue DecodeResidueFromAtomGroupKey(GroupKey group_key, const std::string & class_key) const
    {
        if (class_key == ChemicalDataHelper::GetSimpleAtomClassKey())
        {
            Logger::Log(LogLevel::Error, "Simple class group key is not recording Residue info.");
            return Residue::UNK;
        }
        if (class_key == ChemicalDataHelper::GetComponentAtomClassKey())
        {
            auto unpack_key{ KeyPackerComponentAtomClass::Unpack(group_key) };
            return static_cast<Residue>(std::get<0>(unpack_key));
        }
        if (class_key == ChemicalDataHelper::GetStructureAtomClassKey())
        {
            auto unpack_key{ KeyPackerStructureAtomClass::Unpack(group_key) };
            return static_cast<Residue>(std::get<1>(unpack_key));
        }
        return Residue::UNK;
    }

    Residue DecodeResidueFromBondGroupKey(GroupKey group_key, const std::string & class_key) const
    {
        if (class_key == ChemicalDataHelper::GetSimpleBondClassKey())
        {
            Logger::Log(LogLevel::Error, "Simple class group key is not recording Residue info.");
            return Residue::UNK;
        }
        if (class_key == ChemicalDataHelper::GetComponentBondClassKey())
        {
            auto unpack_key{ KeyPackerComponentBondClass::Unpack(group_key) };
            return static_cast<Residue>(std::get<0>(unpack_key));
        }
        return Residue::UNK;
    }

    Residue GetResidueFromAtomGroupKey(GroupKey group_key, const std::string & class_key) const
    {
        return DecodeResidueFromAtomGroupKey(group_key, class_key);
    }

    Residue GetResidueFromBondGroupKey(GroupKey group_key, const std::string & class_key) const
    {
        return DecodeResidueFromBondGroupKey(group_key, class_key);
    }

private:
    const GroupPotentialBucket * RequireAtomGroup(GroupKey group_key, const std::string & class_key) const
    {
        auto * group{ TryGetAtomGroup(group_key, class_key, true) };
        if (group == nullptr)
        {
            throw std::runtime_error("Atom group key is not available.");
        }
        return group;
    }

    const GroupPotentialBucket * RequireBondGroup(GroupKey group_key, const std::string & class_key) const
    {
        auto * group{ TryGetBondGroup(group_key, class_key, true) };
        if (group == nullptr)
        {
            throw std::runtime_error("Bond group key is not available.");
        }
        return group;
    }

    const GroupPotentialEntry * FindAtomGroupEntry(const std::string & class_key) const
    {
        return ModelAnalysisData::Of(m_model_object).Atoms().FindGroupEntry(class_key);
    }

    const GroupPotentialEntry * FindBondGroupEntry(const std::string & class_key) const
    {
        return ModelAnalysisData::Of(m_model_object).Bonds().FindGroupEntry(class_key);
    }

    static void LogMissingAtomGroup(GroupKey group_key, const std::string & class_key)
    {
        std::ostringstream oss;
        if (class_key == ChemicalDataHelper::GetSimpleAtomClassKey())
        {
            oss << "Atom class group key : "
                << KeyPackerSimpleAtomClass::GetKeyString(group_key)
                << " not found.";
        }
        else if (class_key == ChemicalDataHelper::GetComponentAtomClassKey())
        {
            oss << "Atom class group key : "
                << KeyPackerComponentAtomClass::GetKeyString(group_key)
                << " not found.";
        }
        else if (class_key == ChemicalDataHelper::GetStructureAtomClassKey())
        {
            oss << "Atom class group key : "
                << KeyPackerStructureAtomClass::GetKeyString(group_key)
                << " not found.";
        }
        Logger::Log(LogLevel::Error, oss.str());
    }

    static void LogMissingBondGroup(GroupKey group_key, const std::string & class_key)
    {
        std::ostringstream oss;
        if (class_key == ChemicalDataHelper::GetSimpleBondClassKey())
        {
            oss << "Bond class group key : "
                << KeyPackerSimpleBondClass::GetKeyString(group_key)
                << " not found.";
        }
        else if (class_key == ChemicalDataHelper::GetComponentBondClassKey())
        {
            oss << "Bond class group key : "
                << KeyPackerComponentBondClass::GetKeyString(group_key)
                << " not found.";
        }
        Logger::Log(LogLevel::Error, oss.str());
    }

    static GaussianPosterior BuildPriorPosterior(const GroupPotentialBucket & group)
    {
        GaussianPosterior posterior;
        posterior.estimate = group.prior;
        posterior.variance = group.prior_variance;
        return posterior;
    }
};

} // namespace rhbm_gem
