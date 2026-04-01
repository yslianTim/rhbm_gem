#pragma once

#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/data/object/AtomClassifier.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondClassifier.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/GroupPotentialEntry.hpp>
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

    const ModelObject & GetModelObject() const { return m_model_object; }

    double GetAtomGausEstimateMinimum(int par_id, Element element) const
    {
        std::vector<double> gaus_estimate_list;
        gaus_estimate_list.reserve(m_model_object.GetNumberOfSelectedAtom());
        for (const auto * atom : m_model_object.GetSelectedAtomList())
        {
            if (atom->GetElement() != element) continue;
            auto * entry{ atom->GetLocalPotentialEntry() };
            gaus_estimate_list.emplace_back(entry->GetEstimateMDPDE().GetParameter(par_id));
        }
        return ArrayStats<double>::ComputeMin(gaus_estimate_list.data(), gaus_estimate_list.size());
    }

    double GetBondGausEstimateMinimum(int par_id) const
    {
        std::vector<double> gaus_estimate_list;
        gaus_estimate_list.reserve(m_model_object.GetNumberOfSelectedBond());
        for (const auto * bond : m_model_object.GetSelectedBondList())
        {
            auto * entry{ bond->GetLocalPotentialEntry() };
            gaus_estimate_list.emplace_back(entry->GetEstimateMDPDE().GetParameter(par_id));
        }
        return ArrayStats<double>::ComputeMin(gaus_estimate_list.data(), gaus_estimate_list.size());
    }

    bool HasAtomGroup(GroupKey group_key, const std::string & class_key, bool verbose=false) const
    {
        return CheckAtomGroupKey(group_key, class_key, verbose);
    }

    bool HasBondGroup(GroupKey group_key, const std::string & class_key, bool verbose=false) const
    {
        return CheckBondGroupKey(group_key, class_key, verbose);
    }

    double GetAtomGausEstimateMean(GroupKey group_key, const std::string & class_key, int par_id) const
    {
        return GetAtomGroup(group_key, class_key).mean.GetParameter(par_id);
    }

    double GetAtomGausEstimatePrior(GroupKey group_key, const std::string & class_key, int par_id) const
    {
        return GetAtomGroup(group_key, class_key).prior.GetParameter(par_id);
    }

    double GetBondGausEstimatePrior(GroupKey group_key, const std::string & class_key, int par_id) const
    {
        return GetBondGroup(group_key, class_key).prior.GetParameter(par_id);
    }

    double GetAtomGausVariancePrior(GroupKey group_key, const std::string & class_key, int par_id) const
    {
        return BuildPriorPosterior(GetAtomGroup(group_key, class_key)).GetVariance(par_id);
    }

    double GetBondGausVariancePrior(GroupKey group_key, const std::string & class_key, int par_id) const
    {
        return BuildPriorPosterior(GetBondGroup(group_key, class_key)).GetVariance(par_id);
    }

    const std::vector<AtomObject *> & GetAtomObjectList(GroupKey group_key, const std::string & class_key) const
    {
        return GetAtomGroup(group_key, class_key).atom_members;
    }

    const std::vector<BondObject *> & GetBondObjectList(GroupKey group_key, const std::string & class_key) const
    {
        return GetBondGroup(group_key, class_key).bond_members;
    }

    std::vector<AtomObject *> GetOutlierAtomObjectList(GroupKey group_key, const std::string & class_key) const
    {
        std::vector<AtomObject *> outlier_atom_list;
        const auto & atom_list{ GetAtomObjectList(group_key, class_key) };
        outlier_atom_list.reserve(atom_list.size());
        for (auto * atom : atom_list)
        {
            if (atom->GetLocalPotentialEntry()->GetOutlierTag(class_key))
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
        auto * local_entry{ atom_list.front()->GetLocalPotentialEntry() };
        if (local_entry == nullptr)
        {
            throw std::runtime_error("Local entry of the first atom member is not available.");
        }
        return local_entry->GetAlphaR();
    }

    double GetBondAlphaR(GroupKey group_key, const std::string & class_key) const
    {
        const auto & bond_list{ GetBondObjectList(group_key, class_key) };
        if (bond_list.empty())
        {
            throw std::runtime_error("Bond group has no members.");
        }
        auto * local_entry{ bond_list.front()->GetLocalPotentialEntry() };
        if (local_entry == nullptr)
        {
            throw std::runtime_error("Local entry of the first bond member is not available.");
        }
        return local_entry->GetAlphaR();
    }

    double GetAtomAlphaG(GroupKey group_key, const std::string & class_key) const
    {
        return GetAtomGroup(group_key, class_key).alpha_g;
    }

    double GetBondAlphaG(GroupKey group_key, const std::string & class_key) const
    {
        return GetBondGroup(group_key, class_key).alpha_g;
    }

    Residue GetResidueFromAtomGroupKey(GroupKey group_key, const std::string & class_key) const
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

    Residue GetResidueFromBondGroupKey(GroupKey group_key, const std::string & class_key) const
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

private:
    const GroupPotentialEntry & GetAtomGroupEntry(const std::string & class_key) const
    {
        return *m_model_object.GetAtomGroupPotentialEntry(class_key);
    }

    const GroupPotentialEntry & GetBondGroupEntry(const std::string & class_key) const
    {
        return *m_model_object.GetBondGroupPotentialEntry(class_key);
    }

    const GroupPotentialBucket & GetAtomGroup(GroupKey group_key, const std::string & class_key) const
    {
        if (!CheckAtomGroupKey(group_key, class_key, true))
        {
            throw std::runtime_error("Atom group key is not available.");
        }
        return GetAtomGroupEntry(class_key).GetGroup(group_key);
    }

    const GroupPotentialBucket & GetBondGroup(GroupKey group_key, const std::string & class_key) const
    {
        if (!CheckBondGroupKey(group_key, class_key, true))
        {
            throw std::runtime_error("Bond group key is not available.");
        }
        return GetBondGroupEntry(class_key).GetGroup(group_key);
    }

    bool CheckAtomGroupKey(GroupKey group_key, const std::string & class_key, bool verbose) const
    {
        if (GetAtomGroupEntry(class_key).HasGroup(group_key))
        {
            return true;
        }
        if (verbose)
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
        return false;
    }

    bool CheckBondGroupKey(GroupKey group_key, const std::string & class_key, bool verbose) const
    {
        if (GetBondGroupEntry(class_key).HasGroup(group_key))
        {
            return true;
        }
        if (verbose)
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
        return false;
    }

    GaussianPosterior BuildPriorPosterior(const GroupPotentialBucket & group) const
    {
        GaussianPosterior posterior;
        posterior.estimate = group.prior;
        posterior.variance = group.prior_variance;
        return posterior;
    }
};

} // namespace rhbm_gem
