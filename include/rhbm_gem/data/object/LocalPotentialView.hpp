#pragma once

#include <stdexcept>

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/LocalPotentialEntry.hpp>

namespace rhbm_gem {

class LocalPotentialView
{
    const LocalPotentialEntry & m_local_entry;

public:
    explicit LocalPotentialView(const AtomObject & atom_object) :
        m_local_entry(*RequireLocalEntry(atom_object.GetLocalPotentialEntry(), "Atom local entry")) {}

    explicit LocalPotentialView(const BondObject & bond_object) :
        m_local_entry(*RequireLocalEntry(bond_object.GetLocalPotentialEntry(), "Bond local entry")) {}

    explicit LocalPotentialView(const LocalPotentialEntry & local_entry) :
        m_local_entry(local_entry) {}

    const LocalPotentialEntry & GetLocalEntry() const { return m_local_entry; }

    bool IsOutlier(const std::string & class_key) const
    {
        return m_local_entry.GetOutlierTag(class_key);
    }

    const std::vector<std::tuple<float, float>> & GetDistanceAndMapValueList() const
    {
        return m_local_entry.GetDistanceAndMapValueList();
    }

    const GaussianEstimate & GetEstimateMDPDE() const { return m_local_entry.GetEstimateMDPDE(); }
    const GaussianPosterior & GetPosterior(const std::string & key) const { return m_local_entry.GetPosterior(key); }
    double GetAlphaR() const { return m_local_entry.GetAlphaR(); }

private:
    static const LocalPotentialEntry * RequireLocalEntry(
        const LocalPotentialEntry * entry,
        const char * context)
    {
        if (entry == nullptr)
        {
            throw std::runtime_error(std::string(context) + " is not available.");
        }
        return entry;
    }
};

} // namespace rhbm_gem
