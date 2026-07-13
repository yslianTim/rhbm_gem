#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <vector>

#include "core/detail/LocalFittingHealth.hpp"

namespace rhbm_gem::core::detail {

using LocalFittingAndersonRegimeClusterKey = std::vector<std::size_t>;

struct LocalFittingAndersonRegimeSignature
{
    JointOffsetSolveStatus joint_offset_status{ JointOffsetSolveStatus::SystemBuildFailed };
    double global_ridge_ratio{ 0.0 };
    std::vector<double> effective_ridge_multiplier_list{};
};

using LocalFittingAndersonRegimeSignatureMap =
    std::map<LocalFittingAndersonRegimeClusterKey, LocalFittingAndersonRegimeSignature>;

inline void ValidateLocalFittingAndersonRegimeSignature(
    const LocalFittingAndersonRegimeSignature & signature)
{
    static_cast<void>(
        IsJointOffsetSolveProgressEligible(signature.joint_offset_status));
    if (!std::isfinite(signature.global_ridge_ratio) || signature.global_ridge_ratio <= 0.0)
    {
        throw std::invalid_argument("Anderson regime global ridge ratio must be positive and finite.");
    }
    if (signature.effective_ridge_multiplier_list.empty())
    {
        throw std::invalid_argument("Anderson regime ridge multiplier list must not be empty.");
    }
    for (const auto multiplier : signature.effective_ridge_multiplier_list)
    {
        if (!std::isfinite(multiplier) || multiplier <= 0.0)
        {
            throw std::invalid_argument(
                "Anderson regime ridge multiplier must be positive and finite.");
        }
    }
}

inline bool AreLocalFittingAndersonRegimeSignaturesEqual(
    const LocalFittingAndersonRegimeSignature & lhs,
    const LocalFittingAndersonRegimeSignature & rhs)
{
    return lhs.joint_offset_status == rhs.joint_offset_status &&
        lhs.global_ridge_ratio == rhs.global_ridge_ratio &&
        lhs.effective_ridge_multiplier_list == rhs.effective_ridge_multiplier_list;
}

class LocalFittingAndersonRegimeTracker
{
    LocalFittingAndersonRegimeSignatureMap m_committed_signature_by_key{};

public:
    void Reconcile(const std::vector<LocalFittingAndersonRegimeClusterKey> & cluster_key_list)
    {
        LocalFittingAndersonRegimeSignatureMap reconciled_signature_by_key;
        for (const auto & key : cluster_key_list)
        {
            const auto iter{ m_committed_signature_by_key.find(key) };
            if (iter == m_committed_signature_by_key.end()) continue;
            reconciled_signature_by_key.emplace(key, std::move(iter->second));
        }
        m_committed_signature_by_key = std::move(reconciled_signature_by_key);
    }

    std::vector<LocalFittingAndersonRegimeClusterKey> FindIncompatible(
        const LocalFittingAndersonRegimeSignatureMap & current_signature_by_key) const
    {
        std::vector<LocalFittingAndersonRegimeClusterKey> incompatible_key_list;
        for (const auto & [key, current_signature] : current_signature_by_key)
        {
            ValidateLocalFittingAndersonRegimeSignature(current_signature);
            if (current_signature.effective_ridge_multiplier_list.size() != key.size())
            {
                throw std::invalid_argument(
                    "Anderson regime signature size must match its cluster key.");
            }
            const auto committed_iter{ m_committed_signature_by_key.find(key) };
            if (committed_iter == m_committed_signature_by_key.end()) continue;
            if (!AreLocalFittingAndersonRegimeSignaturesEqual(
                    committed_iter->second, current_signature))
            {
                incompatible_key_list.emplace_back(key);
            }
        }
        return incompatible_key_list;
    }

    void Invalidate(const std::vector<LocalFittingAndersonRegimeClusterKey> & key_list)
    {
        for (const auto & key : key_list)
        {
            m_committed_signature_by_key.erase(key);
        }
    }

    void InvalidateContaining(const std::vector<std::size_t> & atom_index_list)
    {
        for (auto iter = m_committed_signature_by_key.begin();
             iter != m_committed_signature_by_key.end();)
        {
            const auto contains_affected_atom{
                std::any_of(
                    atom_index_list.begin(), atom_index_list.end(),
                    [&](std::size_t atom_index)
                    {
                        return std::binary_search(
                            iter->first.begin(), iter->first.end(), atom_index);
                    })
            };
            if (contains_affected_atom)
            {
                iter = m_committed_signature_by_key.erase(iter);
            }
            else
            {
                ++iter;
            }
        }
    }

    void Commit(
        const std::vector<LocalFittingAndersonRegimeClusterKey> & key_list,
        const LocalFittingAndersonRegimeSignatureMap & current_signature_by_key)
    {
        for (const auto & key : key_list)
        {
            const auto signature_iter{ current_signature_by_key.find(key) };
            if (signature_iter == current_signature_by_key.end())
            {
                throw std::invalid_argument("Anderson regime signature is missing for committed cluster.");
            }
            ValidateLocalFittingAndersonRegimeSignature(signature_iter->second);
            if (signature_iter->second.effective_ridge_multiplier_list.size() != key.size())
            {
                throw std::invalid_argument(
                    "Anderson regime signature size must match its cluster key.");
            }
            m_committed_signature_by_key[key] = signature_iter->second;
        }
    }
};

} // namespace rhbm_gem::core::detail
