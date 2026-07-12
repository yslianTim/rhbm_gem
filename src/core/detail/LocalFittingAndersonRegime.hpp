#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace rhbm_gem::core::detail {

using LocalFittingAndersonRegimeClusterKey = std::vector<std::size_t>;

struct LocalFittingAndersonRegimeSignature
{
    double global_ridge_ratio{ 0.0 };
    std::vector<double> effective_ridge_multiplier_list{};
};

using LocalFittingAndersonRegimeSignatureMap =
    std::map<LocalFittingAndersonRegimeClusterKey, LocalFittingAndersonRegimeSignature>;

inline void ValidateLocalFittingAndersonRegimeSignature(
    const LocalFittingAndersonRegimeSignature & signature)
{
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
    return lhs.global_ridge_ratio == rhs.global_ridge_ratio &&
        lhs.effective_ridge_multiplier_list == rhs.effective_ridge_multiplier_list;
}

inline LocalFittingAndersonRegimeSignatureMap BuildLocalFittingAndersonRegimeSignatureMap(
    const std::vector<LocalFittingAndersonRegimeClusterKey> & cluster_key_list,
    const std::vector<std::size_t> & active_index_list,
    double global_ridge_ratio,
    const std::vector<double> & effective_ridge_multiplier_list)
{
    if (active_index_list.size() != effective_ridge_multiplier_list.size())
    {
        throw std::invalid_argument("Anderson regime active atom and ridge multiplier sizes differ.");
    }
    if (!std::isfinite(global_ridge_ratio) || global_ridge_ratio <= 0.0)
    {
        throw std::invalid_argument("Anderson regime global ridge ratio must be positive and finite.");
    }

    std::unordered_map<std::size_t, std::size_t> active_position_by_index;
    active_position_by_index.reserve(active_index_list.size());
    for (std::size_t position = 0; position < active_index_list.size(); position++)
    {
        const auto atom_index{ active_index_list.at(position) };
        if (!active_position_by_index.emplace(atom_index, position).second)
        {
            throw std::invalid_argument("Anderson regime active atom indexes must be unique.");
        }
        const auto multiplier{ effective_ridge_multiplier_list.at(position) };
        if (!std::isfinite(multiplier) || multiplier <= 0.0)
        {
            throw std::invalid_argument(
                "Anderson regime ridge multiplier must be positive and finite.");
        }
    }

    LocalFittingAndersonRegimeSignatureMap signature_by_key;
    std::size_t covered_atom_count{ 0 };
    std::vector<char> covered_active_position_list(active_index_list.size(), 0);
    for (const auto & key : cluster_key_list)
    {
        if (key.empty() || !std::is_sorted(key.begin(), key.end()) ||
            std::adjacent_find(key.begin(), key.end()) != key.end())
        {
            throw std::invalid_argument("Anderson regime cluster key must be non-empty and canonical.");
        }

        LocalFittingAndersonRegimeSignature signature;
        signature.global_ridge_ratio = global_ridge_ratio;
        signature.effective_ridge_multiplier_list.reserve(key.size());
        for (const auto atom_index : key)
        {
            const auto active_iter{ active_position_by_index.find(atom_index) };
            if (active_iter == active_position_by_index.end())
            {
                throw std::invalid_argument("Anderson regime cluster atom is not active.");
            }
            if (covered_active_position_list.at(active_iter->second) != 0)
            {
                throw std::invalid_argument(
                    "Anderson regime cluster keys must not share active atoms.");
            }
            covered_active_position_list.at(active_iter->second) = 1;
            signature.effective_ridge_multiplier_list.emplace_back(
                effective_ridge_multiplier_list.at(active_iter->second));
            covered_atom_count++;
        }
        if (!signature_by_key.emplace(key, std::move(signature)).second)
        {
            throw std::invalid_argument("Anderson regime cluster keys must be unique.");
        }
    }
    if (covered_atom_count != active_index_list.size())
    {
        throw std::invalid_argument("Anderson regime cluster keys must cover active atoms exactly once.");
    }
    return signature_by_key;
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
