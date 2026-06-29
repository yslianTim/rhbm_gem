#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <rhbm_gem/utils/algorithm/ParameterChange.hpp>

namespace rhbm_gem::algorithm {

class ConvergenceFreezeTracker
{
    std::vector<bool> m_frozen_list;
    std::vector<int> m_stable_count_list;
    double m_freeze_threshold{ 0.0 };
    int m_stable_iteration_size{ 0 };

public:
    ConvergenceFreezeTracker(
        std::size_t value_size,
        double change_tolerance,
        double freeze_change_ratio,
        int stable_iteration_size)
        : m_frozen_list(value_size, false),
          m_stable_count_list(value_size, 0),
          m_freeze_threshold{ std::sqrt(change_tolerance) * freeze_change_ratio },
          m_stable_iteration_size{ stable_iteration_size }
    {
    }

    std::vector<std::size_t> BuildActiveIndexList() const
    {
        std::vector<std::size_t> active_index_list;
        active_index_list.reserve(m_frozen_list.size() - GetFrozenCount());
        for (std::size_t i = 0; i < m_frozen_list.size(); i++)
        {
            if (!m_frozen_list.at(i))
            {
                active_index_list.emplace_back(i);
            }
        }
        return active_index_list;
    }

    void Update(
        const std::vector<ParameterChange> & change_list,
        const std::vector<std::size_t> & active_index_list)
    {
        if (change_list.size() != m_frozen_list.size())
        {
            throw std::invalid_argument("Convergence freeze tracker input size is inconsistent.");
        }

        for (const auto index : active_index_list)
        {
            if (index >= m_frozen_list.size())
            {
                throw std::invalid_argument("Convergence freeze tracker active index is out of range.");
            }
            if (m_frozen_list.at(index)) continue;

            if (GetMaximumParameterChange(change_list.at(index)) < m_freeze_threshold)
            {
                m_stable_count_list.at(index)++;
                if (m_stable_count_list.at(index) >= m_stable_iteration_size)
                {
                    m_frozen_list.at(index) = true;
                }
            }
            else
            {
                m_stable_count_list.at(index) = 0;
            }
        }
    }

    bool IsFrozen(std::size_t index) const
    {
        if (index >= m_frozen_list.size())
        {
            throw std::invalid_argument("Convergence freeze tracker index is out of range.");
        }
        return m_frozen_list.at(index);
    }

    bool Thaw(std::size_t index)
    {
        if (index >= m_frozen_list.size())
        {
            throw std::invalid_argument("Convergence freeze tracker index is out of range.");
        }
        const bool was_frozen{ m_frozen_list.at(index) };
        m_frozen_list.at(index) = false;
        m_stable_count_list.at(index) = 0;
        return was_frozen;
    }

    std::size_t GetFrozenCount() const
    {
        return static_cast<std::size_t>(
            std::count(m_frozen_list.begin(), m_frozen_list.end(), true));
    }

    std::size_t GetActiveCount() const
    {
        return m_frozen_list.size() - GetFrozenCount();
    }
    
};

} // namespace rhbm_gem::algorithm
