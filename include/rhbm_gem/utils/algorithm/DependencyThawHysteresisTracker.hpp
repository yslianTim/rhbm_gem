#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace rhbm_gem::algorithm {

class DependencyThawHysteresisTracker
{
    std::vector<double> m_multiplier_list;
    double m_growth_multiplier{ 1.0 };
    double m_max_multiplier{ 1.0 };
    double m_frozen_decay{ 1.0 };

    void ValidateIndex(std::size_t index) const
    {
        if (index >= m_multiplier_list.size())
        {
            throw std::invalid_argument("Dependency thaw hysteresis index is out of range.");
        }
    }

public:
    DependencyThawHysteresisTracker(
        std::size_t value_size,
        double growth_multiplier,
        double max_multiplier,
        double frozen_decay)
        : m_multiplier_list(value_size, 1.0),
          m_growth_multiplier{ growth_multiplier },
          m_max_multiplier{ max_multiplier },
          m_frozen_decay{ frozen_decay }
    {
        if (m_growth_multiplier < 1.0 || m_max_multiplier < 1.0 ||
            m_frozen_decay < 0.0 || m_frozen_decay > 1.0)
        {
            throw std::invalid_argument("Dependency thaw hysteresis settings are invalid.");
        }
    }

    double GetThreshold(std::size_t index, double base_threshold) const
    {
        ValidateIndex(index);
        return base_threshold * m_multiplier_list.at(index);
    }

    bool ShouldThaw(std::size_t index, double change, double base_threshold) const
    {
        return change >= GetThreshold(index, base_threshold);
    }

    void RecordDependencyThaw(std::size_t index)
    {
        ValidateIndex(index);
        m_multiplier_list.at(index) = std::min(
            m_max_multiplier,
            m_multiplier_list.at(index) * m_growth_multiplier);
    }

    void DecayFrozen(std::size_t index)
    {
        ValidateIndex(index);
        m_multiplier_list.at(index) = std::max(1.0, m_multiplier_list.at(index) * m_frozen_decay);
    }
};

} // namespace rhbm_gem::algorithm
