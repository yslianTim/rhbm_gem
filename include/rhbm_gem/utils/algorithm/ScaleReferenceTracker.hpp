#pragma once

#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <vector>

namespace rhbm_gem::algorithm {

struct ScaleReference
{
    bool has_reference{ false };
    double scale{ std::numeric_limits<double>::infinity() };
};

class ScaleReferenceTracker
{
    std::size_t m_warmup_sample_count{ 0 };
    std::vector<double> m_scale_sample_list{};
    std::size_t m_accepted_scale_sample_count{ 0 };
    bool m_locked{ false };

    bool AddScaleSample(double scale_sample)
    {
        if (!std::isfinite(scale_sample)) return false;

        m_scale_sample_list.emplace_back(scale_sample);
        while (m_scale_sample_list.size() > m_warmup_sample_count)
        {
            m_scale_sample_list.erase(m_scale_sample_list.begin());
        }
        return true;
    }

    static ScaleReference BuildReference(const std::vector<double> & scale_sample_list)
    {
        ScaleReference reference;
        if (scale_sample_list.empty())
        {
            return reference;
        }

        const auto scale_sum{
            std::accumulate(scale_sample_list.begin(), scale_sample_list.end(), 0.0)
        };
        const auto scale{
            scale_sum / static_cast<double>(scale_sample_list.size())
        };
        if (!std::isfinite(scale))
        {
            return reference;
        }

        reference.has_reference = true;
        reference.scale = scale;
        return reference;
    }

public:
    explicit ScaleReferenceTracker(
        std::size_t warmup_sample_count,
        std::optional<double> initial_scale_sample = std::nullopt)
        : m_warmup_sample_count{ warmup_sample_count }
    {
        if (m_warmup_sample_count == 0)
        {
            throw std::invalid_argument("Scale reference warmup sample count must be positive.");
        }
        if (initial_scale_sample.has_value())
        {
            AddScaleSample(*initial_scale_sample);
        }
    }

    bool HasReference() const
    {
        return !m_scale_sample_list.empty();
    }

    bool IsLocked() const
    {
        return m_locked;
    }

    ScaleReference GetCommittedReference() const
    {
        return BuildReference(m_scale_sample_list);
    }

    ScaleReference GetProvisionalReference(double scale_sample) const
    {
        if (m_locked)
        {
            return GetCommittedReference();
        }

        auto scale_sample_list{ m_scale_sample_list };
        if (std::isfinite(scale_sample))
        {
            scale_sample_list.emplace_back(scale_sample);
            while (scale_sample_list.size() > m_warmup_sample_count)
            {
                scale_sample_list.erase(scale_sample_list.begin());
            }
        }
        return BuildReference(scale_sample_list);
    }

    void CommitScaleSample(double scale_sample)
    {
        if (m_locked) return;

        if (!AddScaleSample(scale_sample)) return;
        m_accepted_scale_sample_count++;
        if (m_accepted_scale_sample_count >= m_warmup_sample_count)
        {
            m_locked = true;
        }
    }
};

} // namespace rhbm_gem::algorithm
