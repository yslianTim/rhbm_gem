#pragma once

#include <algorithm>

#include <rhbm_gem/utils/algorithm/ParameterChangeStats.hpp>

namespace rhbm_gem::algorithm {

class AdaptiveRelaxationController
{
    double m_beta;
    double m_beta_min{ 0.0 };
    double m_beta_max{ 1.0 };
    double m_growth_factor{ 1.0 };
    double m_shrink_factor{ 1.0 };
    double m_improvement_ratio{ 0.0 };
    double m_previous_change{ 0.0 };
    int m_improvement_streak{ 0 };
    int m_increase_streak_size{ 0 };
    bool m_has_previous_change{ false };

public:
    AdaptiveRelaxationController(
        double initial_beta,
        double beta_min,
        double beta_max,
        double growth_factor,
        double shrink_factor,
        double improvement_ratio,
        int increase_streak_size)
        : m_beta{ std::clamp(initial_beta, beta_min, beta_max) },
          m_beta_min{ beta_min },
          m_beta_max{ beta_max },
          m_growth_factor{ growth_factor },
          m_shrink_factor{ shrink_factor },
          m_improvement_ratio{ improvement_ratio },
          m_increase_streak_size{ increase_streak_size }
    {
    }

    double GetBeta() const
    {
        return m_beta;
    }

    bool IsAtMinimum() const
    {
        return m_beta <= m_beta_min;
    }

    double Shrink()
    {
        m_beta = std::max(m_beta_min, m_beta * m_shrink_factor);
        m_improvement_streak = 0;
        return m_beta;
    }

    void Update(const ParameterChangeStats & stats)
    {
        Update(GetMaximumParameterChange(stats));
    }

    void Update(double change)
    {
        if (!m_has_previous_change)
        {
            m_previous_change = change;
            m_has_previous_change = true;
            return;
        }

        if (change > m_previous_change * (1.0 + m_improvement_ratio))
        {
            Shrink();
        }
        else if (change < m_previous_change * (1.0 - m_improvement_ratio))
        {
            m_improvement_streak++;
            if (m_improvement_streak >= m_increase_streak_size)
            {
                m_beta = std::min(m_beta_max, m_beta * m_growth_factor);
                m_improvement_streak = 0;
            }
        }
        else
        {
            m_improvement_streak = 0;
        }
        m_previous_change = change;
    }
    
};

} // namespace rhbm_gem::algorithm
