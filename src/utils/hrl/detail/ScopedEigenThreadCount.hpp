#pragma once

#include <Eigen/Core>

namespace rhbm_gem::rhbm_helper::detail
{

// Scoped override for Eigen's global thread count. This changes process-wide
// Eigen behavior temporarily; it is not a per-operation thread hint.
class ScopedEigenThreadCount final
{
public:
    explicit ScopedEigenThreadCount(int requested_thread_count) :
        m_previous_thread_count{ Eigen::nbThreads() },
        m_requested_thread_count{
            (requested_thread_count <= 0) ? 1 : requested_thread_count
        },
        m_changed{ m_requested_thread_count != m_previous_thread_count }
    {
        if (m_changed)
        {
            Eigen::setNbThreads(m_requested_thread_count);
        }
    }

    ScopedEigenThreadCount(const ScopedEigenThreadCount &) = delete;
    ScopedEigenThreadCount & operator=(const ScopedEigenThreadCount &) = delete;
    ScopedEigenThreadCount(ScopedEigenThreadCount &&) = delete;
    ScopedEigenThreadCount & operator=(ScopedEigenThreadCount &&) = delete;

    ~ScopedEigenThreadCount()
    {
        if (m_changed && Eigen::nbThreads() != m_previous_thread_count)
        {
            Eigen::setNbThreads(m_previous_thread_count);
        }
    }

private:
    int m_previous_thread_count;
    int m_requested_thread_count;
    bool m_changed;
};

} // namespace rhbm_gem::rhbm_helper::detail
