#pragma once

#include <Eigen/Core>

namespace rhbm_gem::core::detail {

class ScopedEigenThreadCount
{
public:
    explicit ScopedEigenThreadCount(int requested_thread_count)
        : m_previous_thread_count{ Eigen::nbThreads() },
          m_changed{ requested_thread_count > 0 &&
              m_previous_thread_count != requested_thread_count }
    {
        if (m_changed)
        {
            Eigen::setNbThreads(requested_thread_count);
        }
    }

    ScopedEigenThreadCount(const ScopedEigenThreadCount &) = delete;
    ScopedEigenThreadCount & operator=(const ScopedEigenThreadCount &) = delete;

    ~ScopedEigenThreadCount()
    {
        if (m_changed && Eigen::nbThreads() != m_previous_thread_count)
        {
            Eigen::setNbThreads(m_previous_thread_count);
        }
    }

private:
    int m_previous_thread_count{ 1 };
    bool m_changed{ false };
};

} // namespace rhbm_gem::core::detail
