#include <gtest/gtest.h>

#include "utils/hrl/detail/ScopedEigenThreadCount.hpp"

namespace
{
class ScopedEigenThreadCountRestorer final
{
public:
    explicit ScopedEigenThreadCountRestorer(int thread_count) :
        m_thread_count{ thread_count }
    {
    }

    ~ScopedEigenThreadCountRestorer()
    {
        Eigen::setNbThreads(m_thread_count);
    }

private:
    int m_thread_count;
};
} // namespace

TEST(ScopedEigenThreadCountTest, NonPositiveRequestFallsBackToOneThreadAndRestoresPreviousValue)
{
    const int original_thread_count{ Eigen::nbThreads() };
    ScopedEigenThreadCountRestorer restore_original_thread_count(original_thread_count);

    {
        HRLModelAlgorithms::detail::ScopedEigenThreadCount scoped_thread_count(0);
        EXPECT_EQ(1, Eigen::nbThreads());
    }

    EXPECT_EQ(original_thread_count, Eigen::nbThreads());
}

TEST(ScopedEigenThreadCountTest, MatchingRequestKeepsObservableThreadCountUnchanged)
{
    const int original_thread_count{ Eigen::nbThreads() };
    ScopedEigenThreadCountRestorer restore_original_thread_count(original_thread_count);

    {
        HRLModelAlgorithms::detail::ScopedEigenThreadCount scoped_thread_count(original_thread_count);
        EXPECT_EQ(original_thread_count, Eigen::nbThreads());
    }

    EXPECT_EQ(original_thread_count, Eigen::nbThreads());
}

TEST(ScopedEigenThreadCountTest, ScopeRestoresPreviousThreadCountAfterChange)
{
    const int original_thread_count{ Eigen::nbThreads() };
    ScopedEigenThreadCountRestorer restore_original_thread_count(original_thread_count);
    const int alternate_thread_count{ (original_thread_count == 1) ? 2 : 1 };

    {
        HRLModelAlgorithms::detail::ScopedEigenThreadCount scoped_thread_count(alternate_thread_count);
        EXPECT_EQ(alternate_thread_count, Eigen::nbThreads());
    }

    EXPECT_EQ(original_thread_count, Eigen::nbThreads());
}

TEST(ScopedEigenThreadCountTest, NestedScopesRestoreOuterThenOriginalThreadCount)
{
    const int original_thread_count{ Eigen::nbThreads() };
    ScopedEigenThreadCountRestorer restore_original_thread_count(original_thread_count);
    const int outer_thread_count{ (original_thread_count == 1) ? 2 : 1 };
    const int inner_thread_count{
        (outer_thread_count == 2) ? 3 : 2
    };

    {
        HRLModelAlgorithms::detail::ScopedEigenThreadCount outer_scope(outer_thread_count);
        EXPECT_EQ(outer_thread_count, Eigen::nbThreads());

        {
            HRLModelAlgorithms::detail::ScopedEigenThreadCount inner_scope(inner_thread_count);
            EXPECT_EQ(inner_thread_count, Eigen::nbThreads());
        }

        EXPECT_EQ(outer_thread_count, Eigen::nbThreads());
    }

    EXPECT_EQ(original_thread_count, Eigen::nbThreads());
}
