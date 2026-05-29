#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/domain/SamplingTypes.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

namespace rhbm_gem::core {
struct LocalTestData
{
    GaussianModel3D gaus_true;
    std::vector<LocalPotentialSampleList> replica_sampling_entries;
};

struct GroupTestData
{
    GaussianModel3D gaus_true;
    std::vector<std::vector<LocalPotentialSampleList>> replica_member_sampling_entries;
};

struct LocalScenario
{
    GaussianModel3D gaus_true;
    int sampling_entry_size{ 1 };
    double data_error_sigma{ 1.0 };
    double outlier_ratio{ 0.0 };
    int replica_size{ 1 };
    std::optional<std::uint32_t> random_seed{};
};

struct GaussianParameterDistribution
{
    GaussianModel3D mean;
    GaussianModel3DUncertainty sigma;
};

struct GroupScenario
{
    int member_size{ 1 };
    int sampling_entry_size{ 10 };
    GaussianParameterDistribution inlier_distribution;
    GaussianParameterDistribution outlier_distribution;
    double outlier_ratio{ 0.0 };
    int replica_size{ 1 };
    std::optional<std::uint32_t> random_seed{};
};

struct AtomModelScenario
{
    Spot spot{ Spot::UNK };
    GaussianModel3D gaus_true;
    double data_error_sigma{ 1.0 };
    int replica_size{ 1 };
    std::optional<std::uint32_t> random_seed{};
};

LocalTestData BuildLocalTestData(const LocalScenario & scenario);
LocalTestData BuildLocalTestData(const AtomModelScenario & scenario);
GroupTestData BuildGroupTestData(const GroupScenario & scenario);

} // namespace rhbm_gem::core
