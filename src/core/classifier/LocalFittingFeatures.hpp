#pragma once

#include "core/classifier/LocalFittingFeatureNames.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace rhbm_gem {

class ModelObject;
class GaussianModel3D;

} // namespace rhbm_gem

namespace rhbm_gem::core::classifier {

inline constexpr std::size_t kLocalFittingIdentifierColumnCount{ 3 };
inline constexpr std::size_t kLocalFittingColumnCount{
    kLocalFittingIdentifierColumnCount + kLocalFittingFeatureCount
};

inline constexpr std::array<std::string_view, kLocalFittingIdentifierColumnCount>
    kLocalFittingIdentifierNames{
        "serial id",
        "residue",
        "spot",
    };

struct LocalFittingFeatureRow
{
    int serial_id{ 0 };
    std::string residue{};
    std::string spot{};
    std::array<double, kLocalFittingFeatureCount> features{};
    std::string tail_peeling_ratio_unavailable_reason;
};

struct JointModelResponseFeatures
{
    double center_height;
    double gaussian_integral;
    double charge_integral;
};

JointModelResponseFeatures EvaluateJointModelResponseFeatures(const GaussianModel3D & point);

std::string BuildLocalFittingCsvHeader();

std::vector<LocalFittingFeatureRow> BuildLocalFittingFeatureRows(
    const ModelObject & model_object);

} // namespace rhbm_gem::core::classifier
