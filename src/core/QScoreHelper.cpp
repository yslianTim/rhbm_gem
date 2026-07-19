#include <rhbm_gem/core/QScoreHelper.hpp>

#include <algorithm>

#include <rhbm_gem/data/object/MapObject.hpp>

namespace rhbm_gem::core {

std::tuple<float, float> GetReferenceGaussianParameters(const MapObject & map_object)
{
    const auto reference_high{
        std::min(
            map_object.GetMapValueMean() + 10.0f * map_object.GetMapValueSD(),
            map_object.GetMapValueMax())
    };
    const auto offset{
        std::max(
            map_object.GetMapValueMean() - map_object.GetMapValueSD(),
            map_object.GetMapValueMin())
    };
    return std::make_tuple(reference_high - offset, offset);
}

} // namespace rhbm_gem::core
