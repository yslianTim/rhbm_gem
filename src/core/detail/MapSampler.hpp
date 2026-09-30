#pragma once

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

namespace rhbm_gem {
class MapObject;
namespace core::detail {
LocalPotentialSampleList BuildLocalPotentialSampleList(const MapObject &, const SamplingPointList &);
}
}
