#pragma once

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include <functional>
#endif

namespace rhbm_gem {
class MapObject;
namespace core::detail {
LocalPotentialSampleList BuildLocalPotentialSampleList(const MapObject &, const SamplingPointList &);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
using PotentialSamplingObserver = std::function<void(int)>;
PotentialSamplingObserver & PotentialSamplingObserverForTesting();
#endif
}
}
