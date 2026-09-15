#pragma once
#include <rhbm_gem/utils/domain/SamplingTypes.hpp>
#include <string>
namespace rhbm_gem { class MapObject; }
namespace second_stage_test {
LocalPotentialSampleList SampleExperimentPoints(const rhbm_gem::MapObject &, const SamplingPointList &);
void RunForwardExperiment(const std::string & manifest, const std::string & map,
    const std::string & capture_directory, const std::string & output_directory);
}
