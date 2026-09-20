#pragma once
#include <rhbm_gem/data/object/JointAnalysisResult.hpp>
#include <string_view>

namespace rhbm_gem::joint_result_io {
std::string_view StatusText(JointCheckStatus);
std::string Encode(const JointAnalysisResult &);
JointAnalysisResult Decode(std::string_view);
}
