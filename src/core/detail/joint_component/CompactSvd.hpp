#pragma once
#include "Numerics.hpp"

namespace rhbm_gem::core::joint_component {
struct CompactSvdResult
{
    Vector singular_values;
    Matrix right_vectors;
    Eigen::Index rank{};
    double threshold{};
    bool valid{},used_bdc{},jacobi_retry{};
};
// Preserve the historical distinction: spectrum counts use strict >, while
// solve/derivative evidence uses Eigen's native cutoff and zero-value guard.
enum class RankBoundary {SvdNative,StrictGreater};
struct RankRequest
{
    RankPolicy policy;
    Eigen::Index columns{}; // Decision dimensions may be larger than a component compact.
    double absolute{-1};
    RankBoundary boundary{RankBoundary::SvdNative};
};
enum class CompactSvdVectors {None,Right};
// Thresholds belong to the original problem, not the compact matrix dimensions.
// Right vectors are retained only when explicitly requested.
CompactSvdResult CompactSvd(const Matrix &,double relative,double absolute=-1,
    CompactSvdVectors=CompactSvdVectors::None);
CompactSvdResult EvaluateRank(const Matrix &,const RankRequest &,
    CompactSvdVectors=CompactSvdVectors::None);
bool CompactRankDecisionSafe(const CompactSvdResult &);
}
