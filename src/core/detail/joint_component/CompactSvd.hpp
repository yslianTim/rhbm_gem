#pragma once
#include "Numerics.hpp"
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include <functional>
#endif

namespace rhbm_gem::core::joint_component {
struct CompactSvdResult
{
    Vector singular_values,solution;
    Matrix right_vectors;
    Eigen::Index rank{};
    double threshold{};
    bool valid{},used_bdc{},jacobi_retry{};
};
enum class CompactSvdVectors {None,Right};
// Thresholds belong to the original problem, not the compact matrix dimensions.
// A non-null RHS requests a least-squares solution. Right vectors are retained
// only when explicitly requested; the default needs singular values alone.
CompactSvdResult CompactSvd(const Matrix &,double relative,double absolute=-1,const Vector * rhs=nullptr,
    CompactSvdVectors=CompactSvdVectors::None);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
enum class CompactSvdMode {Automatic,Legacy,ValuesOnly};
CompactSvdMode & CompactSvdModeForTesting();
using CompactSvdCapture=std::function<void(const Matrix &,double,double,const Vector *,const CompactSvdResult &)>;
CompactSvdCapture & CompactSvdCaptureForTesting();
#endif
}
