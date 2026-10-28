#pragma once
#include <Eigen/Dense>
#include <algorithm>
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include <chrono>
#endif

namespace rhbm_gem::core::joint_component {
inline constexpr Eigen::Index derivative_tile_rows=8192;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
struct TiledQRAppendTiming
{
    double dense_assembly_seconds{},householder_factorization_seconds{},qtrhs_transform_seconds{},
        r_extraction_seconds{};
};
#endif
// Consume original rows once. The discarded orthogonal RHS tail is never used
// as an objective: callers retain and evaluate the full residual separately.
struct TiledQR
{
    Eigen::MatrixXd r,target;
    Eigen::MatrixXd scratch_rows,scratch_rhs;
    TiledQR(Eigen::Index columns,Eigen::Index responses):r(0,columns),target(0,responses) {}
    void Append(const Eigen::MatrixXd & rows,const Eigen::MatrixXd & rhs,bool reference_order=false
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        ,TiledQRAppendTiming * timing=nullptr
#endif
    )
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        const auto assembly_started=timing ? std::chrono::steady_clock::now() :
            std::chrono::steady_clock::time_point{};
#endif
        const auto prior=r.rows();
        const auto assembled_rows=prior+rows.rows();
        if(scratch_rows.rows()!=assembled_rows || scratch_rows.cols()!=r.cols())
            scratch_rows.resize(assembled_rows,r.cols());
        if(scratch_rhs.rows()!=assembled_rows || scratch_rhs.cols()!=target.cols())
            scratch_rhs.resize(assembled_rows,target.cols());
        auto & a=scratch_rows;
        auto & b=scratch_rhs;
        a.topRows(prior)=r; a.bottomRows(rows.rows())=rows;
        b.topRows(prior)=target; b.bottomRows(rows.rows())=rhs;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(timing) timing->dense_assembly_seconds+=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-assembly_started).count();
#endif
        if(reference_order)
        {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            const auto factor_started=timing ? std::chrono::steady_clock::now() :
                std::chrono::steady_clock::time_point{};
#endif
            const Eigen::HouseholderQR<Eigen::MatrixXd> qr(a);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            if(timing) timing->householder_factorization_seconds+=std::chrono::duration<double>(
                std::chrono::steady_clock::now()-factor_started).count();
            const auto transform_started=timing ? std::chrono::steady_clock::now() :
                std::chrono::steady_clock::time_point{};
#endif
            const Eigen::MatrixXd transformed=qr.householderQ().adjoint()*b;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            if(timing) timing->qtrhs_transform_seconds+=std::chrono::duration<double>(
                std::chrono::steady_clock::now()-transform_started).count();
            const auto extraction_started=timing ? std::chrono::steady_clock::now() :
                std::chrono::steady_clock::time_point{};
#endif
            const auto keep=std::min(a.rows(),a.cols());
            r=qr.matrixQR().topRows(keep).triangularView<Eigen::Upper>();
            target=transformed.topRows(keep);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            if(timing) timing->r_extraction_seconds+=std::chrono::duration<double>(
                std::chrono::steady_clock::now()-extraction_started).count();
#endif
            return;
        }
        // The assembled tile is disposable. Factor and transform it in place;
        // retain the same Householder arithmetic without two full-size copies.
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        const auto factor_started=timing ? std::chrono::steady_clock::now() :
            std::chrono::steady_clock::time_point{};
#endif
        const Eigen::HouseholderQR<Eigen::Ref<Eigen::MatrixXd>> qr(a);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(timing) timing->householder_factorization_seconds+=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-factor_started).count();
        const auto transform_started=timing ? std::chrono::steady_clock::now() :
            std::chrono::steady_clock::time_point{};
#endif
        b.applyOnTheLeft(qr.householderQ().adjoint());
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(timing) timing->qtrhs_transform_seconds+=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-transform_started).count();
        const auto extraction_started=timing ? std::chrono::steady_clock::now() :
            std::chrono::steady_clock::time_point{};
#endif
        const auto keep=std::min(a.rows(),a.cols());
        r=qr.matrixQR().topRows(keep).triangularView<Eigen::Upper>(); target=b.topRows(keep);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(timing) timing->r_extraction_seconds+=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-extraction_started).count();
#endif
    }
};
}
