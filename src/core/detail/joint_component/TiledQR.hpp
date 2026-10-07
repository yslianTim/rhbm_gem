#pragma once
#include <Eigen/Dense>
#include "ResourceWork.hpp"
#include <algorithm>
#include <cstddef>
#include <string>
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include <chrono>
#endif

namespace rhbm_gem::core::joint_component {
inline constexpr Eigen::Index derivative_tile_rows=8192;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
struct TiledQrTelemetry
{
    std::string role;
    std::size_t append_calls{},rows_processed{},maximum_dense_design_bytes{},maximum_dense_response_bytes{};
    Eigen::Index columns{},responses{},maximum_assembled_rows{};
    double qr_seconds{},assembly_copy_seconds{},householder_seconds{},rhs_transform_seconds{};
};
#endif
// Consume original rows once. The discarded orthogonal RHS tail is never used
// as an objective: callers retain and evaluate the full residual separately.
struct TiledQR
{
    Eigen::MatrixXd r,target;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    TiledQrTelemetry telemetry;
#endif
    TiledQR(Eigen::Index columns,Eigen::Index responses,const char * role=nullptr):r(0,columns),target(0,responses)
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(role) telemetry.role=role;
#else
        (void)role;
#endif
    }
    void Append(const Eigen::MatrixXd & rows,const Eigen::MatrixXd & rhs,bool reference_order=false)
    {
        const auto prior=r.rows();
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        const auto started=std::chrono::steady_clock::now();
        const auto assembled=prior+rows.rows();
        const auto record=[&] {
            ++telemetry.append_calls; telemetry.rows_processed+=static_cast<std::size_t>(rows.rows());
            telemetry.columns=r.cols(); telemetry.responses=target.cols();
            telemetry.maximum_assembled_rows=std::max(telemetry.maximum_assembled_rows,assembled);
            telemetry.maximum_dense_design_bytes=std::max(telemetry.maximum_dense_design_bytes,
                static_cast<std::size_t>(assembled)*static_cast<std::size_t>(r.cols())*sizeof(double));
            telemetry.maximum_dense_response_bytes=std::max(telemetry.maximum_dense_response_bytes,
                static_cast<std::size_t>(assembled)*static_cast<std::size_t>(target.cols())*sizeof(double));
            telemetry.qr_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        };
#endif
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        const auto assembly_started=std::chrono::steady_clock::now();
#endif
        Eigen::MatrixXd a(prior+rows.rows(),r.cols()),b(prior+rows.rows(),target.cols());
        RecordDenseShape("tiled-qr-design",a.rows(),a.cols());
        RecordDenseShape("tiled-qr-response",b.rows(),b.cols());
        a.topRows(prior)=r; a.bottomRows(rows.rows())=rows;
        b.topRows(prior)=target; b.bottomRows(rows.rows())=rhs;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        telemetry.assembly_copy_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-assembly_started).count();
#endif
        if(reference_order)
        {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            const auto householder_started=std::chrono::steady_clock::now();
#endif
            const Eigen::HouseholderQR<Eigen::MatrixXd> qr(a);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            telemetry.householder_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-householder_started).count();
            const auto rhs_started=std::chrono::steady_clock::now();
#endif
            const Eigen::MatrixXd transformed=qr.householderQ().adjoint()*b;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            telemetry.rhs_transform_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-rhs_started).count();
#endif
            const auto keep=std::min(a.rows(),a.cols());
            r=qr.matrixQR().topRows(keep).triangularView<Eigen::Upper>();
            target=transformed.topRows(keep);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            record();
#endif
            return;
        }
        // The assembled tile is disposable. Factor and transform it in place;
        // retain the same Householder arithmetic without two full-size copies.
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        const auto householder_started=std::chrono::steady_clock::now();
#endif
        const Eigen::HouseholderQR<Eigen::Ref<Eigen::MatrixXd>> qr(a);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        telemetry.householder_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-householder_started).count();
        const auto rhs_started=std::chrono::steady_clock::now();
#endif
        b.applyOnTheLeft(qr.householderQ().adjoint());
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        telemetry.rhs_transform_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-rhs_started).count();
#endif
        const auto keep=std::min(a.rows(),a.cols());
        r=qr.matrixQR().topRows(keep).triangularView<Eigen::Upper>(); target=b.topRows(keep);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        record();
#endif
    }
};
}
