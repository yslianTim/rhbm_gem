#pragma once
#include <Eigen/Dense>

namespace rhbm_gem::core::joint_component {
inline constexpr Eigen::Index derivative_tile_rows=8192;
// Consume original rows once. The discarded orthogonal RHS tail is never used
// as an objective: callers retain and evaluate the full residual separately.
struct TiledQR
{
    Eigen::MatrixXd r,target;
    TiledQR(Eigen::Index columns,Eigen::Index responses):r(0,columns),target(0,responses) {}
    void Append(const Eigen::MatrixXd & rows,const Eigen::MatrixXd & rhs)
    {
        const auto prior=r.rows();
        Eigen::MatrixXd a(prior+rows.rows(),r.cols()),b(prior+rows.rows(),target.cols());
        a.topRows(prior)=r; a.bottomRows(rows.rows())=rows;
        b.topRows(prior)=target; b.bottomRows(rows.rows())=rhs;
        // The assembled tile is disposable. Factor and transform it in place;
        // retain the same Householder arithmetic without two full-size copies.
        const Eigen::HouseholderQR<Eigen::Ref<Eigen::MatrixXd>> qr(a);
        b.applyOnTheLeft(qr.householderQ().adjoint());
        const auto keep=std::min(a.rows(),a.cols());
        r=qr.matrixQR().topRows(keep).triangularView<Eigen::Upper>(); target=b.topRows(keep);
    }
};
}
