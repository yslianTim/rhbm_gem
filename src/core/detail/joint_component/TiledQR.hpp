#pragma once
#include <Eigen/Dense>
#include <algorithm>

namespace rhbm_gem::core::joint_component {
inline constexpr Eigen::Index derivative_tile_rows=8192;
// Consume original rows once. The discarded orthogonal RHS tail is never used
// as an objective: callers retain and evaluate the full residual separately.
struct TiledQR
{
    Eigen::MatrixXd r,target;
    Eigen::MatrixXd scratch_rows,scratch_rhs;
    TiledQR(Eigen::Index columns=0,Eigen::Index responses=0):r(0,columns),target(0,responses) {}
    void Reset(Eigen::Index columns,Eigen::Index responses)
    {r.resize(0,columns); target.resize(0,responses);}
    void Append(const Eigen::MatrixXd & rows,const Eigen::MatrixXd & rhs,bool reference_order=false)
    {
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
        if(reference_order)
        {
            const Eigen::HouseholderQR<Eigen::MatrixXd> qr(a);
            const Eigen::MatrixXd transformed=qr.householderQ().adjoint()*b;
            const auto keep=std::min(a.rows(),a.cols());
            r=qr.matrixQR().topRows(keep).triangularView<Eigen::Upper>();
            target=transformed.topRows(keep);
            return;
        }
        // The assembled tile is disposable. Factor and transform it in place;
        // retain the same Householder arithmetic without two full-size copies.
        const Eigen::HouseholderQR<Eigen::Ref<Eigen::MatrixXd>> qr(a);
        b.applyOnTheLeft(qr.householderQ().adjoint());
        const auto keep=std::min(a.rows(),a.cols());
        r=qr.matrixQR().topRows(keep).triangularView<Eigen::Upper>(); target=b.topRows(keep);
    }
};
}
