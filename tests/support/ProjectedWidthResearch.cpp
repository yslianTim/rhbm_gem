#include "ProjectedWidthResearch.hpp"
#include <Eigen/SVD>
#include <cmath>
#include <limits>

namespace second_stage_test::projected_width {
namespace {
Matrix ConnectedDesign()
{
    Matrix z=Matrix::Zero(8,4);
    z(0,0)=1; z(1,0)=1; z(1,1)=1; z(2,1)=1;
    z(2,2)=1; z(3,2)=1; z(3,3)=1; z(4,3)=1;
    return z;
}
Matrix IndependentWidths()
{
    Matrix d=Matrix::Zero(8,3);
    d(5,0)=1; d(6,1)=1; d(7,2)=1;
    return d;
}
Fixture Make(std::string name,Matrix z,Matrix d,Eigen::Index full_columns=0,
    Eigen::Index targets=0,bool active=false,Eigen::Index policy_rows=-1)
{
    const auto rows=policy_rows<0 ? z.rows() : policy_rows;
    Fixture out; out.name=std::move(name); out.free_design=std::move(z); out.width_derivative=std::move(d);
    out.full_design_columns=full_columns ? full_columns : out.free_design.cols();
    out.target_columns=targets ? targets : out.free_design.cols();
    out.active_amplitude_boundary=active;
    out.policy={rows,out.full_design_columns,out.width_derivative.cols()};
    return out;
}
bool NearThreshold(const Vector & values,double threshold)
{
    const double lower=std::nextafter(threshold,0.0);
    const double upper=std::nextafter(threshold,std::numeric_limits<double>::infinity());
    for(Eigen::Index k=0;k<values.size();++k)
        if(values(k)>=lower && values(k)<=upper) return true;
    return false;
}
SpectrumOracle Spectrum(const Vector & values,double threshold,Eigen::Index columns)
{
    SpectrumOracle out; out.singular_values=values; out.threshold=threshold;
    if(values.size()==0) return out;
    out.rank=(values.array()>threshold).count();
    out.minimum=values(values.size()-1); out.full_rank=out.rank==columns;
    out.threshold_sensitive=NearThreshold(values,threshold);
    return out;
}
}

std::vector<Fixture> Fixtures()
{
    std::vector<Fixture> out;
    const Matrix z=ConnectedDesign();
    out.push_back(Make("well-conditioned-full-rank",Matrix::Identity(8,4),IndependentWidths()));

    Matrix ill=Matrix::Zero(8,3); ill(0,0)=1; ill(1,1)=1; ill(2,2)=1;
    Matrix ill_d=Matrix::Zero(8,3); ill_d(3,0)=1; ill_d(4,1)=1; ill_d(5,2)=1e-7;
    out.push_back(Make("ill-conditioned-full-rank",ill,ill_d));

    Matrix deficient_d=Matrix::Zero(8,3);
    deficient_d(5,0)=1; deficient_d(6,1)=1; deficient_d(5,2)=1; deficient_d(6,2)=1;
    out.push_back(Make("exactly-rank-deficient",z,deficient_d));

    Matrix deficient_z=Matrix::Zero(8,4);
    deficient_z(0,0)=1; deficient_z(1,1)=1; deficient_z(2,2)=1;
    deficient_z(0,3)=1; deficient_z(1,3)=1;
    out.push_back(Make("rank-deficient-free-design",deficient_z,IndependentWidths()));

    Matrix threshold_z=Matrix::Zero(8,3); threshold_z(0,0)=1; threshold_z(1,1)=1; threshold_z(2,2)=1;
    Matrix threshold_d=Matrix::Zero(8,3);
    threshold_d(3,0)=1; threshold_d(4,1)=.25;
    threshold_d(5,2)=std::numeric_limits<double>::epsilon()*1000.0;
    out.push_back(Make("near-rank-threshold",threshold_z,threshold_d,3,3,false,1000));

    out.push_back(Make("connected-full-rank",z,IndependentWidths()));
    Matrix local_d=IndependentWidths(); local_d(5,0)=1; local_d(6,0)=.05; local_d(6,1)=1;
    out.push_back(Make("locally-coupled",z,local_d));

    Matrix strong=Matrix::Zero(8,4);
    strong.topLeftCorner(4,4)<<1,.8,.7,.6, .8,1,.9,.7, .7,.9,1,.8, .6,.7,.8,1;
    Matrix strong_d=Matrix::Zero(8,3);
    strong_d.bottomRightCorner(4,3)<<1,.9,.8, .8,1,.9, .9,.8,1, .7,.9,.8;
    out.push_back(Make("strongly-coupled",strong,strong_d));

    Matrix separator=Matrix::Zero(8,4);
    separator(0,0)=1; separator(1,0)=1; separator(1,1)=1; separator(2,1)=1;
    separator(3,1)=1; separator(3,2)=1; separator(4,2)=1; separator(5,2)=1;
    separator(5,3)=1; separator(6,3)=1;
    out.push_back(Make("weak-separator",separator,IndependentWidths()));

    out.push_back(Make("mixed-target-halo",z,IndependentWidths(),4,3));
    out.push_back(Make("active-amplitude-boundary",z.leftCols(3),IndependentWidths(),4,3,true));

    Matrix multiple=Matrix::Zero(8,4); multiple(0,0)=1; multiple(1,1)=1; multiple(4,2)=1; multiple(5,3)=1;
    Matrix multiple_d=Matrix::Zero(8,2); multiple_d(2,0)=1; multiple_d(6,1)=1;
    out.push_back(Make("multiple-components",multiple,multiple_d));

    Matrix cancellation_z=Matrix::Zero(8,4);
    for(Eigen::Index k=0;k<4;++k) cancellation_z(k,k)=1;
    Matrix cancellation_d=Matrix::Zero(8,3);
    cancellation_d(0,0)=1; cancellation_d(4,0)=1e-13;
    cancellation_d(1,1)=1; cancellation_d(5,1)=1e-13;
    cancellation_d(2,2)=1; cancellation_d(6,2)=1e-13;
    out.push_back(Make("cancellation-sensitive",cancellation_z,cancellation_d));
    return out;
}

DenseOracle EvaluateDenseOracle(const Fixture & fixture)
{
    DenseOracle out; out.rows=fixture.free_design.rows(); out.free_columns=fixture.free_design.cols();
    out.width_columns=fixture.width_derivative.cols();
    if(out.rows<=0 || out.free_columns<=0 || out.width_columns<=0 ||
        fixture.width_derivative.rows()!=out.rows)
    {out.reason="invalid-shape"; return out;}
    Matrix z=fixture.free_design; out.design_column_norms.resize(z.cols());
    for(Eigen::Index k=0;k<z.cols();++k)
    {
        out.design_column_norms(k)=z.col(k).norm();
        if(!(out.design_column_norms(k)>0) || !std::isfinite(out.design_column_norms(k)))
        {out.reason="zero-or-nonfinite-design-column"; return out;}
        z.col(k)/=out.design_column_norms(k);
    }
    const runtime::Sparse design_sparse=z.sparseView(); out.structural_nonzeros=design_sparse.nonZeros();
    for(Eigen::Index row=0;row<z.rows();++row) if(z.row(row).squaredNorm()>0) ++out.supported_rows;

    Eigen::JacobiSVD<Matrix> design_svd(z,Eigen::ComputeThinU);
    if(design_svd.info()!=Eigen::Success || !design_svd.singularValues().allFinite())
    {out.reason="design-svd-failed"; return out;}
    const auto design_values=design_svd.singularValues();
    const double design_tau=fixture.policy.Relative(z.cols())*design_values(0);
    out.design=Spectrum(design_values,design_tau,z.cols());
    const Matrix q=design_svd.matrixU().leftCols(out.design.rank);
    out.projected_width=fixture.width_derivative-q*(q.transpose()*fixture.width_derivative);
    if(!out.projected_width.allFinite()) {out.reason="nonfinite-projection"; return out;}
    Eigen::JacobiSVD<Matrix> width_svd(out.projected_width);
    if(width_svd.info()!=Eigen::Success || !width_svd.singularValues().allFinite())
    {out.reason="projected-width-svd-failed"; return out;}
    const auto width_values=width_svd.singularValues();
    const double width_tau=fixture.policy.Relative(out.width_columns)*width_values(0);
    out.width=Spectrum(width_values,width_tau,out.width_columns);
    out.projection_residual=(z.transpose()*out.projected_width).norm();
    out.valid=true; return out;
}
}
