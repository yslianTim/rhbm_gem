#include "MatrixFreeProjectedWidth.hpp"
#include <Eigen/SVD>
#include <algorithm>
#include <new>

namespace second_stage_test::projected_width {
namespace {
bool Add(std::size_t a,std::size_t b,std::size_t & out)
{return runtime::CheckedAddBytes(a,b,out);}
bool Product(std::size_t a,std::size_t b,std::size_t & out)
{return runtime::CheckedMultiplyBytes(a,b,out);}
bool AddProduct(std::size_t & total,std::size_t count,std::size_t width)
{
    std::size_t bytes{}; return Product(count,width,bytes) && Add(total,bytes,total);
}
bool SparseStorage(const runtime::Sparse & x,std::size_t & bytes)
{
    if(x.nonZeros()<0 || x.outerSize()<0) return false;
    return Product(static_cast<std::size_t>(x.nonZeros()),sizeof(double)+sizeof(Eigen::Index),bytes) &&
        AddProduct(bytes,static_cast<std::size_t>(x.outerSize())+1,sizeof(Eigen::Index));
}
}

Vector OperatorResult::Project(const Eigen::Ref<const Vector> & u) const
{
    if(!available || u.size()!=normalized_design.rows()) return {};
    ++projection_calls;
    return basis*(basis.transpose()*u);
}
Vector OperatorResult::Apply(const Eigen::Ref<const Vector> & v) const
{
    if(!available || v.size()!=width_columns) return {};
    ++forward_calls;
    const Vector raw=derivative*v;
    return raw-Project(raw);
}
Vector OperatorResult::ApplyAdjoint(const Eigen::Ref<const Vector> & u) const
{
    if(!available || u.size()!=normalized_design.rows()) return {};
    ++adjoint_calls;
    const Vector residual=u-Project(u);
    return derivative.transpose()*residual;
}
double OperatorResult::OrthogonalityResidual(const Eigen::Ref<const Vector> & u) const
{
    if(!available || u.size()!=normalized_design.rows()) return std::numeric_limits<double>::quiet_NaN();
    return (normalized_design.transpose()*u).norm();
}

OperatorResult BuildExactOracleOperator(const runtime::Sparse & design,const runtime::Sparse & derivative,
    const runtime::RankPolicy & policy,std::size_t resident_bytes,std::size_t budget)
{
    OperatorResult out;
    if(design.rows()<=0 || design.cols()<=0 || derivative.rows()!=design.rows() || derivative.cols()<=0 ||
        policy.rows<0)
    {out.reason="invalid-shape"; return out;}
    const auto rows=static_cast<std::size_t>(design.rows()),columns=static_cast<std::size_t>(design.cols());
    const auto widths=static_cast<std::size_t>(derivative.cols());
    std::size_t input_design{},input_derivative{},dense_workspace{},vector_workspace{},estimate{},resident{},vectors{};
    if(!SparseStorage(design,input_design) || !SparseStorage(derivative,input_derivative) ||
        !Product(rows,columns,dense_workspace) || !Product(dense_workspace,sizeof(double),dense_workspace) ||
        !Product(dense_workspace,5,dense_workspace) ||
        !Add(rows,columns,vectors) || !Add(vectors,widths,vectors) ||
        !Product(vectors,8*sizeof(double),vector_workspace) ||
        !Add(input_design,input_derivative,resident) ||
        !Add(resident_bytes ? resident_bytes : resident,input_design,estimate) ||
        !Add(estimate,input_derivative,estimate) || !Add(estimate,dense_workspace,estimate) ||
        !Add(estimate,vector_workspace,estimate))
    {out.status="ResourceLimit"; out.reason="memory-estimate-overflow"; return out;}
    out.estimated_peak_bytes=estimate;
    if(estimate>budget)
    {out.status="ResourceLimit"; out.reason="estimated-peak-exceeds-planning-budget"; return out;}

    try
    {
        out.normalized_design=design;
        out.column_norms.resize(design.cols());
        for(Eigen::Index k=0;k<design.cols();++k)
        {
            const double norm=design.col(k).norm();
            if(!(norm>0) || !std::isfinite(norm))
            {out.reason="zero-or-nonfinite-free-design-column"; return out;}
            out.column_norms(k)=norm;
            for(runtime::Sparse::InnerIterator entry(out.normalized_design,k);entry;++entry) entry.valueRef()/=norm;
        }
        Matrix dense(out.normalized_design);
        Eigen::JacobiSVD<Matrix> svd(dense,Eigen::ComputeThinU);
        if(svd.info()!=Eigen::Success || svd.singularValues().size()==0 || !svd.singularValues().allFinite())
        {out.reason="dense-oracle-projection-svd-failed"; return out;}
        out.design_threshold=policy.Relative(design.cols())*svd.singularValues()(0);
        out.design_rank=(svd.singularValues().array()>out.design_threshold).count();
        out.basis=svd.matrixU().leftCols(out.design_rank);
        out.derivative=derivative;
        out.width_columns=derivative.cols();
        out.basis_bytes=static_cast<std::size_t>(out.basis.size())*sizeof(double);
        out.available=true; out.status="Completed"; out.reason="exact-dense-oracle-projection";
    }
    catch(const std::bad_alloc &)
    {out.status="ResourceLimit"; out.reason="projection-allocation-failed";}
    return out;
}

BoundedProjectionPlan PreflightBoundedProjection(std::size_t rows,std::size_t free_columns,
    std::size_t width_columns,std::size_t design_nonzeros,std::size_t derivative_nonzeros,
    std::size_t resident_bytes,std::size_t budget)
{
    BoundedProjectionPlan out; out.budget_bytes=budget;
    std::size_t sparse{},vectors{},peak{},columns{},outer{};
    if(!AddProduct(sparse,design_nonzeros,sizeof(double)+sizeof(Eigen::Index)) ||
        !AddProduct(sparse,derivative_nonzeros,sizeof(double)+sizeof(Eigen::Index)) ||
        !Add(free_columns,width_columns,columns) || !Add(columns,2,outer) ||
        !AddProduct(sparse,outer,sizeof(Eigen::Index)) ||
        !AddProduct(vectors,rows,3*sizeof(double)) ||
        !AddProduct(vectors,columns,4*sizeof(double)) ||
        !Add(resident_bytes,sparse,peak) || !Add(peak,vectors,peak))
    {out.reason="memory-estimate-overflow"; return out;}
    out.estimated_action_bytes=vectors; out.estimated_peak_bytes=peak;
    if(peak>budget)
    {
        out.resource_status="ResourceLimit"; out.reason="streaming-action-storage-exceeds-planning-budget";
        return out;
    }
    out.resource_status="Pass";
    out.numerical_status="Unavailable";
    out.reason="no-bounded-global-projector-with-validated-rank-contract";
    return out;
}
}
