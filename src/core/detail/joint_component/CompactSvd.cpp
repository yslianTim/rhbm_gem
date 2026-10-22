#include "CompactSvd.hpp"
#include <algorithm>
#include <cmath>

namespace rhbm_gem::core::joint_component {
namespace {
template<class Svd>
CompactSvdResult Decompose(const Matrix & matrix,double relative,double absolute,const Vector * rhs,bool guard,bool right)
{
    CompactSvdResult out; Svd svd;
    svd.compute(matrix);
    if(svd.info()!=Eigen::Success || !svd.singularValues().allFinite()) return out;
    out.singular_values=svd.singularValues();
    const double maximum=out.singular_values(0);
    svd.setThreshold(absolute>=0 && maximum>0 ? absolute/maximum : relative);
    out.threshold=svd.threshold()*maximum;
    if(guard)
    {
        // This only selects the old arithmetic near a rank boundary; it never
        // replaces the caller's threshold or claims an SVD error bound.
        const double margin=64*std::numeric_limits<double>::epsilon()*
            static_cast<double>(std::max(matrix.rows(),matrix.cols()))*maximum;
        if((out.singular_values.array()-out.threshold).abs().minCoeff()<=margin) return out;
    }
    out.rank=svd.rank();
    if(right)
    {
        out.right_vectors=svd.matrixV().leftCols(out.singular_values.size());
        if(!out.right_vectors.allFinite()) return out;
    }
    if(rhs)
    {
        out.solution=svd.solve(*rhs);
        if(!out.solution.allFinite()) return out;
    }
    out.valid=true; return out;
}
template<int Options>
CompactSvdResult Jacobi(const Matrix & a,double relative,double absolute,const Vector * rhs,bool right)
{return Decompose<Eigen::JacobiSVD<Matrix,Options>>(a,relative,absolute,rhs,false,right);}
}
CompactSvdResult CompactSvd(const Matrix & matrix,double relative,double absolute,const Vector * rhs,CompactSvdVectors vectors)
{
    const bool right=vectors==CompactSvdVectors::Right;
    if(matrix.rows()==0 || matrix.cols()==0 || !matrix.allFinite() ||
        (rhs && (rhs->size()!=matrix.rows() || !rhs->allFinite()))) return {};
    const auto jacobi=[&]() {
        if(rhs) return Jacobi<Eigen::ComputeFullU|Eigen::ComputeFullV>(matrix,relative,absolute,rhs,right);
        if(right) return Jacobi<Eigen::ComputeThinV>(matrix,relative,absolute,rhs,right);
        return Jacobi<0>(matrix,relative,absolute,rhs,right);
    };
    CompactSvdResult out;
    const bool bdc=std::min(matrix.rows(),matrix.cols())>=16;
    if(bdc)
    {
        if(rhs) out=Decompose<Eigen::BDCSVD<Matrix,Eigen::ComputeThinU|Eigen::ComputeThinV>>(matrix,relative,absolute,rhs,true,right);
        else if(right) out=Decompose<Eigen::BDCSVD<Matrix,Eigen::ComputeThinV>>(matrix,relative,absolute,rhs,true,right);
        else out=Decompose<Eigen::BDCSVD<Matrix>>(matrix,relative,absolute,rhs,true,right);
        if(!out.valid)
        {
            out=jacobi(); out.jacobi_retry=true;
        }
    }
    else out=jacobi();
    out.used_bdc=bdc;
    return out;
}
CompactSvdResult EvaluateRank(const Matrix & compact,const RankRequest & request,const Vector * rhs,CompactSvdVectors vectors)
{
    if(request.policy.rows<0 || request.columns<=0 || !std::isfinite(request.absolute)) return {};
    auto out=CompactSvd(compact,request.policy.Relative(request.columns),request.absolute,rhs,vectors);
    if(out.valid && request.boundary==RankBoundary::StrictGreater)
        out.rank=(out.singular_values.array()>out.threshold).count();
    return out;
}
bool CompactRankDecisionSafe(const CompactSvdResult & result)
{
    if(!result.valid || result.singular_values.size()==0 || !std::isfinite(result.threshold)) return false;
    const double lower=std::nextafter(result.threshold,0.);
    const double upper=std::nextafter(result.threshold,std::numeric_limits<double>::infinity());
    for(Eigen::Index k=0;k<result.singular_values.size();++k)
        if(result.singular_values(k)>=lower && result.singular_values(k)<=upper) return false;
    return true;
}

}
