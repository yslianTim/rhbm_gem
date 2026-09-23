#include "CompactSvd.hpp"
#include "SparseFactor.hpp"
#include <algorithm>

namespace rhbm_gem::core::joint_component {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
CompactSvdMode & CompactSvdModeForTesting() {static thread_local auto mode=CompactSvdMode::Automatic; return mode;}
CompactSvdCapture & CompactSvdCaptureForTesting() {static thread_local CompactSvdCapture capture; return capture;}
#endif
namespace {
template<class Svd>
CompactSvdResult Decompose(const Matrix & matrix,double relative,double absolute,const Vector * rhs,bool guard,bool right)
{
    auto & work=SparseWorkForTesting();
    CompactSvdResult out; Svd svd;
    {
        WorkTimer timer(rhs ? work.reference_svd_seconds : work.free_design_svd_seconds);
        svd.compute(matrix);
    }
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
        ++work.reference_solves; WorkTimer timer(work.reference_solve_seconds);
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
    auto & work=SparseWorkForTesting();
    ++(rhs ? work.reference_svds : work.free_design_svds);
    if(matrix.rows()==0 || matrix.cols()==0 || !matrix.allFinite() ||
        (rhs && (rhs->size()!=matrix.rows() || !rhs->allFinite()))) return {};
    bool automatic=true,legacy=false;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    automatic=CompactSvdModeForTesting()==CompactSvdMode::Automatic;
    legacy=CompactSvdModeForTesting()==CompactSvdMode::Legacy;
#endif
    const auto jacobi=[&]() {
        if(rhs) return Jacobi<Eigen::ComputeFullU|Eigen::ComputeFullV>(matrix,relative,absolute,rhs,right);
        if(legacy) return Jacobi<Eigen::ComputeThinU|Eigen::ComputeThinV>(matrix,relative,absolute,rhs,right);
        if(right) return Jacobi<Eigen::ComputeThinV>(matrix,relative,absolute,rhs,right);
        return Jacobi<0>(matrix,relative,absolute,rhs,right);
    };
    CompactSvdResult out;
    const bool bdc=automatic && std::min(matrix.rows(),matrix.cols())>=16;
    if(bdc)
    {
        ++work.bdc_svds;
        if(rhs) out=Decompose<Eigen::BDCSVD<Matrix,Eigen::ComputeThinU|Eigen::ComputeThinV>>(matrix,relative,absolute,rhs,true,right);
        else if(right) out=Decompose<Eigen::BDCSVD<Matrix,Eigen::ComputeThinV>>(matrix,relative,absolute,rhs,true,right);
        else out=Decompose<Eigen::BDCSVD<Matrix>>(matrix,relative,absolute,rhs,true,right);
        if(!out.valid)
        {
            ++work.jacobi_retries; WorkTimer timer(work.jacobi_retry_seconds);
            out=jacobi(); out.jacobi_retry=true;
        }
    }
    else out=jacobi();
    out.used_bdc=bdc;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    if(CompactSvdCaptureForTesting()) CompactSvdCaptureForTesting()(matrix,relative,absolute,rhs,out);
#endif
    return out;
}
}
