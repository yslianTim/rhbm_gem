#include "support/JointTestNumerics.hpp"
#include "support/JointRuntimeJson.hpp"
#include <chrono>

namespace second_stage_test::matched::joint_abc {
namespace {
namespace j=boost::json;
using Matrix=Eigen::MatrixXd;
using Vector=Eigen::VectorXd;
double Seconds(std::chrono::steady_clock::time_point start)
{return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();}

}

Domain::Domain(const unique_grid::Grid & grid,const std::vector<Atom> & geometry)
    :runtime::Domain(static_cast<Eigen::Index>(grid.voxels.size()),std::vector<std::vector<Support>>(geometry.size()))
{
    std::vector<std::vector<Support>> support(geometry.size());
    for (std::size_t a=0;a<geometry.size();++a) for (Eigen::Index p=0;p<rows;++p)
    {
        const double square=SquareDistance(grid.voxels[static_cast<std::size_t>(p)].position,geometry[a].position);
        if (square<=2.5*2.5) support[a].push_back({p,square});
    }
    static_cast<runtime::Domain &>(*this)=runtime::Domain(rows,std::move(support));
}

Evaluation::Evaluation(runtime::Evaluation e):runtime::Evaluation(std::move(e)),certificate(runtime_json::Certificate(runtime::Evaluation::certificate)) {}
Evaluation Evaluate(const Domain & d,const Vector & y,const Vector & eta,bool reference,const EvaluationContext * c,const std::vector<runtime::LinearBlock> * blocks)
{return Evaluation(runtime::EvaluateProfile(d,y,eta,reference,c,blocks));}
Evaluation AtState(const Domain & d,const Vector & y,const Vector & eta,const Vector & beta,const EvaluationContext & c)
{return Evaluation(runtime::EvaluateState(d,y,eta,beta,c));}
j::object MatrixSpectrum(const Sparse & x,const RankPolicy & p,Eigen::Index columns,bool normalize)
{return runtime_json::MatrixSpectrum(runtime::ComputeSpectrum(x,p,columns,normalize));}
j::object MatrixSpectrum(const Matrix & x,const RankPolicy & p,Eigen::Index columns,bool normalize)
{return runtime_json::MatrixSpectrum(runtime::ComputeSpectrum(x,p,columns,normalize));}
Vector LocalCorrection(const Evaluation & e,const Differential & d,const EvaluationContext & c,double threshold)
{return runtime::ComputeLocalCorrection(e,d,c,threshold);}
Differential Differentiate(const Evaluation & e,double scale,const EvaluationContext * c,double threshold)
{return runtime::DifferentiateProfile(e,scale,c,threshold);}
j::object Trust(const Domain & d,const Vector & y,const Evaluation & e,const EvaluationContext * c)
{
    const auto fallback=c ? EvaluationContext{} : MakeContext(y,static_cast<Eigen::Index>(d.atoms.size()));
    return runtime_json::Trust(runtime::CheckTrust(d,y,e,c ? *c : fallback));
}
j::object Assess(const Domain & d,const Vector & y,const Vector & eta,const EvaluationContext & c,const Vector * beta)
{return runtime_json::Assessment(runtime::AssessProfile(d,y,eta,c,beta));}
j::object Fit(const Domain & domain,const Vector & y,const Vector & initial_b,const EvaluationContext * provided)
{
    const auto fallback=provided ? EvaluationContext{} : MakeContext(y,static_cast<Eigen::Index>(domain.atoms.size()));
    const auto & context=provided ? *provided : fallback;
    const auto start=std::chrono::steady_clock::now();
    const auto search=runtime::SearchProfile(domain,y,initial_b,context);
    auto out=runtime_json::Search(search,context,y.size());
    auto assessment=Assess(domain,y,search.eta,context);
    for(auto & field:assessment) out[field.key()]=std::move(field.value());
    out["seconds"]=Seconds(start); return out;
}

} // namespace second_stage_test::matched::joint_abc
