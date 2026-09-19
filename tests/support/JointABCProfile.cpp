#include "support/JointABCProfile.hpp"
#include "support/FixedBOracle.hpp"
#include "support/JointRuntimeJson.hpp"
#include "core/command/detail/SimulationManifest.hpp"
#include <unsupported/Eigen/NonLinearOptimization>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sys/resource.h>

namespace second_stage_test::matched::joint_abc {
namespace {
namespace j=boost::json;
namespace fs=std::filesystem;
using Matrix=Eigen::MatrixXd;
using Vector=Eigen::VectorXd;
j::value Number(double x) {return std::isfinite(x) ? j::value(x) : j::value(nullptr);}
j::array Values(const Vector & v) {j::array out; for (double x:v) out.push_back(Number(x)); return out;}
Vector Parse(const j::value & v)
{
    Vector out(static_cast<Eigen::Index>(v.as_array().size()));
    for (Eigen::Index k=0;k<out.size();++k) out(k)=j::value_to<double>(v.at(static_cast<std::size_t>(k)));
    return out;
}
double Seconds(std::chrono::steady_clock::time_point start)
{return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();}
void Write(const fs::path & path,const j::value & value)
{
    std::ofstream out(path); out.exceptions(std::ios::failbit|std::ios::badbit); out<<j::serialize(value)<<'\n';
}
std::ofstream CSV(const fs::path & path,const std::string & header)
{
    std::ofstream out(path); out.exceptions(std::ios::failbit|std::ios::badbit);
    out<<header<<'\n'<<std::setprecision(17); return out;
}
}

Domain::Domain(const unique_grid::Grid & grid,const std::vector<Atom> & geometry)
    :runtime::Domain(static_cast<Eigen::Index>(grid.voxels.size()),std::vector<std::vector<Support>>(geometry.size()))
{
    for (std::size_t a=0;a<geometry.size();++a) for (Eigen::Index p=0;p<rows;++p)
    {
        const double square=SquareDistance(grid.voxels[static_cast<std::size_t>(p)].position,geometry[a].position);
        if (square<=2.5*2.5) atoms[a].push_back({p,square});
    }
}

Evaluation::Evaluation(runtime::Evaluation e):runtime::Evaluation(std::move(e)),certificate(runtime_json::Certificate(runtime::Evaluation::certificate)) {}
Evaluation Evaluate(const Domain & d,const Vector & y,const Vector & eta,bool reference,const EvaluationContext * c,const std::vector<joint_ac::LinearBlock> * blocks)
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
j::object Fit(const Domain & domain,const Vector & y,const Vector & initial_b,j::object * resources,const std::string & variant,const EvaluationContext * provided)
{
    const auto fallback=provided ? EvaluationContext{} : MakeContext(y,static_cast<Eigen::Index>(domain.atoms.size()));
    const auto & context=provided ? *provided : fallback;
    const auto start=std::chrono::steady_clock::now();
    auto measure=[](auto since) {
        struct rusage usage{}; getrusage(RUSAGE_SELF,&usage);
#ifdef __APPLE__
        const auto bytes=usage.ru_maxrss;
#else
        const auto bytes=usage.ru_maxrss*1024;
#endif
        return j::object{{"seconds",Seconds(since)},{"process_peak_rss_bytes",bytes}};
    };
    const auto kind=variant.empty() ? runtime::SearchVariant::Original : variant=="guarded" ? runtime::SearchVariant::Guarded :
        variant=="guarded-log" ? runtime::SearchVariant::GuardedLog : runtime::SearchVariant::Legacy;
    const auto search=runtime::SearchProfile(domain,y,initial_b,context,kind);
    if(resources) (*resources)["search"]=measure(start);
    const auto audit_start=std::chrono::steady_clock::now();
    auto out=runtime_json::Search(search,context,variant,y.size());
    auto assessment=Assess(domain,y,search.eta,context);
    for(auto & field:assessment) out[field.key()]=std::move(field.value());
    if(resources) (*resources)["endpoint_audit"]=measure(audit_start);
    out["seconds"]=Seconds(start); return out;
}

void Run(const std::string & manifest,const std::string & map,const std::string & checkpoint,const std::string & output_path)
{
    const auto data=fixed_b::Prepare(manifest,map,checkpoint,output_path,"joint-abc-profile");
    const Domain domain(data.grid,data.atoms); const fs::path output(output_path);
    fs::create_directories(output/"weak-directions"); j::array completed;
    for (bool quantized:{false,true}) for (const std::string start:{"checkpoint","narrower","wider","mixed"})
    {
        const std::string name=start+(quantized ? "-float32" : "-double");
        Vector initial=data.checkpoint_b;
        for (Eigen::Index a=0;a<initial.size();++a)
        {
            const int serial=j::value_to<int>(data.identities.at(static_cast<std::size_t>(a)).at("serial_id"));
            if (start=="narrower" || (start=="mixed" && serial%2)) initial(a)*=.8;
            else if (start=="wider" || start=="mixed") initial(a)*=1.2;
        }
        const auto & y=quantized ? data.y32 : data.y64;
        std::cout<<"Starting "<<name<<std::endl;
        auto fit=Fit(domain,y,initial); fit["case"]=name; fit["initial_b"]=Values(initial);
        struct rusage usage{}; getrusage(RUSAGE_SELF,&usage);
#ifdef __APPLE__
        fit["process_peak_rss_bytes"]=usage.ru_maxrss;
#else
        fit["process_peak_rss_bytes"]=usage.ru_maxrss*1024;
#endif
        const auto & endpoint=fit.at("primary");
        if (endpoint.at("valid").as_bool())
        {
            auto model=data.atoms; const auto beta=Parse(endpoint.at("beta")),b=Parse(endpoint.at("b"));
            // Output replay is not an optimizer evaluation or an endpoint certificate.
            for (std::size_t a=0;a<model.size();++a)
            {model[a].amplitude=beta(static_cast<Eigen::Index>(2*a)); model[a].charge=beta(static_cast<Eigen::Index>(2*a+1)); model[a].width=b(static_cast<Eigen::Index>(a));}
            auto residual=CSV(output/"residuals"/(name+".csv"),"row,prediction,residual");
            for (Eigen::Index p=0;p<y.size();++p)
            {
                // Match sparse column accumulation order on the original rows.
                double prediction{};
                for (const auto & atom:model)
                {
                    const auto basis=EvaluateBasis(SquareDistance(data.grid.voxels[static_cast<std::size_t>(p)].position,atom.position),atom.width,2.5);
                    prediction+=atom.amplitude*basis.gaussian; prediction+=atom.charge*basis.charge;
                }
                residual<<p<<','<<prediction<<','<<prediction-y(p)<<'\n';
            }
            residual.close();
            fit["residual_sha256"]=rhbm_gem::core::simulation::FileSha256((output/"residuals"/(name+".csv")).string());
        }
        if (fit.contains("width_spectrum"))
        {
            auto weak=CSV(output/"weak-directions"/(name+".csv"),"serial_id,weakest,second,third");
            const auto & directions=fit.at("width_spectrum").at("weak_directions");
            for (std::size_t a=0;a<data.atoms.size();++a)
            {
                weak<<j::value_to<int>(data.identities[a].at("serial_id"));
                for (std::size_t k=0;k<3;++k) weak<<','<<j::value_to<double>(directions.at(k).at(a));
                weak<<'\n';
            }
        }
        Write(output/"fits"/(name+".json"),fit); completed.emplace_back(name);
        Write(output/"completion.json",j::object{{"complete",completed.size()==8},{"execution_complete",completed.size()==8},{"cases",completed}});
        std::cout<<name<<" qualified="<<fit.at("joint_qualified")<<" seconds="<<fit.at("seconds")<<std::endl;
    }
}
} // namespace second_stage_test::matched::joint_abc
