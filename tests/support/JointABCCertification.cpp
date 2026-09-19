#include "support/JointABCCertification.hpp"
#include "support/JointABCPrecision.hpp"
#include "support/FixedBOracle.hpp"
#include "core/command/detail/SimulationGeometry.hpp"
#include "core/command/detail/SimulationManifest.hpp"
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <sstream>
#include <map>
#include <sys/resource.h>

namespace second_stage_test::matched::certification {
namespace {
namespace j=boost::json;
namespace fs=std::filesystem;
namespace sim=rhbm_gem::core::simulation;
using Vector=Eigen::VectorXd;
using Matrix=Eigen::MatrixXd;
j::value Number(double v) {return std::isfinite(v) ? j::value(v) : j::value(nullptr);}
Vector Parse(const j::value & v) {Vector a(static_cast<Eigen::Index>(v.as_array().size())); for(Eigen::Index k=0;k<a.size();++k) a(k)=j::value_to<double>(v.at(static_cast<std::size_t>(k))); return a;}
j::value Read(const fs::path & path) {std::ifstream in(path); if(!in) throw std::runtime_error("Missing certification file: "+path.string()); j::parse_options o; o.numbers=j::number_precision::precise; return j::parse(std::string(std::istreambuf_iterator<char>(in),{}),{},o);}
void Write(const fs::path & path,const j::value & v) {std::ofstream out(path); out.exceptions(std::ios::failbit|std::ios::badbit); out<<j::serialize(v)<<'\n';}
double Relative(const Vector & a,const Vector & b) {return (a-b).norm()/std::max({1e-12,a.norm(),b.norm()});}
j::object EndpointEvidence(const joint_abc::Evaluation & e)
{
    return {{"valid",e.valid},{"reason",e.reason},{"certificate",e.certificate}};
}
std::vector<std::vector<double>> Table(const fs::path & p)
{
    std::ifstream in(p); if(!in) throw std::runtime_error("Missing snapshot table.");
    std::string line; std::getline(in,line); std::vector<std::vector<double>> rows;
    while(std::getline(in,line)) {std::istringstream s(line); std::string cell; auto & row=rows.emplace_back(); while(std::getline(s,cell,',')) row.push_back(std::stod(cell));}
    return rows;
}
}

joint_abc::Domain Snapshot(const unique_grid::Grid & grid,const std::vector<Atom> & atoms,const fs::path & output)
{
    std::vector<std::vector<joint_abc::Support>> support(atoms.size()); j::array offsets{0};
    std::ofstream table(output/"contributors.csv"); table<<"row,atom,square\n"<<std::setprecision(17); std::size_t entries{};
    for (std::size_t row=0;row<grid.voxels.size();++row)
    {
        for (std::size_t a=0;a<atoms.size();++a)
        {
            const double square=sim::SupportSquare(grid.voxels[row].position,atoms[a].position);
            if(square<=6.25) {support[a].push_back({static_cast<Eigen::Index>(row),square}); table<<row<<','<<a<<','<<square<<'\n'; ++entries;}
        }
        offsets.push_back(entries);
    }
    table.close();
    Write(output/"snapshot.json",j::object{{"schema_version",1},{"support_policy","sphere-fma-v1"},
        {"rows",grid.voxels.size()},{"atoms",atoms.size()},{"memberships",entries},{"row_offsets",offsets},
        {"voxels_sha256",sim::FileSha256(output/"voxels.csv")},{"contributors_sha256",sim::FileSha256(output/"contributors.csv")}});
    return {static_cast<Eigen::Index>(grid.voxels.size()),std::move(support)};
}

void Audit(const joint_abc::Domain & domain,const Vector & y,const j::object & fit,const fs::path & output,const joint_abc::EvaluationContext * provided)
{
    const auto plan=joint_abc::RegisteredAudit(static_cast<Eigen::Index>(domain.atoms.size()),j::value_to<std::string>(fit.at("dataset")),j::value_to<std::string>(fit.at("case")));
    const auto fallback=provided ? joint_abc::EvaluationContext{} : joint_abc::MakeContext(y,static_cast<Eigen::Index>(domain.atoms.size()),"",&plan);
    const auto * context=provided ? provided : &fallback;
    const auto start=std::chrono::steady_clock::now();
    j::object audit{{"schema_version",2},{"case",fit.at("case")},{"dataset",fit.at("dataset")},
        {"legacy_joint_qualified",fit.at("joint_qualified")},{"derivative_status","unavailable"},
        {"derivative_verified",false},{"reason","invalid-endpoint"}};
    if(!fit.contains("primary")) {audit["reason"]="initialization-failed"; Write(output,audit); return;}
    const auto & primary=fit.at("primary");
    if (!primary.at("valid").as_bool()) {Write(output,audit); return;}
    const bool supplied=fit.if_contains("assembled_state_preserved") && fit.at("assembled_state_preserved").as_bool();
    const auto e=supplied ? joint_abc::AtState(domain,y,Parse(primary.at("eta")),Parse(primary.at("beta")),*context) :
        joint_abc::Evaluate(domain,y,Parse(primary.at("eta")),false,context);
    audit["endpoint_trust"]=joint_abc::Trust(domain,y,e,context);
    const double scale=context->scale;
    Vector norms(e.beta.size()); for(Eigen::Index k=0;k<norms.size();++k) norms(k)=e.x.col(k).norm();
    const Vector u=norms.array()*e.beta.array()/scale,dual=(e.x.transpose()*e.residual).array()/norms.array()/scale;
    j::array slack;
    for (Eigen::Index k=0;k<e.beta.size();k+=2) slack.push_back(j::object{{"atom",k/2},{"scaled_a",Number(u(k))},
        {"dual",Number(dual(k))},{"complementarity",Number(u(k)*dual(k))},{"exactly_active",e.beta(k)==0}});
    audit["constraint_evidence"]=slack;
    if (!fit.contains("width_spectrum")) {audit["reason"]="missing-width-spectrum"; Write(output,audit); return;}
    const auto d=joint_abc::Differentiate(e,scale,context);
    if(!d.valid) {audit["reason"]=d.reason; Write(output,audit); return;}
    Matrix directions(e.eta.size(),3); directions.col(0)=Vector::Ones(e.eta.size()).normalized();
    for(Eigen::Index k=0;k<e.eta.size();++k) directions(k,1)=k%2 ? -1 : 1;
    directions.col(1).normalize(); directions.col(2)=Parse(fit.at("width_spectrum").at("weak_directions").at(0));
    if(context->audit.directions.size()) directions=context->audit.directions;
    const bool precision=context->audit.precision;
    if(precision) audit["precision"]=PrecisionAudit(domain,y,e,directions,context,context->audit.block_precision);
    const bool old_verified=fit.at("derivative_verified").as_bool();
    bool verified=old_verified; bool crossed{},missing_interval{};
    j::array ladders;
    const bool expanded_audit=context->audit.expanded_if_unverified && (!old_verified || precision);
    if(expanded_audit)
    {
        verified=true;
        for(Eigen::Index direction=0;direction<directions.cols();++direction)
        {
            std::vector<Vector> finite; std::vector<double> noise; std::vector<bool> valid;
            j::array samples;
            for(int k=0;k<=16;++k)
            {
                const double h=std::ldexp(.01,-k);
                const auto plus=joint_abc::Evaluate(domain,y,e.eta+h*directions.col(direction),false,context);
                const auto minus=joint_abc::Evaluate(domain,y,e.eta-h*directions.col(direction),false,context);
                const auto rp=joint_abc::Evaluate(domain,y,e.eta+h*directions.col(direction),true,context);
                const auto rm=joint_abc::Evaluate(domain,y,e.eta-h*directions.col(direction),true,context);
                const bool solved=plus.valid && minus.valid && rp.valid && rm.valid;
                const bool face=solved && plus.certificate.at("active_atoms")==e.certificate.at("active_atoms") &&
                    minus.certificate.at("active_atoms")==e.certificate.at("active_atoms") &&
                    rp.certificate.at("active_atoms")==e.certificate.at("active_atoms") && rm.certificate.at("active_atoms")==e.certificate.at("active_atoms");
                crossed|=solved && !face; valid.push_back(face);
                finite.push_back(solved ? Vector((rp.residual-rm.residual)/(2*h*scale)) : Vector::Zero(y.size()));
                noise.push_back(solved ? ((plus.residual-rp.residual).norm()+(minus.residual-rm.residual).norm())/(2*h*scale) : std::numeric_limits<double>::infinity());
                const Vector analytic=d.jacobian*directions.col(direction);
                samples.push_back(j::object{{"h",h},{"same_face",face},{"plus",EndpointEvidence(rp)},{"minus",EndpointEvidence(rm)},
                    {"primary_plus",EndpointEvidence(plus)},{"primary_minus",EndpointEvidence(minus)},
                    {"plus_coefficient_difference",solved ? Number(((plus.beta-rp.beta).array().abs()/(1+plus.beta.array().abs().max(rp.beta.array().abs()))).maxCoeff()) : j::value(nullptr)},
                    {"minus_coefficient_difference",solved ? Number(((minus.beta-rm.beta).array().abs()/(1+minus.beta.array().abs().max(rm.beta.array().abs()))).maxCoeff()) : j::value(nullptr)},
                    {"cancellation_ratio",solved ? Number((rp.residual.norm()+rm.residual.norm())/std::max(std::numeric_limits<double>::min(),(rp.residual-rm.residual).norm())) : j::value(nullptr)},
                    {"absolute_error",solved ? Number((finite.back()-analytic).norm()) : j::value(nullptr)},
                    {"relative_error",solved ? Number(Relative(finite.back(),analytic)) : j::value(nullptr)},
                    {"inner_noise",Number(noise.back())},{"sensitivity",Number(analytic.norm())}});
            }
            // Selection uses only successive finite differences and inner disagreement.
            double best=std::numeric_limits<double>::infinity(); int selected=-1; Vector first,second;
            j::array candidates;
            for(int k=0;k<15;++k)
            {
                if(!valid[static_cast<std::size_t>(k)] || !valid[static_cast<std::size_t>(k+1)] || !valid[static_cast<std::size_t>(k+2)]) continue;
                const auto i=static_cast<std::size_t>(k);
                const Vector a=(4*finite[i+1]-finite[i])/3,b=(4*finite[i+2]-finite[i+1])/3;
                const double uncertainty=((a-b).norm()+(noise[i]+5*noise[i+1]+4*noise[i+2])/3)/std::max({1e-12,a.norm(),b.norm()});
                candidates.push_back(j::object{{"index",k},{"estimated_relative_error",Number(uncertainty)}});
                if(uncertainty<best) {best=uncertainty; selected=k; first=a; second=b;}
            }
            const Vector analytic=d.jacobian*directions.col(direction);
            const bool pass=selected>=0 && best<=1e-7 && Relative(first,analytic)<=1e-6 && Relative(second,analytic)<=1e-6;
            missing_interval|=selected<0;
            verified &= pass;
            ladders.push_back(j::object{{"direction",direction},{"samples",samples},{"candidates",candidates},
                {"selected",selected<0 ? j::value(nullptr) : j::value(selected)},{"estimated_relative_error",Number(best)},
                {"first_relative_error",selected>=0 ? Number(Relative(first,analytic)) : j::value(nullptr)},
                {"second_relative_error",selected>=0 ? Number(Relative(second,analytic)) : j::value(nullptr)},{"passed",pass}});
        }
        if(!precision || !audit.at("precision").at("agreement_passed").as_bool() ||
            !audit.at("precision").at("derivative_passed").as_bool()) verified=false;
    }
    audit["ladders"]=ladders; audit["derivative_verified"]=verified;
    audit["derivative_status"]=verified ? "verified" : crossed && missing_interval ? "transition-unverified" : "resolution-unverified";
    audit["reason"]=verified ? (expanded_audit ? "richardson-and-high-precision" : "legacy-two-step-evidence") : "insufficient-derivative-evidence";
    audit["seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    struct rusage usage{}; getrusage(RUSAGE_SELF,&usage);
#ifdef __APPLE__
    audit["process_peak_rss_bytes"]=usage.ru_maxrss;
#else
    audit["process_peak_rss_bytes"]=usage.ru_maxrss*1024;
#endif
    Write(output,audit);
}

void AuditDirectory(const fs::path & directory)
{
    Eigen::setNbThreads(1);
    std::vector<fs::path> datasets;
    for(const auto & entry:fs::directory_iterator(directory/"datasets")) datasets.push_back(entry.path());
    std::sort(datasets.begin(),datasets.end());
    for(const auto & dataset:datasets)
    {
        const auto root=dataset; const auto snapshot=Read(root/"snapshot.json");
        if(sim::FileSha256(root/"voxels.csv")!=j::value_to<std::string>(snapshot.at("voxels_sha256")) ||
            sim::FileSha256(root/"contributors.csv")!=j::value_to<std::string>(snapshot.at("contributors_sha256"))) throw std::runtime_error("Snapshot hash mismatch.");
        const auto voxels=Table(root/"voxels.csv"),contributors=Table(root/"contributors.csv");
        if(j::value_to<std::size_t>(snapshot.at("rows"))!=voxels.size() ||
            j::value_to<std::size_t>(snapshot.at("memberships"))!=contributors.size()) throw std::runtime_error("Invalid snapshot population.");
        std::vector<std::vector<joint_abc::Support>> support(j::value_to<std::size_t>(snapshot.at("atoms")));
        for(const auto & row:contributors) support.at(static_cast<std::size_t>(row.at(1))).push_back({static_cast<Eigen::Index>(row.at(0)),row.at(2)});
        const joint_abc::Domain domain(static_cast<Eigen::Index>(voxels.size()),std::move(support));
        Vector y64(domain.rows),y32(domain.rows);
        for(Eigen::Index k=0;k<domain.rows;++k) {y64(k)=voxels[static_cast<std::size_t>(k)][7]; y32(k)=voxels[static_cast<std::size_t>(k)][8];}
        std::map<std::string,std::pair<j::value,std::string>> audits,boundaries;
        for(const std::string variant:{"legacy","guarded","guarded-log"})
        {
            std::vector<fs::path> files;
            for(const auto & file:fs::directory_iterator(root/variant/"fits")) files.push_back(file.path());
            std::sort(files.begin(),files.end());
            for(const auto & file:files)
            {
                const auto fit=Read(file).as_object(); const auto & y=file.stem().string().ends_with("double") ? y64 : y32;
                if(j::value_to<std::string>(fit.at("observation_snapshot_sha256"))!=sim::FileSha256(root/"snapshot.json"))
                    throw std::runtime_error("Modified observation snapshot.");
                std::cout<<"Audit "<<root.filename()<<'/'<<variant<<'/'<<file.stem()<<std::endl;
                const auto output=root/variant/"audits"/file.filename();
                const std::string key=j::serialize(j::array{file.stem().string().ends_with("double"),
                    fit.if_contains("primary") ? fit.at("primary") : j::value(nullptr),fit.at("joint_qualified"),fit.if_contains("width_spectrum") ? fit.at("width_spectrum") : j::value(nullptr)});
                const std::string source=variant+"/audits/"+file.filename().string();
                if(audits.contains(key))
                {
                    auto value=audits.at(key).first.as_object(); value["case"]=fit.at("case");
                    value["reused_identical_endpoint"]=audits.at(key).second; value["seconds"]=0;
                    if(value.contains("precision")) value.at("precision").as_object()["seconds"]=0;
                    Write(output,value);
                }
                else {Audit(domain,y,fit,output); audits.emplace(key,std::make_pair(Read(output),source));}
                if(root.filename()=="active-a" && fit.contains("primary") && fit.at("primary").at("valid").as_bool())
                {
                    auto boundary=boundaries.contains(key) ? boundaries.at(key).first.as_object() :
                        BoundaryAudit(domain,y,Parse(fit.at("primary").at("eta")),Parse(fit.at("primary").at("beta")));
                    if(boundaries.contains(key)) {boundary["seconds"]=0; boundary["reused_identical_endpoint"]=boundaries.at(key).second;}
                    else boundaries.emplace(key,std::make_pair(boundary,source));
                    Write(root/variant/"audits"/(file.stem().string()+"-boundary.json"),boundary);
                }
            }
        }
        if(root.filename()=="active-a")
        {
            const auto truth=Read(root/"scoring-truth.json");
            for(bool quantized:{false,true}) Write(root/(quantized ? "truth-boundary-float32.json" : "truth-boundary-double.json"),
                BoundaryAudit(domain,quantized ? y32 : y64,Parse(truth.at("b")).array().log(),Parse(truth.at("beta"))));
        }
    }
}
} // namespace second_stage_test::matched::certification
