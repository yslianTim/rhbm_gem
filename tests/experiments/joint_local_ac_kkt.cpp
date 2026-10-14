#include "core/detail/joint_component/Numerics.hpp"
#include <boost/json.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>

namespace {
namespace j=boost::json;
namespace n=rhbm_gem::core::joint_component;

j::value Read(const char * path)
{
    std::ifstream input(path);
    if(!input) throw std::runtime_error(std::string("Cannot read local fixture: ")+path);
    j::parse_options options; options.numbers=j::number_precision::precise;
    return j::parse(std::string(std::istreambuf_iterator<char>(input),{}),{},options);
}
n::Vector VectorFrom(const j::value & value)
{
    const auto & array=value.as_array(); n::Vector out(static_cast<Eigen::Index>(array.size()));
    for(std::size_t k=0;k<array.size();++k) out(static_cast<Eigen::Index>(k))=j::value_to<double>(array[k]);
    return out;
}
std::vector<std::string> StringsFrom(const j::value & value)
{
    std::vector<std::string> out;
    for(const auto & item:value.as_array()) out.push_back(j::value_to<std::string>(item));
    return out;
}
n::Domain DomainFrom(const j::value & fixture)
{
    const auto rows=j::value_to<Eigen::Index>(fixture.at("identity").at("local_design_rows"));
    std::vector<std::vector<n::Support>> support;
    for(const auto & atom:fixture.at("support").as_array())
    {
        std::vector<n::Support> entries;
        for(const auto & item:atom.as_array()) entries.push_back({
            j::value_to<Eigen::Index>(item.at("row")),j::value_to<double>(item.at("squared_distance"))});
        support.push_back(std::move(entries));
    }
    return {rows,std::move(support)};
}
n::EvaluationContext ContextFrom(const j::value & fixture)
{
    const auto & value=fixture.at("numerical_context"); n::EvaluationContext out;
    out.scale=j::value_to<double>(value.at("parent_observation_scale"));
    out.linear.rank_relative=j::value_to<double>(value.at("rank_relative"));
    out.linear.release_factor=j::value_to<double>(value.at("release_factor"));
    out.linear.release_response_norm=j::value_to<double>(value.at("release_response_norm"));
    out.linear.active_set_iteration_factor=j::value_to<int>(value.at("active_set_iteration_factor"));
    out.atom_ids=n::Identities(StringsFrom(fixture.at("atom_ids")));
    out.row_ids=n::Identities(StringsFrom(fixture.at("row_ids")));
    out.independent_search=true; out.update_budget=1;
    return out;
}
j::array NumbersJson(n::VectorRef values);
using MP=boost::multiprecision::cpp_dec_float_50;
MP ExactInput(double value)
{
    if(value==0.0) return 0;
    int exponent{}; const double fraction=std::frexp(std::abs(value),&exponent);
    const auto mantissa=static_cast<std::uint64_t>(std::ldexp(fraction,std::numeric_limits<double>::digits));
    MP exact=boost::multiprecision::ldexp(MP(mantissa),exponent-std::numeric_limits<double>::digits);
    if(std::signbit(value)) exact=-exact;
    return exact;
}
MP HighPrecisionKkt(const n::Evaluation & evaluation,n::VectorRef y,const std::vector<MP> & beta,
    const std::vector<Eigen::Index> & free_columns,double scale)
{
    std::vector<MP> residual(static_cast<std::size_t>(y.size()));
    for(Eigen::Index row=0;row<y.size();++row) residual[static_cast<std::size_t>(row)]=-ExactInput(y(row));
    for(Eigen::Index column=0;column<evaluation.x.cols();++column)
        for(n::Sparse::InnerIterator entry(evaluation.x,column);entry;++entry)
            residual[static_cast<std::size_t>(entry.row())]+=ExactInput(entry.value())*beta[static_cast<std::size_t>(column)];
    std::vector<bool> is_free(static_cast<std::size_t>(evaluation.x.cols()));
    for(const auto column:free_columns) is_free[static_cast<std::size_t>(column)]=true;
    MP kkt{};
    for(Eigen::Index column=0;column<evaluation.x.cols();++column)
    {
        MP norm_squared{},gradient{};
        for(n::Sparse::InnerIterator entry(evaluation.x,column);entry;++entry)
        {
            const MP value=ExactInput(entry.value()); norm_squared+=value*value;
            gradient+=value*residual[static_cast<std::size_t>(entry.row())];
        }
        const MP norm=sqrt(norm_squared),u=norm*beta[static_cast<std::size_t>(column)]/MP(scale);
        const MP scaled_gradient=gradient/norm/MP(scale);
        const MP unprojected=u-scaled_gradient;
        const MP projected=is_free[static_cast<std::size_t>(column)] || unprojected>0 ? unprojected : MP(0);
        const MP error=abs(u-projected);
        if(error>kkt) kkt=error;
    }
    return kkt;
}
MP HighPrecisionObjective(const n::Evaluation & evaluation,n::VectorRef y,const std::vector<MP> & beta)
{
    std::vector<MP> residual(static_cast<std::size_t>(y.size()));
    for(Eigen::Index row=0;row<y.size();++row) residual[static_cast<std::size_t>(row)]=-ExactInput(y(row));
    for(Eigen::Index column=0;column<evaluation.x.cols();++column)
        for(n::Sparse::InnerIterator entry(evaluation.x,column);entry;++entry)
            residual[static_cast<std::size_t>(entry.row())]+=ExactInput(entry.value())*beta[static_cast<std::size_t>(column)];
    MP objective{};
    for(const auto & value:residual) objective+=value*value;
    return objective/2;
}
j::value HighPrecisionReference(const j::value & fixture)
{
    const auto domain=DomainFrom(fixture); const auto y=VectorFrom(fixture.at("effective_response"));
    const auto eta=VectorFrom(fixture.at("old_eta")); const auto context=ContextFrom(fixture);
    const auto evaluation=n::EvaluateProfile(domain,y,eta,false,&context);
    if(evaluation.x.rows()!=y.size() || evaluation.x.cols()!=2*eta.size() ||
        evaluation.beta.size()!=evaluation.x.cols() || !evaluation.beta.allFinite())
        throw std::runtime_error("High-precision reference requires a finite profile design and coefficients.");
    std::vector<Eigen::Index> free_columns;
    for(Eigen::Index column=0;column<evaluation.x.cols();++column)
        if(column%2!=0 || evaluation.beta(column)>0) free_columns.push_back(column);
    constexpr Eigen::Index maximum_dense_elements=1'000'000;
    if(evaluation.x.rows()>maximum_dense_elements/static_cast<Eigen::Index>(free_columns.size()))
        throw std::runtime_error("High-precision reference exceeds the bounded local-face budget.");
    const auto count=static_cast<Eigen::Index>(free_columns.size());
    std::vector<MP> norms(static_cast<std::size_t>(count));
    std::vector<std::vector<MP>> normal(static_cast<std::size_t>(count),
        std::vector<MP>(static_cast<std::size_t>(count)+1));
    for(Eigen::Index column=0;column<count;++column)
    {
        norms[static_cast<std::size_t>(column)]=ExactInput(evaluation.x.col(
            free_columns[static_cast<std::size_t>(column)]).norm());
    }
    for(Eigen::Index row=0;row<evaluation.x.rows();++row)
    {
        std::vector<MP> values(static_cast<std::size_t>(count));
        for(Eigen::Index column=0;column<count;++column)
            values[static_cast<std::size_t>(column)]=ExactInput(evaluation.x.coeff(row,
                free_columns[static_cast<std::size_t>(column)]))/norms[static_cast<std::size_t>(column)];
        for(Eigen::Index left=0;left<count;++left)
        {
            auto & equation=normal[static_cast<std::size_t>(left)];
            equation[static_cast<std::size_t>(count)]+=values[static_cast<std::size_t>(left)]*ExactInput(y(row));
            for(Eigen::Index right=0;right<count;++right)
                equation[static_cast<std::size_t>(right)]+=values[static_cast<std::size_t>(left)]*
                    values[static_cast<std::size_t>(right)];
        }
    }
    for(Eigen::Index pivot=0;pivot<count;++pivot)
    {
        Eigen::Index selected=pivot;
        for(Eigen::Index row=pivot+1;row<count;++row)
            if(abs(normal[static_cast<std::size_t>(row)][static_cast<std::size_t>(pivot)])>
               abs(normal[static_cast<std::size_t>(selected)][static_cast<std::size_t>(pivot)])) selected=row;
        if(normal[static_cast<std::size_t>(selected)][static_cast<std::size_t>(pivot)]==0)
            throw std::runtime_error("High-precision free face is rank deficient.");
        std::swap(normal[static_cast<std::size_t>(pivot)],normal[static_cast<std::size_t>(selected)]);
        const auto diagonal=normal[static_cast<std::size_t>(pivot)][static_cast<std::size_t>(pivot)];
        for(Eigen::Index row=pivot+1;row<count;++row)
        {
            auto & equation=normal[static_cast<std::size_t>(row)];
            const MP factor=equation[static_cast<std::size_t>(pivot)]/diagonal;
            for(Eigen::Index column=pivot+1;column<=count;++column)
                equation[static_cast<std::size_t>(column)]-=factor*
                    normal[static_cast<std::size_t>(pivot)][static_cast<std::size_t>(column)];
        }
    }
    std::vector<MP> theta(static_cast<std::size_t>(count));
    for(Eigen::Index row=count;row-->0;)
    {
        MP value=normal[static_cast<std::size_t>(row)][static_cast<std::size_t>(count)];
        for(Eigen::Index column=row+1;column<count;++column)
            value-=normal[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)]*
                theta[static_cast<std::size_t>(column)];
        theta[static_cast<std::size_t>(row)]=value/normal[static_cast<std::size_t>(row)][static_cast<std::size_t>(row)];
    }
    std::vector<MP> reference_beta(static_cast<std::size_t>(evaluation.x.cols()));
    for(Eigen::Index column=0;column<count;++column)
        reference_beta[static_cast<std::size_t>(free_columns[static_cast<std::size_t>(column)])]=
            theta[static_cast<std::size_t>(column)]/norms[static_cast<std::size_t>(column)];
    std::vector<MP> rounded_beta(static_cast<std::size_t>(evaluation.x.cols()));
    for(Eigen::Index column=0;column<evaluation.beta.size();++column)
        rounded_beta[static_cast<std::size_t>(column)]=ExactInput(evaluation.beta(column));
    auto rounded_reference_beta=reference_beta;
    n::Vector rounded_beta_vector=n::Vector::Zero(evaluation.x.cols());
    for(std::size_t column=0;column<rounded_reference_beta.size();++column)
    {
        const double value=rounded_reference_beta[column].convert_to<double>();
        rounded_beta_vector(static_cast<Eigen::Index>(column))=value;
        rounded_reference_beta[column]=ExactInput(value);
    }
    const auto rounded_certificate=n::CertifyLinear(evaluation.x,y,rounded_beta_vector,context.scale);
    const auto text=[](const MP & value) {return value.str(45,std::ios_base::scientific);};
    j::array rounded_beta_json=NumbersJson(rounded_beta_vector);
    j::array free_json; for(const auto column:free_columns) free_json.push_back(column);
    return j::object{{"precision","boost::multiprecision::cpp_dec_float_50"},
        {"free_columns",std::move(free_json)},
        {"normal_equations_conditioning_note","50-digit normal equations for fixed double design; 23-column bounded local face"},
        {"production_beta_exact_design_kkt",text(HighPrecisionKkt(evaluation,y,rounded_beta,free_columns,context.scale))},
        {"high_precision_solution_kkt",text(HighPrecisionKkt(evaluation,y,reference_beta,free_columns,context.scale))},
        {"rounded_high_precision_solution_kkt",text(HighPrecisionKkt(evaluation,y,rounded_reference_beta,free_columns,context.scale))},
        {"production_beta_exact_design_objective",text(HighPrecisionObjective(evaluation,y,rounded_beta))},
        {"high_precision_solution_objective",text(HighPrecisionObjective(evaluation,y,reference_beta))},
        {"rounded_high_precision_solution_objective",text(HighPrecisionObjective(evaluation,y,rounded_reference_beta))},
        {"rounded_high_precision_solution_double_certificate_kkt",rounded_certificate.projected_kkt},
        {"rounded_high_precision_solution_double_certificate_objective",rounded_certificate.objective},
        {"rounded_high_precision_beta",std::move(rounded_beta_json)}};
}
j::value Run(const j::value & fixture)
{
    if(j::value_to<std::string>(fixture.at("format"))!="rhbm-joint-local-ac-kkt-v1")
        throw std::invalid_argument("Unsupported local fixture format.");
    const auto domain=DomainFrom(fixture); const auto y=VectorFrom(fixture.at("effective_response"));
    const auto eta=VectorFrom(fixture.at("old_eta")); const auto context=ContextFrom(fixture);
    const auto evaluation=n::EvaluateProfile(domain,y,eta,false,&context);
    const auto search=n::SearchProfile(domain,y,eta.array().exp(),context);
    const auto spectrum=n::ComputeSpectrum(evaluation.x,{domain.rows,evaluation.x.cols(),eta.size()},evaluation.x.cols(),true);
    const auto & original=fixture.at("production_solver");
    const auto captured_beta=VectorFrom(original.at("beta"));
    const auto captured_certificate=n::CertifyLinear(evaluation.x,y,captured_beta,context.scale);
    std::vector<Eigen::Index> original_active;
    for(const auto & atom:original.at("active_a_atom_indices").as_array())
        original_active.push_back(j::value_to<Eigen::Index>(atom));
    const bool active_face_matches=evaluation.certificate.active_atoms==original_active;
    const bool recovered=evaluation.valid && evaluation.certificate.linear_solver_valid &&
        evaluation.certificate.available && evaluation.certificate.feasible && evaluation.certificate.kkt_passed &&
        evaluation.certificate.free_rank.value_or(-1)==j::value_to<int>(original.at("linear_rank")) &&
        evaluation.certificate.linear_releases.value_or(-1)==0 && active_face_matches && captured_certificate.available &&
        evaluation.certificate.objective<=captured_certificate.objective;
    double best_neighbor_kkt=evaluation.certificate.projected_kkt;
    Eigen::Index best_neighbor_column=-1; int best_neighbor_direction=0;
    for(Eigen::Index column=0;column<evaluation.beta.size();++column)
    {
        if(column%2==0 && evaluation.beta(column)==0) continue;
        for(const double direction:{-std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity()})
        {
            auto neighbor=evaluation.beta;
            neighbor(column)=std::nextafter(neighbor(column),direction);
            const auto certificate=n::CertifyLinear(evaluation.x,y,neighbor,context.scale);
            if(certificate.available && certificate.feasible && certificate.projected_kkt<best_neighbor_kkt)
            {
                best_neighbor_kkt=certificate.projected_kkt; best_neighbor_column=column;
                best_neighbor_direction=std::signbit(direction) ? -1 : 1;
            }
        }
    }
    return j::object{
        {"recovered",recovered},{"same_active_face",active_face_matches},
        {"search_stop_reason",search.stop_reason},{"initial_profile_reason",evaluation.reason},
        {"affected_row_count",fixture.at("identity").at("affected_row_count")},
        {"local_design_rows",evaluation.x.rows()},{"local_design_columns",evaluation.x.cols()},
        {"local_design_nonzero_entries",evaluation.x.nonZeros()},
        {"local_design_rank",spectrum.available ? j::value(spectrum.rank) : j::value(nullptr)},
        {"condition_indicator",spectrum.available ? j::value(spectrum.condition) : j::value(nullptr)},
        {"free_rank",evaluation.certificate.free_rank ? j::value(*evaluation.certificate.free_rank) : j::value(nullptr)},
        {"linear_solve_count",evaluation.certificate.linear_solves ? j::value(*evaluation.certificate.linear_solves) : j::value(nullptr)},
        {"active_set_releases",evaluation.certificate.linear_releases ? j::value(*evaluation.certificate.linear_releases) : j::value(nullptr)},
        {"scaled_projected_kkt",evaluation.certificate.projected_kkt},
        {"result_objective",evaluation.certificate.objective},
        {"captured_recorded_objective",j::value_to<double>(original.at("objective"))},
        {"captured_beta_current_objective",captured_certificate.objective},
        {"result_beta",NumbersJson(evaluation.beta)},
        {"best_single_ulp_neighbor_kkt",best_neighbor_kkt},
        {"best_single_ulp_neighbor_column",best_neighbor_column},
        {"best_single_ulp_neighbor_direction",best_neighbor_direction},
        {"primal_feasible",evaluation.certificate.feasible},
        {"effective_response_norm",y.norm()},{"parent_observation_scale",context.scale}};
}
j::value Number(double value) {return std::isfinite(value) ? j::value(value) : j::value(nullptr);}
j::array IndicesJson(const std::vector<Eigen::Index> & values)
{j::array out; for(const auto value:values) out.push_back(value); return out;}
j::array NumbersJson(n::VectorRef values)
{j::array out; for(Eigen::Index k=0;k<values.size();++k) out.push_back(Number(values(k))); return out;}
struct LinearSummary
{
    n::LinearResult linear;
    n::Certificate certificate;
    n::Spectrum free_spectrum;
    double free_raw_gradient{},free_scaled_gradient{},active_raw_violation{},active_scaled_violation{},
        active_release_gradient{},release_tolerance{};
};
LinearSummary Summarize(const n::Sparse & x,n::VectorRef y,const n::LinearPolicy & policy,double scale,
    n::LinearResult linear)
{
    LinearSummary out; out.linear=std::move(linear);
    out.certificate=n::CertifyLinear(x,y,out.linear.beta,scale);
    if(out.linear.beta.size()!=x.cols()) return out;
    n::Vector norms(x.cols());
    for(Eigen::Index k=0;k<x.cols();++k) norms(k)=x.col(k).norm();
    std::vector<Eigen::Index> free_columns;
    for(Eigen::Index column=0;column<x.cols();++column)
        if(column%2!=0 || out.linear.beta(column)>0) free_columns.push_back(column);
    n::Sparse free_design(x.rows(),static_cast<Eigen::Index>(free_columns.size()));
    std::vector<Eigen::Triplet<double>> entries;
    const n::Vector residual=x*out.linear.beta-y;
    const n::Vector raw_gradient=x.transpose()*residual;
    for(std::size_t free_index=0;free_index<free_columns.size();++free_index)
    {
        const auto column=free_columns[free_index];
        for(n::Sparse::InnerIterator entry(x,column);entry;++entry)
            entries.emplace_back(entry.row(),static_cast<Eigen::Index>(free_index),entry.value());
        out.free_raw_gradient=std::max(out.free_raw_gradient,std::abs(raw_gradient(column)));
        out.free_scaled_gradient=std::max(out.free_scaled_gradient,
            std::abs(raw_gradient(column))/norms(column)/scale);
    }
    for(Eigen::Index column=0;column<x.cols();column+=2) if(out.linear.beta(column)==0)
    {
        out.active_raw_violation=std::max(out.active_raw_violation,std::max(0.0,-raw_gradient(column)));
        const double scaled_violation=std::max(0.0,-raw_gradient(column))/norms(column)/scale;
        out.active_scaled_violation=std::max(out.active_scaled_violation,scaled_violation);
        out.active_release_gradient=std::max(out.active_release_gradient,
            std::max(0.0,-raw_gradient(column))/norms(column));
    }
    free_design.setFromTriplets(entries.begin(),entries.end());
    out.free_spectrum=n::ComputeSpectrum(free_design,{x.rows(),free_design.cols(),x.cols()/2},free_design.cols(),true);
    n::Sparse scaled=x;
    for(int column=0;column<scaled.outerSize();++column)
        for(n::Sparse::InnerIterator entry(scaled,column);entry;++entry) entry.valueRef()/=norms(column);
    const double znorm=scaled.norm(); const n::Vector scaled_beta=norms.array()*out.linear.beta.array();
    const double response_norm=policy.release_response_norm>=0 ? policy.release_response_norm : y.norm();
    out.release_tolerance=policy.release_factor*std::numeric_limits<double>::epsilon()*
        std::max(1.0,znorm*(response_norm+znorm*scaled_beta.norm()));
    return out;
}
j::object SummaryJson(const std::string & name,const LinearSummary & summary)
{
    const auto & cert=summary.certificate;
    const std::string status=!summary.linear.valid ? summary.linear.reason :
        !cert.available ? "certificate-unavailable" : cert.kkt_passed ? "qualified" : "kkt-failed";
    return {{"solver",name},{"status",status},{"rank",summary.linear.rank},
        {"objective",Number(cert.objective)},{"scaled_kkt",Number(cert.projected_kkt)},
        {"feasible",cert.feasible},{"active_a_atoms",IndicesJson(cert.active_atoms)},
        {"solves",summary.linear.solves},{"active_set_releases",summary.linear.releases},
        {"release_tolerance_raw",Number(summary.release_tolerance)},
        {"active_release_gradient_raw",Number(summary.active_release_gradient)},
        {"active_a_complementarity_violation_raw",Number(summary.active_raw_violation)},
        {"active_a_complementarity_violation_scaled",Number(summary.active_scaled_violation)},
        {"free_rank",summary.free_spectrum.available ? j::value(summary.free_spectrum.rank) : j::value(nullptr)},
        {"free_condition_indicator",summary.free_spectrum.available ? Number(summary.free_spectrum.condition) : j::value(nullptr)},
        {"free_raw_gradient_inf",Number(summary.free_raw_gradient)},
        {"free_scaled_gradient_inf",Number(summary.free_scaled_gradient)}};
}
struct FaceReference
{
    bool valid{};
    n::Vector beta;
    int rank{};
    std::size_t faces{},full_rank_faces{},feasible_faces{},kkt_feasible_faces{};
    std::uint32_t best_active_mask{},kkt_active_mask{};
    bool best_face_kkt_passed{},has_kkt_face{};
    double kkt_objective{},kkt_scaled{};
};
FaceReference EnumerateFaces(const n::Sparse & x,n::VectorRef y,const n::EvaluationContext & context)
{
    const Eigen::Index columns=x.cols(),atoms=columns/2;
    const n::Matrix dense(x);
    n::Vector norms(columns);
    for(Eigen::Index k=0;k<columns;++k) norms(k)=dense.col(k).norm();
    const n::Matrix scaled=dense*norms.cwiseInverse().asDiagonal();
    const Eigen::HouseholderQR<n::Matrix> reduction(scaled);
    const n::Matrix r=reduction.matrixQR().topRows(columns).triangularView<Eigen::Upper>();
    const n::Vector transformed=(reduction.householderQ().adjoint()*y).head(columns);
    FaceReference best; double best_objective=std::numeric_limits<double>::infinity();
    const std::uint32_t face_count=std::uint32_t{1}<<static_cast<unsigned>(atoms);
    for(std::uint32_t mask=0;mask<face_count;++mask)
    {
        ++best.faces; std::vector<Eigen::Index> free_columns;
        for(Eigen::Index atom=0;atom<atoms;++atom)
        {
            if(!(mask&(std::uint32_t{1}<<static_cast<unsigned>(atom)))) free_columns.push_back(2*atom);
            free_columns.push_back(2*atom+1);
        }
        n::Matrix face_r(columns,static_cast<Eigen::Index>(free_columns.size()));
        for(std::size_t k=0;k<free_columns.size();++k)
            face_r.col(static_cast<Eigen::Index>(k))=r.col(free_columns[k]);
        Eigen::ColPivHouseholderQR<n::Matrix> qr(face_r); qr.setThreshold(context.linear.rank_relative);
        if(qr.rank()!=face_r.cols()) continue;
        ++best.full_rank_faces;
        const n::Vector theta=qr.solve(transformed); n::Vector beta=n::Vector::Zero(columns);
        for(std::size_t k=0;k<free_columns.size();++k)
            beta(free_columns[k])=theta(static_cast<Eigen::Index>(k))/norms(free_columns[k]);
        bool feasible=beta.allFinite();
        for(Eigen::Index atom=0;atom<atoms;++atom) feasible &= beta(2*atom)>=0;
        if(!feasible) continue;
        ++best.feasible_faces;
        const double objective=.5*(dense*beta-y).squaredNorm();
        const auto certificate=n::CertifyLinear(x,y,beta,context.scale);
        if(certificate.kkt_passed)
        {
            ++best.kkt_feasible_faces;
            if(!best.has_kkt_face)
            {
                best.has_kkt_face=true; best.kkt_active_mask=mask;
                best.kkt_objective=certificate.objective; best.kkt_scaled=certificate.projected_kkt;
            }
        }
        if(objective<best_objective)
        {best.valid=true; best_objective=objective; best.beta=std::move(beta); best.rank=static_cast<int>(qr.rank());
            best.best_active_mask=mask; best.best_face_kkt_passed=certificate.kkt_passed;}
    }
    return best;
}
struct FaceSolveReference
{
    n::Vector beta;
    int rank{},precision_bits{};
    bool valid{},feasible{};
};
FaceSolveReference SolveFaceDense(const n::Sparse & x,n::VectorRef y,n::VectorRef face_beta,double rank_relative)
{
    if(face_beta.size()!=x.cols()) return {};
    const Eigen::Index atoms=x.cols()/2; const n::Matrix dense(x);
    std::vector<Eigen::Index> free_columns;
    for(Eigen::Index atom=0;atom<atoms;++atom)
    {
        if(face_beta(2*atom)>0) free_columns.push_back(2*atom);
        free_columns.push_back(2*atom+1);
    }
    n::Matrix design(x.rows(),static_cast<Eigen::Index>(free_columns.size()));
    n::Vector norms(static_cast<Eigen::Index>(free_columns.size()));
    for(std::size_t k=0;k<free_columns.size();++k)
    {
        const auto column=free_columns[k]; norms(static_cast<Eigen::Index>(k))=dense.col(column).norm();
        design.col(static_cast<Eigen::Index>(k))=dense.col(column)/norms(static_cast<Eigen::Index>(k));
    }
    Eigen::ColPivHouseholderQR<n::Matrix> qr(design); qr.setThreshold(rank_relative);
    const n::Vector theta=qr.solve(y); FaceSolveReference out;
    out.beta=n::Vector::Zero(x.cols()); out.rank=static_cast<int>(qr.rank()); out.precision_bits=std::numeric_limits<double>::digits;
    for(std::size_t k=0;k<free_columns.size();++k)
        out.beta(free_columns[k])=theta(static_cast<Eigen::Index>(k))/norms(static_cast<Eigen::Index>(k));
    out.valid=theta.allFinite() && out.rank==static_cast<int>(free_columns.size()) && out.beta.allFinite();
    out.feasible=out.valid;
    for(Eigen::Index atom=0;atom<atoms;++atom) out.feasible &= out.beta(2*atom)>=0;
    return out;
}
FaceSolveReference SolveFaceLongDouble(const n::Sparse & x,n::VectorRef y,n::VectorRef face_beta)
{
    using LongMatrix=Eigen::Matrix<long double,Eigen::Dynamic,Eigen::Dynamic>;
    using LongVector=Eigen::Matrix<long double,Eigen::Dynamic,1>;
    if(face_beta.size()!=x.cols()) return {};
    const n::Matrix dense(x); const Eigen::Index atoms=x.cols()/2;
    std::vector<Eigen::Index> free_columns;
    for(Eigen::Index atom=0;atom<atoms;++atom)
    {
        if(face_beta(2*atom)>0) free_columns.push_back(2*atom);
        free_columns.push_back(2*atom+1);
    }
    LongMatrix scaled(x.rows(),static_cast<Eigen::Index>(free_columns.size()));
    LongVector norms(static_cast<Eigen::Index>(free_columns.size())),response(y.size());
    for(Eigen::Index row=0;row<y.size();++row) response(row)=static_cast<long double>(y(row));
    for(std::size_t k=0;k<free_columns.size();++k)
    {
        const auto column=free_columns[k]; long double squared{};
        for(Eigen::Index row=0;row<x.rows();++row)
        {
            const long double value=static_cast<long double>(dense(row,column));
            scaled(row,static_cast<Eigen::Index>(k))=value; squared+=value*value;
        }
        norms(static_cast<Eigen::Index>(k))=std::sqrt(squared);
        scaled.col(static_cast<Eigen::Index>(k))/=norms(static_cast<Eigen::Index>(k));
    }
    Eigen::ColPivHouseholderQR<LongMatrix> qr(scaled);
    const LongVector theta=qr.solve(response); FaceSolveReference out;
    out.beta=n::Vector::Zero(x.cols()); out.rank=static_cast<int>(qr.rank());
    out.precision_bits=std::numeric_limits<long double>::digits;
    for(std::size_t k=0;k<free_columns.size();++k)
        out.beta(free_columns[k])=static_cast<double>(theta(static_cast<Eigen::Index>(k))/norms(static_cast<Eigen::Index>(k)));
    out.valid=theta.allFinite() && out.rank==static_cast<int>(free_columns.size()) && out.beta.allFinite();
    out.feasible=out.valid;
    for(Eigen::Index atom=0;atom<atoms;++atom) out.feasible &= out.beta(2*atom)>=0;
    return out;
}
 j::array RefineFace(const n::Sparse & x,n::VectorRef y,n::VectorRef face_beta,n::VectorRef initial_beta,
    const n::EvaluationContext & context)
{
    const Eigen::Index atoms=x.cols()/2; const n::Matrix dense(x);
    std::vector<Eigen::Index> free_columns;
    for(Eigen::Index atom=0;atom<atoms;++atom)
    {
        if(face_beta(2*atom)>0) free_columns.push_back(2*atom);
        free_columns.push_back(2*atom+1);
    }
    n::Matrix design(x.rows(),static_cast<Eigen::Index>(free_columns.size()));
    n::Vector norms(static_cast<Eigen::Index>(free_columns.size())),theta(norms.size());
    for(std::size_t k=0;k<free_columns.size();++k)
    {
        const auto column=free_columns[k]; norms(static_cast<Eigen::Index>(k))=dense.col(column).norm();
        design.col(static_cast<Eigen::Index>(k))=dense.col(column)/norms(static_cast<Eigen::Index>(k));
        theta(static_cast<Eigen::Index>(k))=initial_beta(column)*norms(static_cast<Eigen::Index>(k));
    }
    Eigen::ColPivHouseholderQR<n::Matrix> qr(design); qr.setThreshold(context.linear.rank_relative);
    j::array out;
    for(int iteration=0;iteration<=10;++iteration)
    {
        n::Vector beta=n::Vector::Zero(x.cols());
        for(std::size_t k=0;k<free_columns.size();++k)
            beta(free_columns[k])=theta(static_cast<Eigen::Index>(k))/norms(static_cast<Eigen::Index>(k));
        const auto certificate=n::CertifyLinear(x,y,beta,context.scale);
        out.push_back(j::object{{"corrections",iteration},{"objective",Number(certificate.objective)},
            {"scaled_kkt",Number(certificate.projected_kkt)},{"kkt_passed",certificate.kkt_passed},
            {"active_a_atoms",IndicesJson(certificate.active_atoms)},
            {"free_face_rank",static_cast<int>(qr.rank())}});
        if(iteration==10) break;
        const n::Vector correction=qr.solve(y-design*theta);
        if(!correction.allFinite()) break;
        theta+=correction;
    }
    return out;
}
double CoefficientDifference(n::VectorRef a,n::VectorRef b)
{return a.size()==b.size() && a.size() ? (a-b).norm()/(1.0+std::max(a.norm(),b.norm())) : std::numeric_limits<double>::infinity();}
j::value Assess(const j::value & fixture)
{
    const auto domain=DomainFrom(fixture); const auto y=VectorFrom(fixture.at("effective_response"));
    const auto eta=VectorFrom(fixture.at("old_eta")); const auto context=ContextFrom(fixture);
    const auto previous_work=n::ResourceWorkForTesting();
    auto & work=n::ResourceWorkForTesting(); work={}; work.enabled=true;
    const auto production=n::EvaluateProfile(domain,y,eta,false,&context);
    const auto production_work=work; work=previous_work;
    const n::Vector weights=n::Vector::Ones(y.size()); const n::Matrix dense(production.x);
    auto dense_qr=n::SolveLinear(dense,y,weights,false,false,nullptr,&context.linear,context.scale);
    auto dense_svd=n::SolveLinear(dense,y,weights,true,false,nullptr,&context.linear,context.scale);
    const auto qr_summary=Summarize(production.x,y,context.linear,context.scale,std::move(dense_qr));
    const auto svd_summary=Summarize(production.x,y,context.linear,context.scale,std::move(dense_svd));
    const auto faces=EnumerateFaces(production.x,y,context);
    const auto face_summary=Summarize(production.x,y,context.linear,context.scale,
        n::LinearResult{faces.beta,faces.valid,faces.valid ? "face-enumeration" : "no-feasible-face",faces.rank,
            static_cast<int>(faces.faces),0,0});
    const auto dense_face=SolveFaceDense(production.x,y,faces.beta,context.linear.rank_relative);
    const auto dense_face_summary=Summarize(production.x,y,context.linear,context.scale,
        n::LinearResult{dense_face.beta,dense_face.valid,dense_face.valid ? "dense-free-face-qr" : "dense-free-face-rank-deficient",
            dense_face.rank,1,0,0});
    const auto alternate_face=SolveFaceLongDouble(production.x,y,faces.beta);
    const auto alternate_face_summary=Summarize(production.x,y,context.linear,context.scale,
        n::LinearResult{alternate_face.beta,alternate_face.valid,alternate_face.valid ? "long-double-qr" : "long-double-rank-deficient",
            alternate_face.rank,1,0,0});
    n::LinearPolicy parent_policy=context.linear;
    parent_policy.release_response_norm=j::value_to<double>(fixture.at("numerical_context").at("parent_response_norm"));
    n::LinearPolicy local_policy=context.linear; local_policy.release_response_norm=y.norm();
    n::LinearPolicy scaled_policy=context.linear;
    const double parent_scale=context.scale;
    const double normalized_local_response_norm=y.norm()/parent_scale;
    scaled_policy.release_response_norm=parent_scale*normalized_local_response_norm;
    const auto parent_summary=Summarize(production.x,y,parent_policy,parent_scale,
        n::SolveLinear(production.x,y,weights,false,true,nullptr,&parent_policy,nullptr,parent_scale));
    const auto local_summary=Summarize(production.x,y,local_policy,parent_scale,
        n::SolveLinear(production.x,y,weights,false,true,nullptr,&local_policy,nullptr,parent_scale));
    const auto scaled_summary=Summarize(production.x,y,scaled_policy,parent_scale,
        n::SolveLinear(production.x,y,weights,false,true,nullptr,&scaled_policy,nullptr,parent_scale));
    n::LinearResult production_result{production.beta,production.certificate.linear_solver_valid,
        production.certificate.linear_solver_reason,production.certificate.free_rank.value_or(0),
        production.certificate.linear_solves.value_or(0),production.certificate.linear_releases.value_or(0),0};
    const auto production_summary=Summarize(production.x,y,context.linear,context.scale,std::move(production_result));
    const auto & original=fixture.at("production_solver");
    const auto original_beta=VectorFrom(original.at("beta"));
    const auto refinement=RefineFace(production.x,y,faces.beta,original_beta,context);
    const double original_objective=j::value_to<double>(original.at("objective"));
    const bool active_face_matches=production_summary.certificate.active_atoms==face_summary.certificate.active_atoms;
    const bool policy_cases_pass=parent_summary.certificate.kkt_passed && local_summary.certificate.kkt_passed &&
        scaled_summary.certificate.kkt_passed && parent_summary.linear.releases==0 &&
        local_summary.linear.releases==0 && scaled_summary.linear.releases==0 &&
        parent_summary.certificate.active_atoms==face_summary.certificate.active_atoms &&
        local_summary.certificate.active_atoms==face_summary.certificate.active_atoms &&
        scaled_summary.certificate.active_atoms==face_summary.certificate.active_atoms;
    const bool assessment_valid=production_summary.linear.valid && production_summary.certificate.available &&
        production_summary.certificate.feasible && production_summary.certificate.kkt_passed &&
        production_summary.linear.rank==face_summary.linear.rank && production_summary.linear.releases==0 &&
        production_summary.certificate.objective<=original_objective &&
        production_summary.certificate.objective<=face_summary.certificate.objective && active_face_matches &&
        faces.valid && faces.faces==4096 && faces.feasible_faces>0 && policy_cases_pass &&
        refinement.size()>2 && refinement[2].at("kkt_passed").as_bool() &&
        refinement[2].at("objective").as_double()<=original_objective;
    j::array refinement_shapes;
    bool bounded_refinement_observed=false;
    for(const auto & shape:production_work.dense_shapes) if(shape.role=="ac-bounded-kkt-refinement")
    {
        refinement_shapes.push_back(j::object{{"phase",shape.phase},{"rows",shape.rows},{"columns",shape.columns},
            {"observations",shape.observations},{"maximum_bytes",static_cast<std::uint64_t>(shape.maximum_bytes)}});
        bounded_refinement_observed |= shape.rows==production.x.rows() &&
            shape.columns==production_summary.linear.rank && shape.maximum_bytes<=8'000'000;
    }
    const bool final_assessment_valid=assessment_valid && bounded_refinement_observed;
    j::array table{SummaryJson("Production SparseQR",production_summary),SummaryJson("Dense QR",qr_summary),
        SummaryJson("Dense SVD",svd_summary),SummaryJson("Independent active-face enumeration",face_summary),
        SummaryJson("Dense free-face QR",dense_face_summary),SummaryJson("Long-double free-face QR",alternate_face_summary)};
    j::array policy_cases{SummaryJson("Parent response norm",parent_summary),
        SummaryJson("Local effective response norm",local_summary),
        SummaryJson("Consistently scaled local RHS",scaled_summary)};
    const auto spectrum=n::ComputeSpectrum(production.x,{production.x.rows(),production.x.cols(),eta.size()},production.x.cols(),true);
    return j::object{{"fixture",fixture.at("format")},{"dimensions",j::object{
            {"affected_rows",fixture.at("identity").at("affected_row_count")},
            {"local_design_rows",production.x.rows()},{"columns",production.x.cols()},
            {"nonzero_entries",production.x.nonZeros()},{"full_design_rank",fixture.at("identity").at("local_design_rank")},
            {"full_design_condition_indicator",spectrum.available ? Number(spectrum.condition) : j::value(nullptr)}}},
        {"response_scales",j::object{{"parent_norm",parent_policy.release_response_norm},
            {"local_effective_norm",y.norm()},{"parent_to_local_ratio",parent_policy.release_response_norm/y.norm()},
            {"observation_scale",parent_scale},{"normalized_local_norm",normalized_local_response_norm}}},
        {"solver_comparison",table},{"release_policy_comparison",policy_cases},
        {"bounded_refinement_shapes",refinement_shapes},
        {"face_enumeration",j::object{{"faces_examined",faces.faces},{"full_rank_faces",faces.full_rank_faces},
            {"primal_feasible_faces",faces.feasible_faces},{"best_active_mask",faces.best_active_mask},
            {"best_objective_face_kkt_passed",face_summary.certificate.kkt_passed},
            {"objective",Number(face_summary.certificate.objective)},
            {"scaled_kkt",Number(face_summary.certificate.projected_kkt)},
            {"active_a_atoms",IndicesJson(face_summary.certificate.active_atoms)},
            {"best_face_kkt_passed",faces.best_face_kkt_passed},
            {"kkt_feasible_face_count",faces.kkt_feasible_faces},
            {"kkt_face_active_mask",faces.kkt_active_mask},
            {"kkt_face_objective",Number(faces.kkt_objective)},
            {"kkt_face_scaled_kkt",Number(faces.kkt_scaled)},
            {"valid",faces.valid}}},
        {"singular_values",j::object{{"full_design",NumbersJson(spectrum.singular_values)},
            {"free_face",NumbersJson(production_summary.free_spectrum.singular_values)}}},
        {"alternate_free_face_qr",j::object{{"long_double_digits",alternate_face.precision_bits},
            {"rank",alternate_face.rank},{"valid",alternate_face.valid},{"feasible",alternate_face.feasible},
            {"objective",Number(alternate_face_summary.certificate.objective)},
            {"scaled_kkt",Number(alternate_face_summary.certificate.projected_kkt)},
            {"active_a_atoms",IndicesJson(alternate_face_summary.certificate.active_atoms)}}},
        {"coefficient_difference_from_face_reference",j::object{
            {"production",Number(CoefficientDifference(production.beta,faces.beta))},
            {"dense_qr",Number(CoefficientDifference(qr_summary.linear.beta,faces.beta))},
            {"dense_svd",Number(CoefficientDifference(svd_summary.linear.beta,faces.beta))},
            {"dense_free_face",Number(CoefficientDifference(dense_face.beta,faces.beta))},
            {"long_double",Number(CoefficientDifference(alternate_face.beta,faces.beta))}}},
        {"free_face_iterative_refinement",refinement},
        {"assessment_valid",final_assessment_valid},
        {"long_double_digits",std::numeric_limits<long double>::digits}};
}
n::Domain TrialDomainFrom(const j::value & fixture)
{
    const auto rows=j::value_to<Eigen::Index>(fixture.at("identity").at("local_observation_count"));
    std::vector<std::vector<n::Support>> support;
    for(const auto & atom:fixture.at("local_support").as_array())
    {
        std::vector<n::Support> entries;
        for(const auto & item:atom.as_array()) entries.push_back({
            j::value_to<Eigen::Index>(item.at("row")),j::value_to<double>(item.at("squared_distance"))});
        support.push_back(std::move(entries));
    }
    return {rows,std::move(support)};
}
n::EvaluationContext TrialContextFrom(const j::value & fixture)
{
    const auto & value=fixture.at("local_numerical_context"); n::EvaluationContext out;
    out.scale=j::value_to<double>(value.at("local_observation_scale"));
    out.linear.rank_relative=j::value_to<double>(value.at("rank_relative"));
    out.linear.release_factor=j::value_to<double>(value.at("release_factor"));
    out.linear.release_response_norm=j::value_to<double>(value.at("release_response_norm"));
    out.linear.active_set_iteration_factor=j::value_to<int>(value.at("active_set_iteration_factor"));
    out.atom_ids=n::Identities(std::make_shared<const std::vector<std::string>>(StringsFrom(fixture.at("atom_ids"))));
    out.row_ids=n::Identities(std::make_shared<const std::vector<std::string>>(StringsFrom(fixture.at("row_ids"))));
    return out;
}
double RecordedDifference(n::VectorRef value,const j::value & recorded)
{
    const auto expected=VectorFrom(recorded);
    return value.size()==expected.size() ? (value-expected).lpNorm<Eigen::Infinity>() :
        std::numeric_limits<double>::infinity();
}
bool ReplayCutoffRegression()
{
    const double width=0.7,cutoff=std::nextafter(2.5*2.5,std::numeric_limits<double>::infinity());
    n::Domain domain{3,{{{0,0.0},{1,1.0},{2,cutoff}}}};
    n::Vector eta(1),beta(2),y=n::Vector::Zero(3); eta(0)=std::log(width); beta<<1.0,0.0;
    const double matched_width=std::exp(eta(0));
    y(0)=n::EvaluateKernel(0.0,matched_width,2.5).gaussian;
    y(1)=n::EvaluateKernel(1.0,matched_width,2.5).gaussian;
    n::EvaluationContext context; context.scale=std::max(1.0,y.norm());
    const auto primary=n::EvaluateState(domain,y,eta,beta,context);
    n::ReplayDetailsForDiagnostics replay;
    const auto trust=n::CheckReplay(domain,y,primary,context,&replay);
    return primary.valid && primary.certificate.kkt_passed && trust.passed &&
        replay.prediction(2)==0.0 && trust.prediction_difference==0.0;
}
j::value ReplayTrialFixture(const j::value & fixture)
{
    if(j::value_to<std::string>(fixture.at("format"))!="rhbm-joint-local-replay-trial-v1")
        throw std::invalid_argument("Unsupported local replay trial fixture format.");
    const auto domain=TrialDomainFrom(fixture); const auto y=VectorFrom(fixture.at("effective_response"));
    const auto eta=VectorFrom(fixture.at("trial_eta")),beta=VectorFrom(fixture.at("trial_beta"));
    const auto context=TrialContextFrom(fixture);
    const auto primary=n::EvaluateState(domain,y,eta,beta,context);
    n::ReplayDetailsForDiagnostics replay;
    const auto trust=n::CheckReplay(domain,y,primary,context,&replay);
    auto wrong_gradient=primary; wrong_gradient.gradient(0)+=1e-3;
    const auto wrong_gradient_trust=n::CheckReplay(domain,y,wrong_gradient,context);
    auto wrong_kkt=primary; wrong_kkt.certificate.projected_kkt+=1e-6;
    const auto wrong_kkt_trust=n::CheckReplay(domain,y,wrong_kkt,context);
    n::Vector wrong_beta=beta; wrong_beta(0)*=1.1;
    const auto wrong_coefficients=n::EvaluateState(domain,y,eta,wrong_beta,context);
    const auto wrong_coefficients_trust=n::CheckReplay(domain,y,wrong_coefficients,context);
    const bool cutoff_regression=ReplayCutoffRegression();
    const auto primary_prediction=(primary.x*primary.beta).eval();
    n::Vector stable_prediction=n::Vector::Zero(primary.x.rows()),stable_compensation=stable_prediction;
    for(Eigen::Index column=0;column<primary.x.cols();++column)
        for(n::Sparse::InnerIterator entry(primary.x,column);entry;++entry)
        {
            const auto row=entry.row(); const double term=entry.value()*primary.beta(column);
            const double adjusted=term-stable_compensation(row),next=stable_prediction(row)+adjusted;
            stable_compensation(row)=(next-stable_prediction(row))-adjusted; stable_prediction(row)=next;
        }
    const auto reproduced_prediction=RecordedDifference(primary_prediction,fixture.at("primary_prediction"));
    const auto reproduced_residual=RecordedDifference(primary.residual,fixture.at("primary_residual"));
    const auto reproduced_gradient=RecordedDifference(primary.gradient,fixture.at("primary_width_gradient"));
    const auto reproduced_replay_prediction=RecordedDifference(replay.prediction,fixture.at("replay_prediction"));
    const auto reproduced_replay_residual=RecordedDifference(replay.residual,fixture.at("replay_residual"));
    const auto reproduced_replay_gradient=RecordedDifference(replay.width_gradient,fixture.at("replay_width_gradient"));
    return j::object{{"fixture_loaded",true},{"identity",fixture.at("identity")},
        {"local_support_entries",static_cast<std::uint64_t>(std::accumulate(domain.atoms.begin(),domain.atoms.end(),
            std::size_t{},[](std::size_t total,const auto & atom){return total+atom.size();}))},
        {"primary_valid",primary.valid},{"primary_kkt_passed",primary.certificate.kkt_passed},
        {"primary_scaled_kkt",Number(primary.certificate.projected_kkt)},
        {"replay_reason",trust.reason},{"replay_passed",trust.passed},
        {"cutoff_regression_passed",cutoff_regression},
        {"wrong_gradient_rejected",!wrong_gradient_trust.passed},
        {"wrong_kkt_rejected",!wrong_kkt_trust.passed},
        {"wrong_coefficients_rejected",!wrong_coefficients_trust.passed},
        {"prediction_difference",Number(trust.prediction_difference)},
        {"stable_sparse_prediction_replay_difference",(stable_prediction-replay.prediction).lpNorm<Eigen::Infinity>()},
        {"gradient_difference",Number(trust.gradient_difference)},
        {"kkt_difference",Number(trust.kkt_difference)},
        {"recomputed_primary_max_differences",j::object{{"prediction",Number(reproduced_prediction)},
            {"residual",Number(reproduced_residual)},{"width_gradient",Number(reproduced_gradient)}}},
        {"recomputed_replay_max_differences",j::object{{"prediction",Number(reproduced_replay_prediction)},
            {"residual",Number(reproduced_replay_residual)},{"width_gradient",Number(reproduced_replay_gradient)}}},
        {"recorded_trust",fixture.at("trust")}};
}
}

int main(int argc,char ** argv)
{
    try
    {
        if(argc!=3) throw std::invalid_argument("Usage: joint_local_ac_kkt replay|assess|replay-trial FIXTURE.json");
        const std::string action(argv[1]);
        if(action=="replay")
        {
            const auto result=Run(Read(argv[2])); std::cout<<j::serialize(result)<<'\n';
            return result.at("recovered").as_bool() ? 0 : 1;
        }
        if(action=="expect-kkt-rejection")
        {
            const auto result=Run(Read(argv[2]));
            const bool rejected=!result.at("recovered").as_bool() &&
                j::value_to<std::string>(result.at("initial_profile_reason"))=="kkt-failed" &&
                j::value_to<double>(result.at("scaled_projected_kkt"))>1e-10 &&
                result.at("primal_feasible").as_bool();
            std::cout<<j::serialize(j::object{{"expected_kkt_rejection",rejected},
                {"initial_profile_reason",result.at("initial_profile_reason")},
                {"scaled_projected_kkt",result.at("scaled_projected_kkt")},
                {"primal_feasible",result.at("primal_feasible")}})<<'\n';
            return rejected ? 0 : 1;
        }
        if(action=="reference-kkt")
        {
            std::cout<<j::serialize(HighPrecisionReference(Read(argv[2])))<<'\n';
            return 0;
        }
        if(action=="assess")
        {
            const auto result=Assess(Read(argv[2])); std::cout<<j::serialize(result)<<'\n';
            return result.at("assessment_valid").as_bool() ? 0 : 1;
        }
        if(action=="replay-trial")
        {
            const auto result=ReplayTrialFixture(Read(argv[2])); std::cout<<j::serialize(result)<<'\n';
            return result.at("primary_valid").as_bool() && result.at("replay_passed").as_bool() &&
                result.at("cutoff_regression_passed").as_bool() && result.at("wrong_gradient_rejected").as_bool() &&
                result.at("wrong_kkt_rejected").as_bool() && result.at("wrong_coefficients_rejected").as_bool() ? 0 : 1;
        }
        throw std::invalid_argument("Usage: joint_local_ac_kkt replay|assess|replay-trial FIXTURE.json");
    }
    catch(const std::exception & error)
    {std::cerr<<error.what()<<'\n'; return 2;}
}
