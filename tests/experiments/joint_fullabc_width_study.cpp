#include "core/detail/joint_component/Numerics.hpp"
#include <boost/json.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <Eigen/SVD>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <string>

namespace {
namespace j=boost::json;
namespace n=rhbm_gem::core::joint_component;
using MP=boost::multiprecision::cpp_dec_float_50;
using MPMatrix=std::vector<std::vector<MP>>;

j::value Read(const char * path)
{
    std::ifstream input(path);
    if(!input) throw std::runtime_error(std::string("Cannot read fixture: ")+path);
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
    const auto rows=j::value_to<Eigen::Index>(fixture.at("identity").at("local_design_rows"));
    const auto atoms=static_cast<Eigen::Index>(out.atom_ids.size());
    out.rank={rows,2*atoms,atoms};
    out.independent_search=true; out.update_budget=1;
    return out;
}
MP ExactInput(double value)
{
    if(value==0.0) return 0;
    int exponent{}; const double fraction=std::frexp(std::abs(value),&exponent);
    const auto mantissa=static_cast<std::uint64_t>(std::ldexp(fraction,std::numeric_limits<double>::digits));
    MP exact=boost::multiprecision::ldexp(MP(mantissa),exponent-std::numeric_limits<double>::digits);
    return std::signbit(value) ? -exact : exact;
}
std::string MpText(const MP & value)
{return value.str(42,std::ios_base::scientific);}
j::value Number(double value)
{return std::isfinite(value) ? j::value(value) : j::value(nullptr);}
j::array Numbers(const n::Vector & values)
{j::array out; for(Eigen::Index i=0;i<values.size();++i) out.push_back(Number(values(i))); return out;}
j::array Indices(const std::vector<Eigen::Index> & values)
{j::array out; for(const auto value:values) out.push_back(value); return out;}
Eigen::Index FindAtom(const std::vector<std::string> & ids,const std::string & id)
{
    const auto found=std::find(ids.begin(),ids.end(),id);
    if(found==ids.end()) throw std::runtime_error("Required atom is absent from local fixture: "+id);
    return static_cast<Eigen::Index>(found-ids.begin());
}
j::object PairMetrics(const n::Sparse & x,Eigen::Index atom)
{
    n::Vector g=x.col(2*atom),k=x.col(2*atom+1);
    const double gaussian_solver_norm=g.norm(),charge_solver_norm=k.norm();
    const double gn=g.stableNorm(),kn=k.stableNorm();
    double rho=std::numeric_limits<double>::quiet_NaN(),sigma_min=std::numeric_limits<double>::quiet_NaN(),
        condition=std::numeric_limits<double>::quiet_NaN();
    if(gn>0 && kn>0)
    {
        g/=gn; k/=kn; rho=g.dot(k);
        const double sign=rho<0 ? -1.0 : 1.0;
        sigma_min=(g-sign*k).norm()/std::sqrt(2.0);
        const double sigma_max=(g+sign*k).norm()/std::sqrt(2.0);
        condition=sigma_min>0 ? sigma_max/sigma_min : std::numeric_limits<double>::infinity();
    }
    return {{"gaussian_norm",Number(gn)},{"charge_norm",Number(kn)},
        {"gaussian_solver_norm",Number(gaussian_solver_norm)},{"charge_solver_norm",Number(charge_solver_norm)},
        {"gaussian_max_abs",Number(g.size() ? g.cwiseAbs().maxCoeff()*gn : 0.0)},
        {"gaussian_nonzero_entries",x.col(2*atom).nonZeros()},
        {"normalized_correlation",Number(rho)},
        {"normalized_pair_sigma_min",Number(sigma_min)},{"normalized_pair_condition",Number(condition)}};
}
std::vector<Eigen::Index> FreeColumns(const n::Evaluation & evaluation)
{
    std::vector<Eigen::Index> columns;
    for(Eigen::Index col=0;col<evaluation.beta.size();++col)
        if(col%2!=0 || evaluation.beta(col)>0) columns.push_back(col);
    return columns;
}
n::Matrix SelectColumns(const n::Sparse & x,const std::vector<Eigen::Index> & columns,bool normalize)
{
    n::Matrix selected(x.rows(),static_cast<Eigen::Index>(columns.size()));
    for(std::size_t k=0;k<columns.size();++k)
    {
        selected.col(static_cast<Eigen::Index>(k))=x.col(columns[k]);
        if(normalize)
        {
            const double norm=selected.col(static_cast<Eigen::Index>(k)).norm();
            if(norm>0) selected.col(static_cast<Eigen::Index>(k))/=norm;
        }
    }
    return selected;
}
j::object SpectrumJson(const n::Spectrum & spectrum)
{
    return {{"available",spectrum.available},{"rank",spectrum.available ? j::value(spectrum.rank) : j::value(nullptr)},
        {"condition",Number(spectrum.condition)},{"sigma_min",Number(spectrum.minimum)},
        {"singular_values",Numbers(spectrum.singular_values)}};
}
Eigen::Index RelativeRank(const n::Vector & singular_values,const n::RankPolicy & policy,Eigen::Index columns)
{
    if(!singular_values.size()) return 0;
    const double threshold=policy.Relative(columns)*singular_values(0);
    Eigen::Index rank{};
    while(rank<singular_values.size() && singular_values(rank)>threshold) ++rank;
    return rank;
}
struct KernelPair {MP g{},k{};};
KernelPair FunctionalKernel(double square,double width)
{
    if(square>2.5*2.5) return {};
    const MP s=ExactInput(square),b=ExactInput(width),pi=boost::math::constants::pi<MP>();
    const MP b2=b*b,exponential=exp(-s/(MP(2)*b2));
    const MP normalizer=MP(2)*pi*b2;
    KernelPair result; result.g=exponential/(normalizer*sqrt(normalizer));
    if(s<MP("1e-10")) result.k=sqrt(MP(2)/pi)/b;
    else result.k=boost::math::erf(sqrt(s)/(sqrt(MP(2))*b))/sqrt(s);
    return result;
}
std::vector<double> ProfileLogGrid(double minimum,double maximum,double step,double base,double extra)
{
    std::vector<double> values;
    for(double point=minimum;point<=maximum+step*1e-10;point+=step) values.push_back(point);
    values.push_back(base); values.push_back(base-step); values.push_back(base-step/2);
    values.push_back(base+step/2); values.push_back(base+step); values.push_back(extra);
    std::sort(values.begin(),values.end());
    values.erase(std::unique(values.begin(),values.end(),[](double a,double b) {return std::abs(a-b)<1e-12;}),values.end());
    return values;
}
struct HpResult
{
    std::vector<MP> beta;
    std::vector<Eigen::Index> active_atoms;
    MP kkt{},objective{};
    bool converged{};
    std::size_t face_solves{};
};
MP NormColumn(const MPMatrix & x,std::size_t col)
{
    MP sum{}; for(const auto & row:x) sum+=row[col]*row[col]; return sqrt(sum);
}
std::vector<MP> HouseholderSolve(const MPMatrix & full,const std::vector<MP> & y,
    const std::vector<std::size_t> & columns)
{
    const std::size_t rows=full.size(),count=columns.size();
    if(rows<count || !rows) throw std::runtime_error("Invalid bounded multiprecision QR shape.");
    MPMatrix a(rows,std::vector<MP>(count)); std::vector<MP> norms(count),rhs=y;
    for(std::size_t j=0;j<count;++j)
    {
        norms[j]=NormColumn(full,columns[j]);
        if(norms[j]==0) throw std::runtime_error("Multiprecision free face has a zero column.");
        for(std::size_t i=0;i<rows;++i) a[i][j]=full[i][columns[j]]/norms[j];
    }
    for(std::size_t j=0;j<count;++j)
    {
        MP norm2{}; for(std::size_t i=j;i<rows;++i) norm2+=a[i][j]*a[i][j];
        const MP norm=sqrt(norm2); if(norm==0) throw std::runtime_error("Multiprecision QR found a zero pivot.");
        std::vector<MP> v(rows-j);
        for(std::size_t i=j;i<rows;++i) v[i-j]=a[i][j];
        v[0]+=a[j][j]>=0 ? norm : -norm;
        MP vnorm2{}; for(const auto & value:v) vnorm2+=value*value;
        if(vnorm2==0) throw std::runtime_error("Multiprecision QR reflector vanished.");
        const MP tau=MP(2)/vnorm2;
        for(std::size_t col=j;col<count;++col)
        {
            MP dot{}; for(std::size_t i=j;i<rows;++i) dot+=v[i-j]*a[i][col];
            for(std::size_t i=j;i<rows;++i) a[i][col]-=tau*v[i-j]*dot;
        }
        MP dot{}; for(std::size_t i=j;i<rows;++i) dot+=v[i-j]*rhs[i];
        for(std::size_t i=j;i<rows;++i) rhs[i]-=tau*v[i-j]*dot;
    }
    std::vector<MP> theta(count),beta(count);
    for(std::size_t row=count;row-->0;)
    {
        MP value=rhs[row];
        for(std::size_t col=row+1;col<count;++col) value-=a[row][col]*theta[col];
        if(a[row][row]==0) throw std::runtime_error("Multiprecision QR is rank deficient.");
        theta[row]=value/a[row][row]; beta[row]=theta[row]/norms[row];
    }
    return beta;
}
MP Objective(const MPMatrix & x,const std::vector<MP> & y,const std::vector<MP> & beta)
{
    MP objective{};
    for(std::size_t row=0;row<x.size();++row)
    {
        MP residual=-y[row];
        for(std::size_t col=0;col<beta.size();++col) residual+=x[row][col]*beta[col];
        objective+=residual*residual/2;
    }
    return objective;
}
MP Kkt(const MPMatrix & x,const std::vector<MP> & y,const std::vector<MP> & beta,double scale,
    const std::vector<Eigen::Index> & active)
{
    std::vector<bool> is_active(beta.size()/2);
    for(const auto atom:active) is_active[static_cast<std::size_t>(atom)]=true;
    std::vector<MP> residual(y.size());
    for(std::size_t row=0;row<x.size();++row)
    {
        residual[row]=-y[row];
        for(std::size_t col=0;col<beta.size();++col) residual[row]+=x[row][col]*beta[col];
    }
    MP maximum{};
    for(std::size_t col=0;col<beta.size();++col)
    {
        const MP norm=NormColumn(x,col); if(norm==0) continue;
        MP gradient{}; for(std::size_t row=0;row<x.size();++row) gradient+=x[row][col]*residual[row];
        const MP scaled=gradient/norm/ExactInput(scale);
        const MP candidate=col%2==0 && is_active[col/2] ? MP(-scaled) : MP(abs(scaled));
        const MP error=col%2==0 && is_active[col/2] && candidate<0 ? MP(0) : candidate;
        if(error>maximum) maximum=error;
    }
    return maximum;
}
HpResult ConstrainedQr(const MPMatrix & x,const std::vector<MP> & y,std::size_t atom_count,double scale,
    const std::vector<Eigen::Index> & initial_active)
{
    const std::size_t columns=x.front().size();
    std::vector<bool> active(atom_count);
    for(const auto atom:initial_active) if(atom>=0 && static_cast<std::size_t>(atom)<atom_count)
        active[static_cast<std::size_t>(atom)]=true;
    for(std::size_t atom=0;atom<atom_count;++atom) if(NormColumn(x,2*atom)==0) active[atom]=true;
    HpResult result;
    for(std::size_t iteration=0;iteration<4*atom_count+8;++iteration)
    {
        std::vector<std::size_t> free;
        for(std::size_t col=0;col<columns;++col)
            if(col%2!=0 || !active[col/2]) free.push_back(col);
        const auto free_beta=HouseholderSolve(x,y,free); ++result.face_solves;
        result.beta.assign(columns,MP(0));
        for(std::size_t k=0;k<free.size();++k) result.beta[free[k]]=free_beta[k];
        std::size_t most_negative=atom_count; MP worst{};
        for(std::size_t atom=0;atom<atom_count;++atom)
            if(!active[atom] && result.beta[2*atom]<0 &&
               (most_negative==atom_count || result.beta[2*atom]<worst))
            {most_negative=atom; worst=result.beta[2*atom];}
        if(most_negative!=atom_count) {active[most_negative]=true; continue;}

        std::vector<MP> residual(y.size());
        for(std::size_t row=0;row<x.size();++row)
        {
            residual[row]=-y[row];
            for(std::size_t col=0;col<columns;++col) residual[row]+=x[row][col]*result.beta[col];
        }
        MP most_violated{}; std::size_t release=atom_count; result.kkt=0;
        for(std::size_t col=0;col<columns;++col)
        {
            const MP norm=NormColumn(x,col); if(norm==0) continue;
            MP gradient{}; for(std::size_t row=0;row<x.size();++row) gradient+=x[row][col]*residual[row];
            const MP scaled=gradient/norm/ExactInput(scale);
            if(col%2==0 && active[col/2])
            {
                const MP signed_violation=-scaled;
                const MP violation=signed_violation>0 ? signed_violation : MP(0);
                if(violation>result.kkt) result.kkt=violation;
                if(violation>most_violated && violation>MP("1e-40"))
                {most_violated=violation; release=col/2;}
            }
            else {const MP error=abs(scaled); if(error>result.kkt) result.kkt=error;}
        }
        if(release!=atom_count) {active[release]=false; continue;}
        result.objective=0;
        for(const auto & value:residual) result.objective+=value*value/2;
        result.active_atoms.clear();
        for(std::size_t atom=0;atom<atom_count;++atom) if(active[atom])
            result.active_atoms.push_back(static_cast<Eigen::Index>(atom));
        result.converged=true; return result;
    }
    return result;
}
MPMatrix DoubleDesign(const n::Sparse & design)
{
    MPMatrix x(static_cast<std::size_t>(design.rows()),std::vector<MP>(static_cast<std::size_t>(design.cols())));
    for(Eigen::Index col=0;col<design.cols();++col)
        for(n::Sparse::InnerIterator value(design,col);value;++value)
            x[static_cast<std::size_t>(value.row())][static_cast<std::size_t>(col)]=ExactInput(value.value());
    return x;
}
MPMatrix FunctionalDesign(const n::Domain & domain,const n::Vector & eta)
{
    const std::size_t rows=static_cast<std::size_t>(domain.rows),atoms=domain.atoms.size();
    MPMatrix x(rows,std::vector<MP>(2*atoms));
    for(std::size_t atom=0;atom<atoms;++atom)
    {
        const double width=std::exp(eta(static_cast<Eigen::Index>(atom)));
        for(const auto & support:domain.atoms[atom])
        {
            const auto pair=FunctionalKernel(support.square,width);
            x[static_cast<std::size_t>(support.row)][2*atom]=pair.g;
            x[static_cast<std::size_t>(support.row)][2*atom+1]=pair.k;
        }
    }
    return x;
}
MP WidthGradient(const MPMatrix & x,const std::vector<MP> & y,const std::vector<MP> & beta,
    const std::vector<MP> & derivative,double scale)
{
    MP gradient{}; const MP scale_squared=ExactInput(scale)*ExactInput(scale);
    for(std::size_t row=0;row<x.size();++row)
    {
        MP residual=-y[row];
        for(std::size_t col=0;col<beta.size();++col) residual+=x[row][col]*beta[col];
        gradient+=residual*derivative[row]/scale_squared;
    }
    return gradient;
}
MP FixedDesignWidthGradient(const n::Evaluation & evaluation,const MPMatrix & design,
    const std::vector<MP> & y,const std::vector<MP> & beta,Eigen::Index atom,double scale)
{
    std::vector<MP> derivative(static_cast<std::size_t>(evaluation.x.rows()));
    for(Eigen::Index column=2*atom;column<2*atom+2;++column)
        for(n::Sparse::InnerIterator value(evaluation.derivative,column);value;++value)
            derivative[static_cast<std::size_t>(value.row())]+=
                beta[static_cast<std::size_t>(column)]*ExactInput(value.value());
    return WidthGradient(design,y,beta,derivative,scale);
}
MP FunctionalWidthGradient(const n::Domain & domain,const n::Vector & eta,const MPMatrix & design,
    const std::vector<MP> & y,const std::vector<MP> & beta,Eigen::Index atom,double scale)
{
    const std::size_t rows=static_cast<std::size_t>(domain.rows);
    const double width=std::exp(eta(atom)); const MP b=ExactInput(width),b2=b*b;
    const MP pi=boost::math::constants::pi<MP>(),center=sqrt(MP(2)/pi)/b;
    std::vector<MP> derivative(rows);
    for(const auto & support:domain.atoms[static_cast<std::size_t>(atom)])
    {
        const auto pair=FunctionalKernel(support.square,width);
        const MP square=ExactInput(support.square),exponential=exp(-square/(MP(2)*b2));
        const MP dg=pair.g*(square/b2-MP(3)),dk=-center*exponential;
        derivative[static_cast<std::size_t>(support.row)]+=
            beta[static_cast<std::size_t>(2*atom)]*dg+beta[static_cast<std::size_t>(2*atom+1)]*dk;
    }
    return WidthGradient(design,y,beta,derivative,scale);
}
std::vector<MP> MpVector(const n::Vector & values)
{
    std::vector<MP> out(static_cast<std::size_t>(values.size()));
    for(Eigen::Index i=0;i<values.size();++i) out[static_cast<std::size_t>(i)]=ExactInput(values(i));
    return out;
}
j::object HighPrecisionReference(const n::Domain & domain,const n::Vector & y,const n::Vector & eta,
    const n::Evaluation & production,const n::EvaluationContext & context)
{
    if(domain.rows>2000 || production.x.cols()>64 || domain.rows*production.x.cols()>100000)
        throw std::runtime_error("High-precision reference exceeds its bounded 2,000 x 64 design budget.");
    const auto hp_y=MpVector(y); const auto active=production.certificate.active_atoms;
    const auto rounded_design=DoubleDesign(production.x);
    const auto fixed=ConstrainedQr(rounded_design,hp_y,domain.atoms.size(),context.scale,active);
    const auto functional_design=FunctionalDesign(domain,eta);
    const auto functional=ConstrainedQr(functional_design,hp_y,domain.atoms.size(),context.scale,fixed.active_atoms);
    if(!fixed.converged || !functional.converged)
        throw std::runtime_error("50-digit Householder QR did not satisfy constrained KKT.");
    const auto production_beta=MpVector(production.beta);
    std::vector<MP> rounded_beta(fixed.beta.size());
    for(std::size_t col=0;col<fixed.beta.size();++col)
        rounded_beta[col]=ExactInput(fixed.beta[col].convert_to<double>());
    n::Vector rounded_hp_beta(static_cast<Eigen::Index>(rounded_beta.size()));
    for(std::size_t col=0;col<rounded_beta.size();++col)
        rounded_hp_beta(static_cast<Eigen::Index>(col))=rounded_beta[col].convert_to<double>();
    const auto rounded_hp_certificate=n::CertifyLinear(production.x,y,rounded_hp_beta,context.scale);
    MP max_kernel_abs{},max_kernel_rel{},max_production_prediction_difference{};
    for(std::size_t row=0;row<rounded_design.size();++row)
    {
        MP fixed_prediction{},functional_prediction{};
        for(std::size_t col=0;col<rounded_design[row].size();++col)
        {
            const MP difference=abs(functional_design[row][col]-rounded_design[row][col]);
            MP denominator=abs(functional_design[row][col]); if(denominator<MP("1e-300")) denominator=MP("1e-300");
            const MP relative=difference/denominator;
            if(difference>max_kernel_abs) max_kernel_abs=difference;
            if(relative>max_kernel_rel) max_kernel_rel=relative;
            fixed_prediction+=rounded_design[row][col]*production_beta[col];
            functional_prediction+=functional_design[row][col]*production_beta[col];
        }
        const MP difference=abs(functional_prediction-fixed_prediction);
        if(difference>max_production_prediction_difference) max_production_prediction_difference=difference;
    }
    MP max_fixed_scaled_coefficient_delta{},max_functional_scaled_coefficient_delta{};
    for(std::size_t col=0;col<fixed.beta.size();++col)
    {
        const MP norm=NormColumn(rounded_design,col);
        const MP fixed_delta=norm*abs(production_beta[col]-fixed.beta[col]);
        const MP functional_delta=norm*abs(functional.beta[col]-fixed.beta[col]);
        if(fixed_delta>max_fixed_scaled_coefficient_delta) max_fixed_scaled_coefficient_delta=fixed_delta;
        if(functional_delta>max_functional_scaled_coefficient_delta) max_functional_scaled_coefficient_delta=functional_delta;
    }
    const auto text=[](const MP & value){return MpText(value);};
    j::object fixed_json{{"face_solves",fixed.face_solves},
        {"active_a_atom_indices",Indices(fixed.active_atoms)},
        {"high_precision_kkt",text(fixed.kkt)},
        {"production_beta_high_precision_kkt",text(Kkt(rounded_design,hp_y,production_beta,context.scale,active))},
        {"rounded_hp_solution_double_beta_high_precision_kkt",text(Kkt(rounded_design,hp_y,rounded_beta,context.scale,fixed.active_atoms))},
        {"rounded_hp_solution_double_beta_binary64_kkt",rounded_hp_certificate.projected_kkt},
        {"high_precision_objective",text(fixed.objective)},
        {"production_beta_objective",text(Objective(rounded_design,hp_y,production_beta))},
        {"rounded_hp_solution_objective",text(Objective(rounded_design,hp_y,rounded_beta))},
        {"max_column_norm_scaled_production_beta_delta",text(max_fixed_scaled_coefficient_delta)}};
    j::object functional_json{{"face_solves",functional.face_solves},
        {"active_a_atom_indices",Indices(functional.active_atoms)},
        {"high_precision_kkt",text(functional.kkt)},
        {"high_precision_objective",text(functional.objective)},
        {"production_beta_high_precision_kkt",text(Kkt(functional_design,hp_y,production_beta,context.scale,functional.active_atoms))},
        {"fixed_design_reference_beta_high_precision_kkt",text(Kkt(functional_design,hp_y,fixed.beta,context.scale,functional.active_atoms))},
        {"max_kernel_absolute_difference_from_rounded_double_design",text(max_kernel_abs)},
        {"max_kernel_relative_difference_from_rounded_double_design",text(max_kernel_rel)},
        {"max_production_prediction_difference_from_rounded_double_design",text(max_production_prediction_difference)},
        {"max_column_norm_scaled_functional_beta_delta_from_fixed_design_reference",text(max_functional_scaled_coefficient_delta)}};
    return {{"precision","boost::multiprecision::cpp_dec_float_50"},
        {"method","column-normalized Householder QR; active-set constrained KKT checked"},
        {"fixed_rounded_double_design",std::move(fixed_json)},
        {"functional_kernel_design",std::move(functional_json)}};
}
j::array ProfileRows(const n::Domain & domain,const n::Vector & y,const n::Vector & base_eta,
    const n::EvaluationContext & context,const std::vector<std::string> & atom_ids,Eigen::Index atom,
    double log_min,double log_max,double log_step,double extra_width,
    const std::vector<double> * prescribed_widths=nullptr)
{
    j::array rows;
    const auto hp_y=MpVector(y);
    std::vector<double> points;
    if(prescribed_widths)
    {
        for(const auto width:*prescribed_widths)
            if(width>0 && std::isfinite(width)) points.push_back(std::log(width));
        std::sort(points.begin(),points.end());
        points.erase(std::unique(points.begin(),points.end(),[](double a,double b) {
            return std::abs(a-b)<1e-12;
        }),points.end());
    }
    else points=ProfileLogGrid(log_min,log_max,log_step,base_eta(atom),std::log(extra_width));
    for(const double log_width:points)
    {
        auto eta=base_eta; eta(atom)=log_width;
        const auto evaluation=n::EvaluateProfile(domain,y,eta,false,&context);
        const auto spectrum=n::ComputeSpectrum(evaluation.x,context.rank,evaluation.x.cols(),true);
        const double a=evaluation.beta.size() ? evaluation.beta(2*atom) : n::unavailable;
        const double c=evaluation.beta.size() ? evaluation.beta(2*atom+1) : n::unavailable;
        const auto pair=PairMetrics(evaluation.x,atom);
        const double normalized_a=std::isfinite(a) ? a*j::value_to<double>(pair.at("gaussian_norm")) : n::unavailable;
        const double normalized_c=std::isfinite(c) ? c*j::value_to<double>(pair.at("charge_norm")) : n::unavailable;
        const auto fixed_design=DoubleDesign(evaluation.x);
        const auto fixed_reference=ConstrainedQr(fixed_design,hp_y,domain.atoms.size(),context.scale,
            evaluation.certificate.active_atoms);
        const auto functional_design=FunctionalDesign(domain,eta);
        const auto functional_reference=ConstrainedQr(functional_design,hp_y,domain.atoms.size(),context.scale,
            fixed_reference.active_atoms);
        const auto production_beta=evaluation.beta.allFinite() ? MpVector(evaluation.beta) : std::vector<MP>{};
        const MP objective_scale=ExactInput(context.scale)*ExactInput(context.scale);
        rows.push_back(j::object{{"atom_id",atom_ids[static_cast<std::size_t>(atom)]},
            {"log_width",log_width},{"width",Number(std::exp(log_width))},
            {"profile_loss",Number(evaluation.certificate.available ? evaluation.certificate.objective/(context.scale*context.scale) : n::unavailable)},
            {"inner_kkt",Number(evaluation.certificate.projected_kkt)},
            {"profile_gradient_magnitude",Number(evaluation.gradient.size()>atom ? std::abs(evaluation.gradient(atom)) : n::unavailable)},
            {"profile_valid",evaluation.valid},{"profile_reason",evaluation.reason},
            {"active_a_atom_indices",Indices(evaluation.certificate.active_atoms)},
            {"local_full_rank",spectrum.available ? j::value(spectrum.rank) : j::value(nullptr)},
            {"local_full_condition",Number(spectrum.condition)},
            {"A",Number(a)},{"C",Number(c)},{"normalized_A",Number(normalized_a)},
            {"normalized_C",Number(normalized_c)},{"pair",pair},
            {"fixed_double_design_hp_loss",fixed_reference.converged ? MpText(fixed_reference.objective/objective_scale) : "unavailable"},
            {"fixed_double_design_hp_kkt",fixed_reference.converged ? MpText(fixed_reference.kkt) : "unavailable"},
            {"fixed_double_design_hp_active_a",Indices(fixed_reference.active_atoms)},
            {"fixed_design_production_beta_hp_kkt",fixed_reference.converged && !production_beta.empty() ?
                MpText(Kkt(fixed_design,hp_y,production_beta,context.scale,fixed_reference.active_atoms)) : "unavailable"},
            {"fixed_double_design_hp_width_gradient",fixed_reference.converged ?
                MpText(FixedDesignWidthGradient(evaluation,fixed_design,hp_y,fixed_reference.beta,atom,context.scale)) : "unavailable"},
            {"functional_kernel_hp_loss",functional_reference.converged ? MpText(functional_reference.objective/objective_scale) : "unavailable"},
            {"functional_kernel_hp_kkt",functional_reference.converged ? MpText(functional_reference.kkt) : "unavailable"},
            {"functional_kernel_hp_width_gradient",functional_reference.converged ?
                MpText(FunctionalWidthGradient(domain,eta,functional_design,hp_y,functional_reference.beta,atom,context.scale)) : "unavailable"},
            {"functional_kernel_hp_active_a",Indices(functional_reference.active_atoms)}});
    }
    return rows;
}
j::object ProfileCurvature(const j::array & rows,double base_log_width,double step)
{
    const double h=step/2;
    const auto at=[&](double log_width)->const j::object & {
        const auto found=std::find_if(rows.begin(),rows.end(),[&](const auto & value) {
            return std::abs(j::value_to<double>(value.at("log_width"))-log_width)<1e-12;
        });
        if(found==rows.end()) throw std::runtime_error("Profile curvature grid point is absent.");
        return found->as_object();
    };
    const auto & minus=at(base_log_width-h); const auto & center=at(base_log_width);
    const auto & plus=at(base_log_width+h);
    const auto high_precision=[](const j::object & value,const char * field) {
        return MP(value.at(field).as_string().c_str());
    };
    const MP denominator=MP(h)*MP(h);
    const auto fixed_minus=high_precision(minus,"fixed_double_design_hp_loss");
    const auto fixed_center=high_precision(center,"fixed_double_design_hp_loss");
    const auto fixed_plus=high_precision(plus,"fixed_double_design_hp_loss");
    const auto functional_minus=high_precision(minus,"functional_kernel_hp_loss");
    const auto functional_center=high_precision(center,"functional_kernel_hp_loss");
    const auto functional_plus=high_precision(plus,"functional_kernel_hp_loss");
    return {{"half_step",h},
        {"fixed_double_design_centered_slope",MpText((fixed_plus-fixed_minus)/(MP(2)*MP(h)))},
        {"fixed_double_design_centered_curvature",MpText((fixed_plus-MP(2)*fixed_center+fixed_minus)/denominator)},
        {"functional_kernel_centered_slope",MpText((functional_plus-functional_minus)/(MP(2)*MP(h)))},
        {"functional_kernel_centered_curvature",MpText((functional_plus-MP(2)*functional_center+functional_minus)/denominator)},
        {"functional_kernel_analytic_center_gradient",center.at("functional_kernel_hp_width_gradient")},
        {"fixed_double_design_analytic_center_gradient",center.at("fixed_double_design_hp_width_gradient")}};
}
j::object WidthInformation(const n::Evaluation & state,const n::EvaluationContext & context,
    const std::vector<std::string> & atom_ids)
{
    const auto free=FreeColumns(state); const auto x=SelectColumns(state.x,free,true);
    Eigen::JacobiSVD<n::Matrix> ac_svd(x,Eigen::ComputeThinU|Eigen::ComputeThinV);
    const Eigen::Index rank=RelativeRank(ac_svd.singularValues(),context.rank,x.cols());
    n::Matrix d(state.x.rows(),state.eta.size()); d.setZero();
    for(Eigen::Index atom=0;atom<state.eta.size();++atom)
    {
        for(n::Sparse::InnerIterator entry(state.derivative,2*atom);entry;++entry)
            d(entry.row(),atom)+=state.beta(2*atom)*entry.value();
        for(n::Sparse::InnerIterator entry(state.derivative,2*atom+1);entry;++entry)
            d(entry.row(),atom)+=state.beta(2*atom+1)*entry.value();
    }
    n::Matrix projected=d;
    if(rank>0) projected.noalias()-=ac_svd.matrixU().leftCols(rank)*
        (ac_svd.matrixU().leftCols(rank).transpose()*d);
    n::Matrix normalized_d=d,normalized_projected=projected;
    n::Vector d_norm(state.eta.size()),projected_norm(state.eta.size()),retained_fraction(state.eta.size());
    for(Eigen::Index atom=0;atom<state.eta.size();++atom)
    {
        d_norm(atom)=d.col(atom).norm(); projected_norm(atom)=projected.col(atom).norm();
        retained_fraction(atom)=d_norm(atom)>0 ? projected_norm(atom)/d_norm(atom) : 0;
        if(d_norm(atom)>0) normalized_d.col(atom)/=d_norm(atom);
        if(projected_norm(atom)>0) normalized_projected.col(atom)/=projected_norm(atom);
    }
    Eigen::JacobiSVD<n::Matrix> projected_svd(normalized_projected,Eigen::ComputeThinU|Eigen::ComputeThinV);
    Eigen::JacobiSVD<n::Matrix> raw_projected_svd(projected,Eigen::ComputeThinU|Eigen::ComputeThinV);
    n::Matrix combined_free(x.rows(),x.cols()+normalized_d.cols()); combined_free<<x,normalized_d;
    Eigen::JacobiSVD<n::Matrix> combined_free_svd(combined_free,Eigen::ComputeThinU|Eigen::ComputeThinV);
    const auto full_x=SelectColumns(state.x,[&] {
        std::vector<Eigen::Index> all(static_cast<std::size_t>(state.x.cols()));
        std::iota(all.begin(),all.end(),Eigen::Index{0}); return all;
    }(),true);
    n::Matrix combined_full(full_x.rows(),full_x.cols()+normalized_d.cols()); combined_full<<full_x,normalized_d;
    Eigen::JacobiSVD<n::Matrix> combined_full_svd(combined_full,Eigen::ComputeThinU|Eigen::ComputeThinV);
    const auto rank_and_condition=[&](const n::Vector & singular_values,Eigen::Index columns) {
        const Eigen::Index numerical_rank=RelativeRank(singular_values,context.rank,columns);
        const double condition=singular_values.size() && singular_values.tail(1)(0)>0 ?
            singular_values(0)/singular_values.tail(1)(0) : n::unavailable;
        return std::pair<Eigen::Index,double>{numerical_rank,condition};
    };
    const auto projected_spectrum=rank_and_condition(raw_projected_svd.singularValues(),projected.cols());
    const auto normalized_projected_spectrum=rank_and_condition(projected_svd.singularValues(),normalized_projected.cols());
    const auto free_jacobian_spectrum=rank_and_condition(combined_free_svd.singularValues(),combined_free.cols());
    const auto full_jacobian_spectrum=rank_and_condition(combined_full_svd.singularValues(),combined_full.cols());
    j::array per_atom;
    for(Eigen::Index atom=0;atom<state.eta.size();++atom)
        per_atom.push_back(j::object{{"atom_id",atom_ids[static_cast<std::size_t>(atom)]},
            {"width_derivative_norm",Number(d_norm(atom))},{"projected_width_derivative_norm",Number(projected_norm(atom))},
            {"projected_fraction",Number(retained_fraction(atom))},
            {"profile_gradient_magnitude",Number(state.gradient.size()>atom ? std::abs(state.gradient(atom)) : n::unavailable)}});
    return {{"free_amplitude_columns",free.size()},{"normalized_free_amplitude_rank",rank},
        {"normalized_free_amplitude_sigma_min",ac_svd.singularValues().size() ? ac_svd.singularValues().tail(1)(0) : 0.0},
        {"width_derivative_column_norms",Numbers(d_norm)},
        {"projected_width_derivative_column_norms",Numbers(projected_norm)},
        {"projected_width_raw_singular_values",Numbers(raw_projected_svd.singularValues())},
        {"projected_width_column_normalized_singular_values",Numbers(projected_svd.singularValues())},
        {"combined_free_face_normalized_jacobian_singular_values",Numbers(combined_free_svd.singularValues())},
        {"combined_full_face_normalized_jacobian_singular_values",Numbers(combined_full_svd.singularValues())},
        {"projected_width_raw_numerical_rank",projected_spectrum.first},
        {"projected_width_raw_condition",Number(projected_spectrum.second)},
        {"projected_width_column_normalized_numerical_rank",normalized_projected_spectrum.first},
        {"projected_width_column_normalized_condition",Number(normalized_projected_spectrum.second)},
        {"combined_free_face_normalized_jacobian_numerical_rank",free_jacobian_spectrum.first},
        {"combined_free_face_normalized_jacobian_condition",Number(free_jacobian_spectrum.second)},
        {"combined_full_face_normalized_jacobian_numerical_rank",full_jacobian_spectrum.first},
        {"combined_full_face_normalized_jacobian_condition",Number(full_jacobian_spectrum.second)},
        {"per_atom",std::move(per_atom)}};
}
void Run(const char * fixture_path,const std::filesystem::path & output)
{
    if(std::filesystem::exists(output)) throw std::runtime_error("Use a new output directory for this bounded study.");
    std::filesystem::create_directories(output);
    const auto fixture=Read(fixture_path);
    if(j::value_to<std::string>(fixture.at("format"))!="rhbm-joint-local-ac-kkt-v1")
        throw std::runtime_error("Unsupported Block 21 A/C fixture.");
    const auto domain=DomainFrom(fixture); const auto y=VectorFrom(fixture.at("effective_response"));
    const auto eta=VectorFrom(fixture.at("old_eta")); const auto context=ContextFrom(fixture);
    const auto atom_ids=StringsFrom(fixture.at("atom_ids"));
    const auto profile=n::EvaluateProfile(domain,y,eta,false,&context);
    const auto production_beta=VectorFrom(fixture.at("production_solver").at("beta"));
    const auto state=n::EvaluateState(domain,y,eta,production_beta,context);
    const auto full=n::ComputeSpectrum(profile.x,context.rank,profile.x.cols(),true);
    const auto free_columns=FreeColumns(profile);
    const auto free_design=SelectColumns(profile.x,free_columns,false);
    const auto free_spectrum=n::ComputeSpectrum(free_design,context.rank,free_design.cols(),true);
    const auto atom_11721=FindAtom(atom_ids,"11721"),atom_11650=FindAtom(atom_ids,"11650");
    j::array atom_pairs;
    for(Eigen::Index atom=0;atom<eta.size();++atom)
    {
        const double a=profile.beta.size() ? profile.beta(2*atom) : n::unavailable;
        const double c=profile.beta.size() ? profile.beta(2*atom+1) : n::unavailable;
        auto pair=PairMetrics(profile.x,atom);
        pair["atom_id"]=atom_ids[static_cast<std::size_t>(atom)];
        pair["width"]=std::exp(eta(atom)); pair["A"]=Number(a); pair["C"]=Number(c);
        pair["normalized_A"]=Number(a*j::value_to<double>(pair.at("gaussian_norm")));
        pair["normalized_C"]=Number(c*j::value_to<double>(pair.at("charge_norm")));
        atom_pairs.push_back(std::move(pair));
    }
    const auto profile_11721=ProfileRows(domain,y,eta,context,atom_ids,atom_11721,
        std::log(0.1),std::log(1.0e6),0.5,0.49423856328453153);
    const auto profile_11650=ProfileRows(domain,y,eta,context,atom_ids,atom_11650,
        std::log(0.01),std::log(0.8),0.25,0.4101723088203541);
    const auto high_precision=HighPrecisionReference(domain,y,eta,profile,context);
    const auto min_square=[&](Eigen::Index atom) {
        double minimum=std::numeric_limits<double>::infinity();
        for(const auto & value:domain.atoms[static_cast<std::size_t>(atom)]) minimum=std::min(minimum,value.square);
        return minimum;
    };
    j::array profiles; for(const auto & row:profile_11721) profiles.push_back(row);
    for(const auto & row:profile_11650) profiles.push_back(row);
    std::ofstream csv(output/"profile-slices.csv");
    csv.exceptions(std::ios::failbit|std::ios::badbit);
    csv<<"atom_id,log_width,width,profile_loss,inner_kkt,profile_valid,profile_reason,local_full_rank,local_full_condition,A,C,normalized_A,normalized_C,correlation,pair_sigma_min,fixed_double_design_hp_loss,fixed_double_design_hp_kkt,fixed_design_production_beta_hp_kkt,fixed_double_design_hp_width_gradient,functional_kernel_hp_loss,functional_kernel_hp_kkt,functional_kernel_hp_width_gradient\n";
    for(const auto & row:profiles)
    {
        const auto & value=row.as_object(); const auto & pair=value.at("pair").as_object();
        const auto print=[&](const char * key) {return value.at(key).is_null() ? std::string("nan") : j::serialize(value.at(key));};
        csv<<j::serialize(value.at("atom_id"))<<','<<print("log_width")<<','<<print("width")<<','
            <<print("profile_loss")<<','<<print("inner_kkt")<<','<<j::serialize(value.at("profile_valid"))<<','
            <<j::serialize(value.at("profile_reason"))<<','<<print("local_full_rank")<<','<<print("local_full_condition")<<','
            <<print("A")<<','<<print("C")<<','<<print("normalized_A")<<','<<print("normalized_C")<<','
            <<(pair.at("normalized_correlation").is_null() ? "nan" : j::serialize(pair.at("normalized_correlation")))<<','
            <<(pair.at("normalized_pair_sigma_min").is_null() ? "nan" : j::serialize(pair.at("normalized_pair_sigma_min")))<<','
            <<j::serialize(value.at("fixed_double_design_hp_loss"))<<','
            <<j::serialize(value.at("fixed_double_design_hp_kkt"))<<','
            <<j::serialize(value.at("fixed_design_production_beta_hp_kkt"))<<','
            <<j::serialize(value.at("fixed_double_design_hp_width_gradient"))<<','
            <<j::serialize(value.at("functional_kernel_hp_loss"))<<','
            <<j::serialize(value.at("functional_kernel_hp_kkt"))<<','
            <<j::serialize(value.at("functional_kernel_hp_width_gradient"))<<'\n';
    }
    csv.close();
    j::object summary{{"production_active_a_atom_indices",fixture.at("production_solver").at("active_a_atom_indices")},
        {"production_scaled_kkt",fixture.at("production_solver").at("scaled_projected_kkt")},
        {"profile_scaled_kkt",Number(profile.certificate.projected_kkt)},
        {"profile_reason",profile.reason},{"profile_valid",profile.valid},
        {"profile_loss",Number(profile.certificate.objective/(context.scale*context.scale))},
        {"local_rows",domain.rows},{"local_columns",profile.x.cols()},{"parent_observation_scale",context.scale},
        {"response_norm",y.norm()},{"full_normalized_design",SpectrumJson(full)},
        {"free_normalized_design",SpectrumJson(free_spectrum)},{"atom_pair_metrics",std::move(atom_pairs)},
        {"atom_11650_minimum_squared_distance",min_square(atom_11650)},
        {"atom_11650_minimum_radius",std::sqrt(min_square(atom_11650))},
        {"atom_11650_gaussian_column_norm",profile.x.col(2*atom_11650).norm()},
        {"atom_11650_gaussian_norm_over_observation_scale",profile.x.col(2*atom_11650).norm()/context.scale},
        {"fixed_ac_profile_width_11721",profile_11721},{"fixed_ac_profile_width_11650",profile_11650},
        {"atom_11721_profile_curvature_at_captured_width",ProfileCurvature(profile_11721,eta(atom_11721),0.5)},
        {"atom_11650_profile_curvature_at_captured_width",ProfileCurvature(profile_11650,eta(atom_11650),0.25)},
        {"width_information_at_captured_production_beta",WidthInformation(state,context,atom_ids)},
        {"high_precision_constrained_reference",high_precision},
        {"scope","Conditional slices: vary one atom log-width; hold other widths fixed; re-solve constrained A/C."}};
    std::ofstream out(output/"summary.json"); out.exceptions(std::ios::failbit|std::ios::badbit);
    out<<j::serialize(summary)<<'\n'; out.close();
}
void RunTarget(const char * fixture_path,const std::filesystem::path & output,const std::string & atom_id)
{
    if(std::filesystem::exists(output)) throw std::runtime_error("Use a new output directory for this bounded study.");
    std::filesystem::create_directories(output);
    const auto fixture=Read(fixture_path);
    if(j::value_to<std::string>(fixture.at("format"))!="rhbm-joint-local-ac-kkt-v1")
        throw std::runtime_error("Unsupported local A/C fixture.");
    const auto domain=DomainFrom(fixture); const auto y=VectorFrom(fixture.at("effective_response"));
    const auto eta=VectorFrom(fixture.at("old_eta")); const auto context=ContextFrom(fixture);
    const auto atom_ids=StringsFrom(fixture.at("atom_ids")); const auto atom=FindAtom(atom_ids,atom_id);
    auto widths=std::vector<double>{.1,.2,.3,.5,1,2,5,10,100,1000,std::exp(eta(atom)),
        std::exp(eta(atom)-.25),std::exp(eta(atom)+.25)};
    const auto profile=ProfileRows(domain,y,eta,context,atom_ids,atom,.0,.0,.0,.5,&widths);
    std::ofstream csv(output/"target-width-profile.csv");
    csv.exceptions(std::ios::failbit|std::ios::badbit);
    csv<<"atom_id,width,profile_loss,inner_kkt,profile_valid,profile_reason,local_full_rank,local_full_condition,A,C,normalized_A,normalized_C,normalized_pair_condition,functional_kernel_hp_kkt\n";
    for(const auto & row:profile)
    {
        const auto & value=row.as_object(); const auto & pair=value.at("pair").as_object();
        const auto print=[&](const char * key) {return value.at(key).is_null() ? std::string("nan") : j::serialize(value.at(key));};
        csv<<j::serialize(value.at("atom_id"))<<','<<print("width")<<','<<print("profile_loss")<<','
            <<print("inner_kkt")<<','<<j::serialize(value.at("profile_valid"))<<','
            <<j::serialize(value.at("profile_reason"))<<','<<print("local_full_rank")<<','
            <<print("local_full_condition")<<','<<print("A")<<','<<print("C")<<','
            <<print("normalized_A")<<','<<print("normalized_C")<<','
            <<(pair.at("normalized_pair_condition").is_null() ? "nan" : j::serialize(pair.at("normalized_pair_condition")))<<','
            <<j::serialize(value.at("functional_kernel_hp_kkt"))<<'\n';
    }
    csv.close();
    j::array points; for(const auto & row:profile) points.push_back(row);
    const auto identity=fixture.at("identity");
    std::ofstream summary(output/"target-width-profile.json");
    summary.exceptions(std::ios::failbit|std::ios::badbit);
    summary<<j::serialize(j::object{{"atom_id",atom_id},{"component_id",identity.at("component_id")},
        {"block_index",identity.at("block_index")},{"fixture_width",std::exp(eta(atom))},
        {"specified_widths",j::array{.1,.2,.3,.5,1,2,5,10,100,1000}},
        {"half_step_profile_curvature",ProfileCurvature(profile,eta(atom),.5)},
        {"profile_points",std::move(points)},
        {"scope","Conditional target width slice; other block widths held fixed and constrained A/C re-solved at each point."}})<<'\n';
    summary.close();
}
}
int main(int argc,char ** argv)
{
    try
    {
        if(argc==3) Run(argv[1],argv[2]);
        else if(argc==5 && std::string(argv[3])=="--target") RunTarget(argv[1],argv[2],argv[4]);
        else throw std::invalid_argument("Usage: joint_fullabc_width_study FIXTURE.json NEW_OUTPUT_DIR [--target ATOM_ID]");
        return 0;
    }
    catch(const std::exception & error) {std::cerr<<error.what()<<'\n'; return 1;}
}
