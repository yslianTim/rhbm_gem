#include "core/detail/joint_component/Numerics.hpp"
#include <boost/json.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <numbers>
#include <numeric>

namespace {
namespace j=boost::json;
namespace n=rhbm_gem::core::joint_component;
using MP=boost::multiprecision::cpp_dec_float_50;
using V=n::Vector;
struct Kernel {double g{},k{},dg{},dk{};};
struct MpKernel {MP g{},k{},dg{},dk{};};
struct Entry {Eigen::Index row{}; std::size_t atom{}; double square{}; Kernel primary,replay; MpKernel reference;};

j::value Read(const char * path)
{
    std::ifstream input(path);
    if(!input) throw std::runtime_error(std::string("Cannot read replay fixture: ")+path);
    j::parse_options options; options.numbers=j::number_precision::precise;
    return j::parse(std::string(std::istreambuf_iterator<char>(input),{}),{},options);
}
V VectorFrom(const j::value & value)
{
    const auto & array=value.as_array(); V out(static_cast<Eigen::Index>(array.size()));
    for(std::size_t k=0;k<array.size();++k) out(static_cast<Eigen::Index>(k))=j::value_to<double>(array[k]);
    return out;
}
double Number(double value) {return std::isfinite(value) ? value : std::numeric_limits<double>::quiet_NaN();}
double MaxAbs(n::VectorRef value) {return value.size() ? value.cwiseAbs().maxCoeff() : 0.0;}
double Difference(n::VectorRef lhs,n::VectorRef rhs)
{return lhs.size()==rhs.size() && lhs.size() ? (lhs-rhs).lpNorm<Eigen::Infinity>() : 0.0;}
double Gamma(std::size_t operations)
{
    const double n=static_cast<double>(operations),ne=n*std::numeric_limits<double>::epsilon();
    return ne<1 ? ne/(1-ne) : std::numeric_limits<double>::infinity();
}
Kernel ReplayKernel(double square,double width,bool cutoff=false)
{
    Kernel out;
    if(cutoff && square>2.5*2.5) return out;
    const double radius=std::sqrt(square),width_squared=width*width;
    const double exponential=std::exp(-square/(2*width_squared));
    out.g=std::pow(2*M_PI*width*width,-1.5)*exponential;
    out.k=radius<1e-5 ? std::sqrt(2/M_PI)/width : std::erf(radius/width/std::sqrt(2.0))/radius;
    out.dg=out.g*(square/width_squared-3.0);
    out.dk=-std::sqrt(2/M_PI)/width*(square<1e-10 ? 1 : exponential);
    return out;
}
MpKernel ReferenceKernel(double square,double width)
{
    const MP s(square),b(width),cutoff("6.25"); MpKernel out;
    if(s>cutoff) return out;
    const MP pi=boost::math::constants::pi<MP>(),root2=sqrt(MP(2)),root= sqrt(s);
    const MP b2=b*b,exponential=exp(-s/(MP(2)*b2));
    const MP norm=MP(2)*pi*b2;
    out.g=exponential/(norm*sqrt(norm));
    out.dg=out.g*(s/b2-MP(3));
    const MP center=sqrt(MP(2)/pi)/b;
    if(root<MP("1e-5")) {out.k=center; out.dk=-center;}
    else {out.k=boost::math::erf(root/(root2*b))/root; out.dk=-center*exponential;}
    return out;
}
double Pairwise(const std::vector<double> & values,std::size_t begin,std::size_t end)
{
    if(begin==end) return 0;
    if(end-begin==1) return values[begin];
    const auto middle=begin+(end-begin)/2;
    return Pairwise(values,begin,middle)+Pairwise(values,middle,end);
}
double Kahan(const std::vector<double> & values)
{
    double sum{},correction{};
    for(const double value:values)
    {const double adjusted=value-correction,next=sum+adjusted; correction=(next-sum)-adjusted; sum=next;}
    return sum;
}
double Neumaier(const std::vector<double> & values)
{
    double sum{},correction{};
    for(const double value:values)
    {
        const double next=sum+value;
        correction+=std::abs(sum)>=std::abs(value) ? (sum-next)+value : (value-next)+sum;
        sum=next;
    }
    return sum+correction;
}
struct MaxError {MP absolute{},relative{};};
void Observe(MaxError & error,const MP & actual,const MP & reference)
{
    const MP absolute=abs(actual-reference),relative=absolute/max(abs(reference),MP("1e-300"));
    error.absolute=std::max(error.absolute,absolute); error.relative=std::max(error.relative,relative);
}
j::value KernelErrors(const std::vector<Entry> & entries)
{
    MaxError pg,pk,pdg,pdk,rg,rk,rdg,rdk,gr,kr,gdr,kdr;
    for(const auto & e:entries)
    {
        Observe(pg,MP(e.primary.g),e.reference.g); Observe(pk,MP(e.primary.k),e.reference.k);
        Observe(pdg,MP(e.primary.dg),e.reference.dg); Observe(pdk,MP(e.primary.dk),e.reference.dk);
        Observe(rg,MP(e.replay.g),e.reference.g); Observe(rk,MP(e.replay.k),e.reference.k);
        Observe(rdg,MP(e.replay.dg),e.reference.dg); Observe(rdk,MP(e.replay.dk),e.reference.dk);
        Observe(gr,MP(e.primary.g),MP(e.replay.g)); Observe(kr,MP(e.primary.k),MP(e.replay.k));
        Observe(gdr,MP(e.primary.dg),MP(e.replay.dg)); Observe(kdr,MP(e.primary.dk),MP(e.replay.dk));
    }
    const auto item=[](const MaxError & e) {return j::object{{"max_absolute",e.absolute.convert_to<double>()},
        {"max_relative",e.relative.convert_to<double>()}};};
    return j::object{{"gaussian_primary",item(pg)},{"charge_primary",item(pk)},
        {"gaussian_log_width_primary",item(pdg)},{"charge_log_width_primary",item(pdk)},
        {"gaussian_replay",item(rg)},{"charge_replay",item(rk)},
        {"gaussian_log_width_replay",item(rdg)},{"charge_log_width_replay",item(rdk)},
        {"primary_replay_gaussian",item(gr)},{"primary_replay_charge",item(kr)},
        {"primary_replay_gaussian_log_width",item(gdr)},{"primary_replay_charge_log_width",item(kdr)}};
}
j::value Run(const j::value & fixture)
{
    if(j::value_to<std::string>(fixture.at("format"))!="rhbm-joint-local-replay-trial-v1")
        throw std::invalid_argument("Unsupported local replay trial fixture format.");
    const auto rows=j::value_to<Eigen::Index>(fixture.at("identity").at("local_observation_count"));
    std::vector<std::vector<n::Support>> support;
    for(const auto & atom:fixture.at("local_support").as_array())
    {
        std::vector<n::Support> values;
        for(const auto & item:atom.as_array()) values.push_back({j::value_to<Eigen::Index>(item.at("row")),
            j::value_to<double>(item.at("squared_distance"))});
        support.push_back(std::move(values));
    }
    n::Domain domain{rows,std::move(support)}; const V y=VectorFrom(fixture.at("effective_response"));
    const V eta=VectorFrom(fixture.at("trial_eta")),beta=VectorFrom(fixture.at("trial_beta"));
    const auto & saved_context=fixture.at("local_numerical_context"); n::EvaluationContext context;
    context.scale=j::value_to<double>(saved_context.at("local_observation_scale"));
    context.linear.rank_relative=j::value_to<double>(saved_context.at("rank_relative"));
    context.atom_ids=n::Identities(std::vector<std::string>{}); context.row_ids=n::Identities(std::vector<std::string>{});
    const auto primary=n::EvaluateState(domain,y,eta,beta,context);
    n::ReplayDetailsForDiagnostics replay_details;
    const auto trust=n::CheckReplay(domain,y,primary,context,&replay_details);
    const V primary_prediction=(primary.x*primary.beta).eval();
    const std::size_t atoms=domain.atoms.size();
    std::vector<Entry> entries; entries.reserve(static_cast<std::size_t>(primary.x.nonZeros()));
    std::vector<std::vector<double>> primary_terms(static_cast<std::size_t>(rows)),replay_terms(static_cast<std::size_t>(rows));
    std::vector<double> primary_absolute(static_cast<std::size_t>(rows)),replay_absolute(static_cast<std::size_t>(rows));
    std::vector<MP> reference_prediction(static_cast<std::size_t>(rows),MP(0)),reference_absolute(static_cast<std::size_t>(rows),MP(0));
    std::vector<MP> primary_kernel_error(static_cast<std::size_t>(rows),MP(0)),replay_kernel_error(static_cast<std::size_t>(rows),MP(0));
    std::vector<MP> primary_rounded_sum(static_cast<std::size_t>(rows),MP(0)),replay_rounded_sum(static_cast<std::size_t>(rows),MP(0));
    std::vector<MP> primary_kernel_delta(static_cast<std::size_t>(rows),MP(0)),replay_kernel_delta(static_cast<std::size_t>(rows),MP(0));
    std::vector<MP> primary_product_delta(static_cast<std::size_t>(rows),MP(0)),replay_product_delta(static_cast<std::size_t>(rows),MP(0));
    std::vector<MP> primary_kernel_abs(static_cast<std::size_t>(rows),MP(0)),replay_kernel_abs(static_cast<std::size_t>(rows),MP(0));
    std::vector<MP> primary_product_abs(static_cast<std::size_t>(rows),MP(0)),replay_product_abs(static_cast<std::size_t>(rows),MP(0));
    std::size_t above_cutoff{};
    for(std::size_t atom=0;atom<atoms;++atom)
    {
        const double width=std::exp(eta(static_cast<Eigen::Index>(atom)));
        const double a=beta(static_cast<Eigen::Index>(2*atom)),c=beta(static_cast<Eigen::Index>(2*atom+1));
        for(const auto & item:domain.atoms[atom])
        {
            const auto basis=n::EvaluateKernel(item.square,width,2.5);
            const Kernel p{basis.gaussian,basis.charge,basis.gaussian_log_width,basis.charge_log_width};
            const auto r=ReplayKernel(item.square,width);
            const auto ref=ReferenceKernel(item.square,width); Entry entry{item.row,atom,item.square,p,r,ref}; entries.push_back(entry);
            if(item.square>2.5*2.5) ++above_cutoff;
            auto & pt=primary_terms[static_cast<std::size_t>(item.row)];
            auto & rt=replay_terms[static_cast<std::size_t>(item.row)];
            const double pag=a*p.g,pck=c*p.k,rag=a*r.g,rck=c*r.k;
            pt.push_back(pag); pt.push_back(pck); rt.push_back(rag); rt.push_back(rck);
            primary_absolute[static_cast<std::size_t>(item.row)]+=std::abs(pag)+std::abs(pck);
            replay_absolute[static_cast<std::size_t>(item.row)]+=std::abs(rag)+std::abs(rck);
            const MP exact_a=MP(a)*ref.g,exact_c=MP(c)*ref.k;
            const MP exact=exact_a+exact_c;
            const MP primary_a=MP(a)*MP(p.g),primary_c=MP(c)*MP(p.k);
            const MP replay_a=MP(a)*MP(r.g),replay_c=MP(c)*MP(r.k);
            const MP primary_rounded_a=MP(pag),primary_rounded_c=MP(pck);
            const MP replay_rounded_a=MP(rag),replay_rounded_c=MP(rck);
            const MP primary_kernel_abs_term=abs(primary_a-exact_a)+abs(primary_c-exact_c);
            const MP replay_kernel_abs_term=abs(replay_a-exact_a)+abs(replay_c-exact_c);
            const MP primary_product_abs_term=abs(primary_rounded_a-primary_a)+abs(primary_rounded_c-primary_c);
            const MP replay_product_abs_term=abs(replay_rounded_a-replay_a)+abs(replay_rounded_c-replay_c);
            reference_prediction[static_cast<std::size_t>(item.row)]+=exact;
            reference_absolute[static_cast<std::size_t>(item.row)]+=abs(exact_a)+abs(exact_c);
            const auto row=static_cast<std::size_t>(item.row);
            primary_rounded_sum[row]+=primary_rounded_a+primary_rounded_c;
            replay_rounded_sum[row]+=replay_rounded_a+replay_rounded_c;
            primary_kernel_delta[row]+=(primary_a-exact_a)+(primary_c-exact_c);
            replay_kernel_delta[row]+=(replay_a-exact_a)+(replay_c-exact_c);
            primary_product_delta[row]+=(primary_rounded_a-primary_a)+(primary_rounded_c-primary_c);
            replay_product_delta[row]+=(replay_rounded_a-replay_a)+(replay_rounded_c-replay_c);
            primary_kernel_abs[row]+=primary_kernel_abs_term; replay_kernel_abs[row]+=replay_kernel_abs_term;
            primary_product_abs[row]+=primary_product_abs_term; replay_product_abs[row]+=replay_product_abs_term;
            primary_kernel_error[row]+=primary_kernel_abs_term+primary_product_abs_term;
            replay_kernel_error[row]+=replay_kernel_abs_term+replay_product_abs_term;
        }
    }
    const auto row_count=static_cast<std::size_t>(rows);
    std::vector<double> naive_primary(row_count),kahan_primary(row_count),neumaier_primary(row_count),pairwise_primary(row_count),
        kahan_replay(row_count),neumaier_replay(row_count),pairwise_replay(row_count),row_diff(row_count);
    std::vector<MP> residual_reference(static_cast<std::size_t>(rows));
    std::vector<MP> primary_accumulation_delta(row_count),replay_accumulation_delta(row_count);
    double max_prediction_error_primary{},max_prediction_error_replay{},max_prediction_difference{},
        max_residual_formation_error_primary{},max_residual_formation_error_replay{},max_residual_error_primary{},
        max_residual_error_replay{},max_objective_error_primary{},max_objective_error_replay{};
    MP reference_objective=0; double primary_objective=0,replay_objective=0;
    const double delta=64*std::numeric_limits<double>::epsilon()*std::max(1.0,MaxAbs(primary_prediction));
    std::vector<MP> primary_bound(static_cast<std::size_t>(rows)),replay_bound(static_cast<std::size_t>(rows));
    for(Eigen::Index row=0;row<rows;++row)
    {
        const auto idx=static_cast<std::size_t>(row); const auto & pt=primary_terms[idx]; const auto & rt=replay_terms[idx];
        naive_primary[idx]=std::accumulate(pt.begin(),pt.end(),0.0);
        kahan_primary[idx]=Kahan(pt); neumaier_primary[idx]=Neumaier(pt); pairwise_primary[idx]=Pairwise(pt,0,pt.size());
        kahan_replay[idx]=Kahan(rt); neumaier_replay[idx]=Neumaier(rt); pairwise_replay[idx]=Pairwise(rt,0,rt.size());
        const MP exact_residual=reference_prediction[idx]-MP(y(row)); residual_reference[idx]=exact_residual;
        const MP primary_error=MP(primary_prediction(row))-reference_prediction[idx];
        const MP replay_error=MP(replay_details.prediction(row))-reference_prediction[idx];
        primary_accumulation_delta[idx]=MP(primary_prediction(row))-primary_rounded_sum[idx];
        replay_accumulation_delta[idx]=MP(replay_details.prediction(row))-replay_rounded_sum[idx];
        max_prediction_error_primary=std::max(max_prediction_error_primary,abs(primary_error).convert_to<double>());
        max_prediction_error_replay=std::max(max_prediction_error_replay,abs(replay_error).convert_to<double>());
        row_diff[idx]=std::abs(primary_prediction(row)-replay_details.prediction(row));
        max_prediction_difference=std::max(max_prediction_difference,row_diff[idx]);
        const MP p_resid_error=MP(primary.residual(row))- (MP(primary_prediction(row))-MP(y(row)));
        const MP r_resid_error=MP(replay_details.residual(row))-(MP(replay_details.prediction(row))-MP(y(row)));
        max_residual_formation_error_primary=std::max(max_residual_formation_error_primary,abs(p_resid_error).convert_to<double>());
        max_residual_formation_error_replay=std::max(max_residual_formation_error_replay,abs(r_resid_error).convert_to<double>());
        const MP p_resid=MP(primary.residual(row)),r_resid=MP(replay_details.residual(row));
        max_residual_error_primary=std::max(max_residual_error_primary,abs(p_resid-exact_residual).convert_to<double>());
        max_residual_error_replay=std::max(max_residual_error_replay,abs(r_resid-exact_residual).convert_to<double>());
        const MP p_obj=p_resid*p_resid/MP(2),r_obj=r_resid*r_resid/MP(2),ref_obj=exact_residual*exact_residual/MP(2);
        primary_objective+=.5*primary.residual(row)*primary.residual(row);
        replay_objective+=.5*replay_details.residual(row)*replay_details.residual(row);
        reference_objective+=ref_obj;
        max_objective_error_primary=std::max(max_objective_error_primary,abs(p_obj-ref_obj).convert_to<double>());
        max_objective_error_replay=std::max(max_objective_error_replay,abs(r_obj-ref_obj).convert_to<double>());
        const auto nterms=pt.size();
        primary_bound[idx]=primary_kernel_error[idx]+MP(Gamma(3*nterms))*MP(primary_absolute[idx]);
        replay_bound[idx]=replay_kernel_error[idx]+MP(Gamma(4*nterms))*MP(replay_absolute[idx]);
    }

    const double scale=context.scale;
    std::array<std::vector<double>,4> grouped_gradients,separate_dot_gradients;
    for(auto & value:grouped_gradients) value.assign(static_cast<std::size_t>(eta.size()),0.0);
    for(auto & value:separate_dot_gradients) value.assign(static_cast<std::size_t>(eta.size()),0.0);
    std::vector<MP> reference_gradient(atoms,MP(0));
    for(std::size_t atom=0;atom<atoms;++atom)
    {
        double primary_g_dot_p{},primary_k_dot_p{},primary_g_dot_r{},primary_k_dot_r{};
        double replay_g_dot_p{},replay_k_dot_p{},replay_g_dot_r{},replay_k_dot_r{};
        double combined_pp{},combined_pr{},combined_rp{},combined_rr{};
        const double a=beta(static_cast<Eigen::Index>(2*atom)),c=beta(static_cast<Eigen::Index>(2*atom+1));
        MP exact_sum=0;
        for(const auto & e:entries) if(e.atom==atom)
        {
            const auto row=static_cast<std::size_t>(e.row);
            const double rp=primary.residual(e.row),rr=replay_details.residual(e.row);
            primary_g_dot_p+=e.primary.dg*rp; primary_k_dot_p+=e.primary.dk*rp;
            primary_g_dot_r+=e.primary.dg*rr; primary_k_dot_r+=e.primary.dk*rr;
            replay_g_dot_p+=e.replay.dg*rp; replay_k_dot_p+=e.replay.dk*rp;
            replay_g_dot_r+=e.replay.dg*rr; replay_k_dot_r+=e.replay.dk*rr;
            combined_pp+=(a*e.primary.dg+c*e.primary.dk)*rp/scale/scale;
            combined_pr+=(a*e.primary.dg+c*e.primary.dk)*rr/scale/scale;
            combined_rp+=(a*e.replay.dg+c*e.replay.dk)*rp/scale/scale;
            combined_rr+=(a*e.replay.dg+c*e.replay.dk)*rr/scale/scale;
            const MP residual=residual_reference[row];
            exact_sum+=(MP(a)*e.reference.dg+MP(c)*e.reference.dk)*residual/(MP(scale)*MP(scale));
        }
        const auto separate=[](double a,double c,double dg,double dk,double scale)
        {double sum{}; sum+=a*dg/scale/scale; sum+=c*dk/scale/scale; return sum;};
        separate_dot_gradients[0][atom]=separate(a,c,primary_g_dot_p,primary_k_dot_p,scale);
        separate_dot_gradients[1][atom]=separate(a,c,primary_g_dot_r,primary_k_dot_r,scale);
        separate_dot_gradients[2][atom]=separate(a,c,replay_g_dot_p,replay_k_dot_p,scale);
        separate_dot_gradients[3][atom]=separate(a,c,replay_g_dot_r,replay_k_dot_r,scale);
        grouped_gradients[0][atom]=combined_pp; grouped_gradients[1][atom]=combined_pr;
        grouped_gradients[2][atom]=combined_rp; grouped_gradients[3][atom]=combined_rr;
        reference_gradient[atom]=exact_sum;
    }
    const V primary_width_gradient=primary.gradient,replay_width_gradient=replay_details.width_gradient;
    V primary_cross=Eigen::Map<V>(separate_dot_gradients[0].data(),static_cast<Eigen::Index>(atoms));
    V replay_cross=Eigen::Map<V>(separate_dot_gradients[3].data(),static_cast<Eigen::Index>(atoms));
    std::array<V,4> separate_vectors,grouped_vectors;
    for(std::size_t variant=0;variant<4;++variant)
    {
        separate_vectors[variant]=Eigen::Map<V>(separate_dot_gradients[variant].data(),static_cast<Eigen::Index>(atoms));
        grouped_vectors[variant]=Eigen::Map<V>(grouped_gradients[variant].data(),static_cast<Eigen::Index>(atoms));
    }
    V high_gradient(static_cast<Eigen::Index>(atoms));
    for(std::size_t atom=0;atom<atoms;++atom) high_gradient(static_cast<Eigen::Index>(atom))=reference_gradient[atom].convert_to<double>();

    const V replay_residual=replay_details.residual;
    const V primary_ac_gradient=(primary.x.transpose()*primary.residual).eval();
    const V replay_ac_gradient=(primary.x.transpose()*replay_residual).eval();
    V column_norms(primary.x.cols());
    for(Eigen::Index column=0;column<primary.x.cols();++column) column_norms(column)=primary.x.col(column).norm();
    const auto kkt=[&](n::VectorRef raw)
    {
        const V u=column_norms.array()*beta.array()/scale;
        const V gradient=raw.array()/column_norms.array()/scale;
        V projected=u-gradient;
        for(Eigen::Index k=0;k<beta.size();k+=2) projected(k)=std::max(0.0,projected(k));
        return Difference(u,projected);
    };
    std::array<MaxError,4> gradient_errors,separate_gradient_errors,grouped_gradient_errors;
    for(Eigen::Index atom=0;atom<high_gradient.size();++atom)
    {
        Observe(gradient_errors[0],MP(primary_width_gradient(atom)),reference_gradient[static_cast<std::size_t>(atom)]);
        Observe(gradient_errors[1],MP(replay_width_gradient(atom)),reference_gradient[static_cast<std::size_t>(atom)]);
        Observe(gradient_errors[2],MP(primary_cross(atom)),reference_gradient[static_cast<std::size_t>(atom)]);
        Observe(gradient_errors[3],MP(replay_cross(atom)),reference_gradient[static_cast<std::size_t>(atom)]);
        for(std::size_t variant=0;variant<4;++variant)
        {
            Observe(separate_gradient_errors[variant],MP(separate_vectors[variant](atom)),reference_gradient[static_cast<std::size_t>(atom)]);
            Observe(grouped_gradient_errors[variant],MP(grouped_vectors[variant](atom)),reference_gradient[static_cast<std::size_t>(atom)]);
        }
    }
    std::vector<std::size_t> top(static_cast<std::size_t>(rows)); std::iota(top.begin(),top.end(),0);
    std::sort(top.begin(),top.end(),[&](auto a,auto b){return row_diff[a]>row_diff[b];});
    j::array row_errors;
    const double objective_reference=reference_objective.convert_to<double>();
    const auto max_mp_abs=[](const std::vector<MP> & values)
    {MP result=0; for(const auto & value:values) {const MP magnitude=abs(value); if(magnitude>result) result=magnitude;}
        return result.convert_to<double>();};
    for(std::size_t rank=0;rank<std::min<std::size_t>(10,top.size());++rank)
    {
        const auto row=top[rank]; const auto index=static_cast<Eigen::Index>(row); const auto ncon=primary_terms[row].size();
        const double relative=row_diff[row]/std::max(std::abs(primary_prediction(index)),delta);
        const double condition=reference_absolute[row].convert_to<double>()/
            std::max(std::abs(reference_prediction[row].convert_to<double>()),delta);
        row_errors.push_back(j::object{{"local_row",row},{"row_id",fixture.at("row_ids").as_array()[row]},
            {"number_of_contributions",ncon},{"sum_absolute_contributions",Number(primary_absolute[row])},
            {"primary_prediction",Number(primary_prediction(index))},{"replay_prediction",Number(replay_details.prediction(index))},
            {"high_precision_prediction",reference_prediction[row].convert_to<double>()},
            {"absolute_difference",Number(row_diff[row])},{"relative_difference",Number(relative)},
            {"primary_signed_error_decomposition",j::object{
                {"kernel_evaluation",primary_kernel_delta[row].convert_to<double>()},
                {"product_rounding",primary_product_delta[row].convert_to<double>()},
                {"accumulation",primary_accumulation_delta[row].convert_to<double>()}}},
            {"replay_signed_error_decomposition",j::object{
                {"kernel_evaluation",replay_kernel_delta[row].convert_to<double>()},
                {"product_rounding",replay_product_delta[row].convert_to<double>()},
                {"accumulation",replay_accumulation_delta[row].convert_to<double>()}}},
            {"cancellation_indicator",Number(condition)}});
    }
    const double primary_error_bound=std::max_element(primary_bound.begin(),primary_bound.end())->convert_to<double>();
    const double replay_error_bound=std::max_element(replay_bound.begin(),replay_bound.end())->convert_to<double>();
    const double max_bound_violation_primary=[&]
    {
        double value{};
        for(Eigen::Index row=0;row<rows;++row)
        {
            value=std::max(value,std::max(0.0,std::abs(primary_prediction(row)-reference_prediction[static_cast<std::size_t>(row)].convert_to<double>())-
                primary_bound[static_cast<std::size_t>(row)].convert_to<double>()));
        }
        return value;
    }();
    const double max_bound_violation_replay=[&]
    {
        double value{};
        for(Eigen::Index row=0;row<rows;++row)
        {
            value=std::max(value,std::max(0.0,std::abs(replay_details.prediction(row)-reference_prediction[static_cast<std::size_t>(row)].convert_to<double>())-
                replay_bound[static_cast<std::size_t>(row)].convert_to<double>()));
        }
        return value;
    }();
    const double cutoff_next=std::nextafter(6.25,std::numeric_limits<double>::infinity());
    const auto cutoff_primary=n::EvaluateKernel(cutoff_next,.7,2.5);
    const auto cutoff_replay=ReplayKernel(cutoff_next,.7);
    const auto cutoff_endpoint=ReplayKernel(cutoff_next,.7,true);
    const auto as_vector=[](const std::vector<double> & values)
    {return V(Eigen::Map<const V>(values.data(),static_cast<Eigen::Index>(values.size())));};
    const V naive_primary_vector=as_vector(naive_primary),kahan_primary_vector=as_vector(kahan_primary),
        neumaier_primary_vector=as_vector(neumaier_primary),pairwise_primary_vector=as_vector(pairwise_primary),
        kahan_replay_vector=as_vector(kahan_replay),neumaier_replay_vector=as_vector(neumaier_replay),
        pairwise_replay_vector=as_vector(pairwise_replay);
    std::array<MaxError,9> accumulation_errors;
    const std::array<const V *,9> accumulators{&primary_prediction,&replay_details.prediction,&naive_primary_vector,
        &kahan_primary_vector,&neumaier_primary_vector,&pairwise_primary_vector,&kahan_replay_vector,
        &neumaier_replay_vector,&pairwise_replay_vector};
    for(Eigen::Index row=0;row<rows;++row)
        for(std::size_t variant=0;variant<accumulators.size();++variant)
            Observe(accumulation_errors[variant],MP((*accumulators[variant])(row)),
                reference_prediction[static_cast<std::size_t>(row)]);
    double primary_residual_roundtrip{},replay_residual_roundtrip{},primary_prediction_replay_residual{},
        replay_prediction_primary_residual{},max_residual_path_difference{};
    for(Eigen::Index row=0;row<rows;++row)
    {
        const double pp=primary_prediction(row)-y(row),rp=replay_details.prediction(row)-y(row);
        primary_residual_roundtrip=std::max(primary_residual_roundtrip,std::abs(pp-primary.residual(row)));
        replay_residual_roundtrip=std::max(replay_residual_roundtrip,std::abs(rp-replay_details.residual(row)));
        primary_prediction_replay_residual=std::max(primary_prediction_replay_residual,
            std::abs(primary_prediction(row)-(replay_details.residual(row)+y(row))));
        replay_prediction_primary_residual=std::max(replay_prediction_primary_residual,
            std::abs(replay_details.prediction(row)-(primary.residual(row)+y(row))));
        max_residual_path_difference=std::max(max_residual_path_difference,
            std::abs(primary.residual(row)-replay_details.residual(row)));
    }
    std::vector<MP> exact_norm_squared(static_cast<std::size_t>(beta.size()),MP(0));
    std::vector<MP> exact_ac_gradient(static_cast<std::size_t>(beta.size()),MP(0));
    std::vector<MP> exact_width_bound_primary(atoms,MP(0)),exact_width_bound_replay(atoms,MP(0));
    std::vector<MP> width_abs_primary(atoms,MP(0)),width_abs_replay(atoms,MP(0));
    std::vector<std::size_t> atom_terms(atoms);
    for(const auto & entry:entries)
    {
        const auto row=static_cast<std::size_t>(entry.row),atom=entry.atom;
        const MP residual=residual_reference[row],a=MP(beta(static_cast<Eigen::Index>(2*atom))),
            c=MP(beta(static_cast<Eigen::Index>(2*atom+1))),s=MP(scale);
        const std::array<MP,2> refvalues{entry.reference.g,entry.reference.k};
        const std::array<MP,2> pderivatives{MP(entry.primary.dg),MP(entry.primary.dk)},
            rderivatives{MP(entry.replay.dg),MP(entry.replay.dk)},refderivatives{entry.reference.dg,entry.reference.dk};
        const std::array<MP,2> coefficients{a,c};
        for(std::size_t kind=0;kind<2;++kind)
        {
            const auto column=2*atom+kind;
            exact_norm_squared[column]+=refvalues[kind]*refvalues[kind];
            exact_ac_gradient[column]+=refvalues[kind]*residual;
            const MP pterm=coefficients[kind]*pderivatives[kind]*MP(primary.residual(entry.row))/(s*s);
            const MP rterm=coefficients[kind]*rderivatives[kind]*MP(replay_details.residual(entry.row))/(s*s);
            width_abs_primary[atom]+=abs(pterm); width_abs_replay[atom]+=abs(rterm);
            exact_width_bound_primary[atom]+=abs(coefficients[kind]*(pderivatives[kind]-refderivatives[kind])*residual/(s*s))+
                abs(coefficients[kind]*pderivatives[kind]*(MP(primary.residual(entry.row))-residual)/(s*s));
            exact_width_bound_replay[atom]+=abs(coefficients[kind]*(rderivatives[kind]-refderivatives[kind])*residual/(s*s))+
                abs(coefficients[kind]*rderivatives[kind]*(MP(replay_details.residual(entry.row))-residual)/(s*s));
        }
        ++atom_terms[atom];
    }
    double high_precision_kkt{},max_column_norm_error{},max_raw_gradient_error_primary{},max_raw_gradient_error_replay{},
        max_primary_scaled_ac_gradient_error{},max_replay_scaled_ac_gradient_error{};
    for(Eigen::Index column=0;column<beta.size();++column)
    {
        const auto index=static_cast<std::size_t>(column),atom=index/2;
        const MP exact_norm=sqrt(exact_norm_squared[index]);
        const MP exact_scaled_gradient=exact_ac_gradient[index]/exact_norm/MP(scale);
        const MP negative_gradient=-exact_scaled_gradient;
        const MP violation=(column%2==0 && beta(column)==0) ?
            (negative_gradient>0 ? negative_gradient : MP(0)) : abs(exact_scaled_gradient);
        high_precision_kkt=std::max(high_precision_kkt,violation.convert_to<double>());
        max_column_norm_error=std::max(max_column_norm_error,
            abs(MP(column_norms(column))-exact_norm).convert_to<double>());
        max_raw_gradient_error_primary=std::max(max_raw_gradient_error_primary,
            abs(MP(primary_ac_gradient(column))-exact_ac_gradient[index]).convert_to<double>());
        max_raw_gradient_error_replay=std::max(max_raw_gradient_error_replay,
            abs(MP(replay_ac_gradient(column))-exact_ac_gradient[index]).convert_to<double>());
        const MP exact_scaled=exact_ac_gradient[index]/exact_norm/MP(scale);
        max_primary_scaled_ac_gradient_error=std::max(max_primary_scaled_ac_gradient_error,
            abs(MP(primary_ac_gradient(column))/MP(column_norms(column))/MP(scale)-exact_scaled).convert_to<double>());
        max_replay_scaled_ac_gradient_error=std::max(max_replay_scaled_ac_gradient_error,
            abs(MP(replay_ac_gradient(column))/MP(column_norms(column))/MP(scale)-exact_scaled).convert_to<double>());
        const double operation_bound=Gamma(4*std::max<std::size_t>(1,atom_terms[atom]));
        exact_width_bound_primary[atom]+=MP(operation_bound)*width_abs_primary[atom];
        exact_width_bound_replay[atom]+=MP(operation_bound)*width_abs_replay[atom];
    }
    double max_width_bound_violation_primary{},max_width_bound_violation_replay{},
        max_width_bound_primary{},max_width_bound_replay{};
    for(std::size_t atom=0;atom<atoms;++atom)
    {
        max_width_bound_primary=std::max(max_width_bound_primary,exact_width_bound_primary[atom].convert_to<double>());
        max_width_bound_replay=std::max(max_width_bound_replay,exact_width_bound_replay[atom].convert_to<double>());
        max_width_bound_violation_primary=std::max(max_width_bound_violation_primary,
            std::max(0.0,std::abs(primary_width_gradient(static_cast<Eigen::Index>(atom))-high_gradient(static_cast<Eigen::Index>(atom)))-
                exact_width_bound_primary[atom].convert_to<double>()));
        max_width_bound_violation_replay=std::max(max_width_bound_violation_replay,
            std::max(0.0,std::abs(replay_width_gradient(static_cast<Eigen::Index>(atom))-high_gradient(static_cast<Eigen::Index>(atom)))-
                exact_width_bound_replay[atom].convert_to<double>()));
    }
    double max_square_roundtrip{},max_fma_path_difference{},max_fma_error_reference{};
    for(const auto & entry:entries)
    {
        const double a=beta(static_cast<Eigen::Index>(2*entry.atom)),c=beta(static_cast<Eigen::Index>(2*entry.atom+1));
        const double ordinary=a*entry.primary.g+c*entry.primary.k;
        const double fused=std::fma(a,entry.primary.g,c*entry.primary.k);
        const MP exact=MP(a)*entry.reference.g+MP(c)*entry.reference.k;
        max_fma_path_difference=std::max(max_fma_path_difference,std::abs(fused-ordinary));
        max_fma_error_reference=std::max(max_fma_error_reference,abs(MP(fused)-exact).convert_to<double>());
        max_square_roundtrip=std::max(max_square_roundtrip,
            std::abs(std::sqrt(entry.square)*std::sqrt(entry.square)-entry.square));
    }
    std::array<MaxError,2> center_errors;
    for(const double square:{0.0,std::nextafter(1e-10,0.0),1e-10,
        std::nextafter(1e-10,std::numeric_limits<double>::infinity())})
    {
        const auto basis=n::EvaluateKernel(square,.55,2.5);
        const Kernel primary_kernel{basis.gaussian,basis.charge,basis.gaussian_log_width,basis.charge_log_width};
        const auto replay_kernel=ReplayKernel(square,.55);
        const auto reference_kernel=ReferenceKernel(square,.55);
        Observe(center_errors[0],MP(primary_kernel.k),reference_kernel.k);
        Observe(center_errors[0],MP(primary_kernel.dk),reference_kernel.dk);
        Observe(center_errors[1],MP(replay_kernel.k),reference_kernel.k);
        Observe(center_errors[1],MP(replay_kernel.dk),reference_kernel.dk);
    }
    const MP pi=boost::math::constants::pi<MP>();
    const double pi_primary=std::numbers::pi,pi_replay=M_PI;
    const auto error_json=[](const MaxError & value) {return j::object{
        {"max_absolute",value.absolute.convert_to<double>()},{"max_relative",value.relative.convert_to<double>()}};};
    j::object accumulations{{"sparse_primary",error_json(accumulation_errors[0])},
        {"compensated_replay",error_json(accumulation_errors[1])},
        {"sequential_primary_order",error_json(accumulation_errors[2])},
        {"kahan_primary_order",error_json(accumulation_errors[3])},
        {"neumaier_primary_order",error_json(accumulation_errors[4])},
        {"pairwise_primary_order",error_json(accumulation_errors[5])},
        {"kahan_replay_order",error_json(accumulation_errors[6])},
        {"neumaier_replay_order",error_json(accumulation_errors[7])},
        {"pairwise_replay_order",error_json(accumulation_errors[8])}};
    j::object crossed_gradients{{"primary_derivatives_primary_residual_separate_dots",error_json(gradient_errors[2])},
        {"primary_derivatives_replay_residual_separate_dots",error_json(separate_gradient_errors[1])},
        {"replay_derivatives_primary_residual_separate_dots",error_json(separate_gradient_errors[2])},
        {"replay_derivatives_replay_residual_separate_dots",error_json(gradient_errors[3])},
        {"primary_derivatives_primary_residual_combined_terms",error_json(grouped_gradient_errors[0])},
        {"primary_derivatives_replay_residual_combined_terms",error_json(grouped_gradient_errors[1])},
        {"replay_derivatives_primary_residual_combined_terms",error_json(grouped_gradient_errors[2])},
        {"replay_derivatives_replay_residual_combined_terms",error_json(grouped_gradient_errors[3])},
        {"residual_swap_effect_same_primary_derivatives",Number(Difference(grouped_vectors[0],grouped_vectors[1]))},
        {"derivative_swap_effect_same_primary_residual",Number(Difference(grouped_vectors[0],grouped_vectors[2]))},
        {"primary_production_vs_separate_dot_order",Number(Difference(primary_width_gradient,separate_vectors[0]))},
        {"replay_production_vs_combined_term_order",Number(Difference(replay_width_gradient,grouped_vectors[3]))},
        {"separate_vs_combined_same_kernels_residual",Number(Difference(separate_vectors[0],grouped_vectors[0]))}};
    j::object kernel_checks{{"support_entries_above_cutoff",above_cutoff},
        {"sqrt_then_square_max_absolute_error",max_square_roundtrip},
        {"pi_constants_equal",pi_primary==pi_replay},{"pi_reference",pi.convert_to<double>()},
        {"center_and_threshold_kernel_errors",j::object{{"primary",error_json(center_errors[0])},
            {"local_replay",error_json(center_errors[1])}}},
        {"nextafter_above_cutoff",cutoff_next},{"production_gaussian_above_cutoff",Number(cutoff_primary.gaussian)},
        {"production_charge_above_cutoff",Number(cutoff_primary.charge)},
        {"local_replay_gaussian_above_cutoff",Number(cutoff_replay.g)},
        {"local_replay_charge_above_cutoff",Number(cutoff_replay.k)},
        {"endpoint_replay_charge_above_cutoff",Number(cutoff_endpoint.k)},
        {"fma_vs_separate_max_contribution_difference",max_fma_path_difference},
        {"fma_max_abs_error_vs_reference",max_fma_error_reference}};
    return j::object{
        {"format","rhbm-joint-local-replay-analysis-v1"},{"identity",fixture.at("identity")},
        {"high_precision_reference",j::object{{"type","boost::multiprecision::cpp_dec_float_50"},
            {"decimal_digits",std::numeric_limits<MP>::digits10},{"long_double_digits",std::numeric_limits<long double>::digits},
            {"source_doubles_preserved",true},{"erf","boost::math::erf(cpp_dec_float_50)"}}},
        {"kernel_errors",KernelErrors(entries)},{"kernel_checks",kernel_checks},
        {"prediction",j::object{{"primary_max_abs_error_vs_reference",max_prediction_error_primary},
            {"replay_max_abs_error_vs_reference",max_prediction_error_replay},
            {"primary_replay_max_difference",max_prediction_difference},
            {"max_abs_net_error_source",j::object{
                {"primary_kernel_evaluation",max_mp_abs(primary_kernel_delta)},
                {"primary_product_rounding",max_mp_abs(primary_product_delta)},
                {"primary_accumulation",max_mp_abs(primary_accumulation_delta)},
                {"replay_kernel_evaluation",max_mp_abs(replay_kernel_delta)},
                {"replay_product_rounding",max_mp_abs(replay_product_delta)},
                {"replay_accumulation",max_mp_abs(replay_accumulation_delta)}}},
            {"max_abs_sum_of_component_errors",j::object{
                {"primary_kernel_evaluation",max_mp_abs(primary_kernel_abs)},
                {"primary_product_rounding",max_mp_abs(primary_product_abs)},
                {"replay_kernel_evaluation",max_mp_abs(replay_kernel_abs)},
                {"replay_product_rounding",max_mp_abs(replay_product_abs)}}},
            {"accumulation_variants_error_vs_reference",accumulations},{"top_differing_rows",row_errors},
            {"max_conservative_roundoff_bound_primary",primary_error_bound},
            {"max_conservative_roundoff_bound_replay",replay_error_bound},
            {"max_bound_violation_primary",max_bound_violation_primary},
            {"max_bound_violation_replay",max_bound_violation_replay},{"cancellation_delta",delta}}},
        {"residual",j::object{{"primary_max_formation_error",max_residual_formation_error_primary},
            {"replay_max_formation_error",max_residual_formation_error_replay},
            {"primary_max_abs_error_vs_reference",max_residual_error_primary},
            {"replay_max_abs_error_vs_reference",max_residual_error_replay},
            {"primary_replay_max_difference",max_residual_path_difference},
            {"primary_prediction_residual_roundtrip",primary_residual_roundtrip},
            {"replay_prediction_residual_roundtrip",replay_residual_roundtrip},
            {"primary_prediction_vs_replay_residual_roundtrip",primary_prediction_replay_residual},
            {"replay_prediction_vs_primary_residual_roundtrip",replay_prediction_primary_residual},
            {"primary_objective",Number(primary_objective)},{"replay_objective",Number(replay_objective)},
            {"high_precision_objective",objective_reference},
            {"primary_total_objective_error",abs(MP(primary_objective)-reference_objective).convert_to<double>()},
            {"replay_total_objective_error",abs(MP(replay_objective)-reference_objective).convert_to<double>()},
            {"max_row_objective_term_error_primary",max_objective_error_primary},
            {"max_row_objective_term_error_replay",max_objective_error_replay}}},
        {"width_gradient",j::object{{"primary_max_abs_error_vs_reference",gradient_errors[0].absolute.convert_to<double>()},
            {"replay_max_abs_error_vs_reference",gradient_errors[1].absolute.convert_to<double>()},
            {"primary_max_relative_error_vs_reference",gradient_errors[0].relative.convert_to<double>()},
            {"replay_max_relative_error_vs_reference",gradient_errors[1].relative.convert_to<double>()},
            {"primary_gradient_inf",Number(MaxAbs(primary_width_gradient))},
            {"replay_gradient_inf",Number(MaxAbs(replay_width_gradient))},
            {"primary_replay_difference",Number(Difference(primary_width_gradient,replay_width_gradient))},
            {"crossed_derivative_residual_variants",crossed_gradients},
            {"high_precision_gradient_inf",Number(MaxAbs(high_gradient))},
            {"max_conservative_gradient_roundoff_bound_primary",max_width_bound_primary},
            {"max_conservative_gradient_roundoff_bound_replay",max_width_bound_replay},
            {"max_gradient_bound_violation_primary",max_width_bound_violation_primary},
            {"max_gradient_bound_violation_replay",max_width_bound_violation_replay}}},
        {"ac_kkt",j::object{{"primary_scaled_kkt",Number(primary.certificate.projected_kkt)},
            {"replay_scaled_kkt",Number(replay_details.scaled_kkt)},
            {"high_precision_scaled_kkt",high_precision_kkt},
            {"primary_kkt_passed_1e-10",primary.certificate.projected_kkt<=1e-10},
            {"replay_kkt_passed_1e-10",replay_details.scaled_kkt<=1e-10},
            {"primary_replay_difference",Number(trust.kkt_difference)},
            {"raw_gradient_max_difference",Number(Difference(primary_ac_gradient,replay_ac_gradient))},
            {"primary_raw_gradient_max_error_vs_reference",max_raw_gradient_error_primary},
            {"replay_raw_gradient_max_error_vs_reference",max_raw_gradient_error_replay},
            {"scaled_gradient_max_error_primary_vs_reference",max_primary_scaled_ac_gradient_error},
            {"scaled_gradient_max_error_replay_vs_reference",max_replay_scaled_ac_gradient_error},
            {"column_norm_max_error_vs_reference",max_column_norm_error},
            {"primary_kkt_recomputed",Number(kkt(primary_ac_gradient))},
            {"replay_kkt_recomputed",Number(kkt(replay_ac_gradient))}}},
        {"trust",j::object{{"prediction_passed",trust.prediction_passed},{"gradient_passed",trust.gradient_passed},
            {"passed",trust.passed},{"reason",trust.reason}}}};
}
}

int main(int argc,char ** argv)
{
    try
    {
        if(argc!=2) throw std::invalid_argument("Usage: joint_local_replay_analysis FIXTURE.json");
        Eigen::setNbThreads(1);
        std::cout<<j::serialize(Run(Read(argv[1])))<<'\n';
        return 0;
    }
    catch(const std::exception & error)
    {std::cerr<<error.what()<<'\n'; return 2;}
}
