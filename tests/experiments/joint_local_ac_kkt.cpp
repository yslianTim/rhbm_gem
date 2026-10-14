#include "core/detail/joint_component/Numerics.hpp"
#include <boost/json.hpp>
#include <fstream>
#include <iostream>
#include <iterator>

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
    const bool solver_matches=evaluation.certificate.linear_solver_valid==original.at("valid").as_bool() &&
        evaluation.certificate.free_rank.value_or(-1)==j::value_to<int>(original.at("linear_rank")) &&
        evaluation.certificate.linear_solves.value_or(-1)==j::value_to<int>(original.at("linear_solve_count")) &&
        evaluation.certificate.linear_releases.value_or(-1)==j::value_to<int>(original.at("active_set_releases"));
    const bool reproduced=!evaluation.valid && evaluation.reason=="kkt-failed" &&
        search.stop_reason=="inner-kkt-failed" && solver_matches;
    return j::object{
        {"reproduced",reproduced},{"solver_diagnostics_match",solver_matches},
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
        {"primal_feasible",evaluation.certificate.feasible},
        {"effective_response_norm",y.norm()},{"parent_observation_scale",context.scale}};
}
}

int main(int argc,char ** argv)
{
    try
    {
        if(argc!=3 || std::string(argv[1])!="replay")
            throw std::invalid_argument("Usage: joint_local_ac_kkt replay FIXTURE.json");
        const auto result=Run(Read(argv[2])); std::cout<<j::serialize(result)<<'\n';
        return result.at("reproduced").as_bool() ? 0 : 1;
    }
    catch(const std::exception & error)
    {std::cerr<<error.what()<<'\n'; return 2;}
}
