#include "core/detail/joint_component/FixedBBlockCoordinate.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "support/JointOperatorWorkload.hpp"
#include <boost/json.hpp>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>

namespace {
namespace n=rhbm_gem::core::joint_component;
namespace j=boost::json;
using Clock=std::chrono::steady_clock;
using Input=rhbm_gem::core::JointProblemInput;
double Seconds(Clock::time_point start) {return std::chrono::duration<double>(Clock::now()-start).count();}
n::Vector Select(n::VectorRef v,const std::vector<std::size_t> & indices)
{
    n::Vector out(static_cast<Eigen::Index>(indices.size()));
    for(std::size_t k=0;k<indices.size();++k) out(static_cast<Eigen::Index>(k))=v(static_cast<Eigen::Index>(indices[k]));
    return out;
}
j::array BlockTelemetry(const std::vector<n::FixedBBlockRecord> & records)
{
    j::array out;
    for(const auto & record:records)
        out.push_back({{"sweep",record.sweep},{"block",record.block},
            {"objective_before",record.objective_before},{"objective_after",record.objective_after},
            {"objective_reduction",record.objective_reduction},
            {"local_objective_before",record.local_objective_before},
            {"local_objective_after",record.local_objective_after},
            {"local_objective_reduction",record.local_objective_reduction},
            {"local_kkt",record.local_kkt},{"local_feasibility",record.local_feasibility},
            {"global_kkt_before",record.global_kkt_before},{"global_kkt_after",record.global_kkt_after},
            {"max_scaled_ac_change",record.max_scaled_ac_change},{"raw_ac_change",record.raw_ac_change},
            {"solver_beta_equal_old",record.solver_beta_equal_old},
            {"solver_beta_norm_difference",record.solver_beta_norm_difference},
            {"accepted",record.accepted},{"status",record.status},{"reason",record.reason},
            {"global_replay_delta",record.global_replay_delta},
            {"local_global_delta_error",record.local_global_delta_error},
            {"objective_replay_enclosure",record.objective_replay_enclosure},
            {"factor_seconds",record.factor_seconds},{"linear_solves",record.linear_solves},
            {"certificate",{{"available",record.local_certificate_available},
                {"feasible",record.local_certificate_feasible},
                {"kkt_passed",record.local_certificate_kkt_passed},
                {"projected_kkt",record.local_certificate_projected_kkt},
                {"rss",record.local_certificate_rss},{"objective",record.local_certificate_objective},
                {"linear_solves",record.linear_solves},{"rank",record.rank},
                {"active_atoms",record.active_atoms},{"free_columns",record.free_columns}}}});
    return out;
}
struct Outcome
{
    std::string topology,order;
    int atoms{};
    std::size_t blocks{};
    double global_objective{},block_objective{},prediction_inf{},residual_inf{},ac_scaled_inf{},ac_raw_inf{};
    double global_kkt{},block_kkt{},global_feasibility{},block_feasibility{},search_seconds{},initial_objective{};
    std::size_t sweeps{},sweeps_to_global_kkt{},sweeps_to_objective_parity{},active_face_mismatches{},global_free_columns{},free_columns{};
    bool success{},objective_parity{},prediction_parity{},ac_parity{},active_face_equal{};
    std::string reason;
    j::object individual;
    std::vector<double> sweep_objectives,sweep_kkt;
    std::vector<double> block_scaled_parameters;
    j::array block_telemetry;
};
Outcome Run(const std::string & topology,int atoms,const std::string & order)
{
    auto input=std::make_shared<Input>(second_stage_test::OperatorWorkload(topology,atoms));
    const auto & layout=n::BuildParameterLayout(*input); const n::Domain parent_domain(input);
    const auto profile_domain=n::ProfileDomain(parent_domain,layout);
    const auto parent=n::CreateContext(input); const auto context=n::ProfileContext(parent,layout,parent_domain.rows);
    const n::Vector observations=Eigen::Map<const n::Vector>(input->observations.data(),static_cast<Eigen::Index>(input->observations.size()));
    const n::Vector eta=n::Vector::Constant(static_cast<Eigen::Index>(input->atom_ids.size()),std::log(.5));
    const auto profile_y=Select(observations,layout.informative_rows),profile_eta=Select(eta,layout.full_atoms);
    const auto global=n::EvaluateProfile(profile_domain,profile_y,profile_eta,false,&context);
    if(!global.valid) throw std::runtime_error("Global fixed-B profile failed: "+global.reason);
    n::FixedBBlockPolicy policy; policy.core_atoms=128; policy.maximum_sweeps=12; policy.capture_diagnostics=true;
    policy.order=order=="forward" ? n::FixedBBlockOrder::Forward : n::FixedBBlockOrder::Reverse;
    const auto partition=n::BuildStructuralBlockPartition(*input,layout,policy.core_atoms);
    const auto search_started=Clock::now();
    const auto block=n::SearchFixedBBlocks(*input,layout,observations,eta,context,policy);
    const double search_seconds=Seconds(search_started);
    const auto global_prediction=(global.x*global.beta).eval();
    const auto global_residual=(global_prediction-profile_y).eval();
    n::Vector block_prediction(static_cast<Eigen::Index>(layout.informative_rows.size())),block_residual(block_prediction.size());
    for(std::size_t k=0;k<layout.informative_rows.size();++k)
    {
        const auto row=static_cast<Eigen::Index>(layout.informative_rows[k]);
        block_prediction(static_cast<Eigen::Index>(k))=block.state.prediction(row);
        block_residual(static_cast<Eigen::Index>(k))=block.state.residual(row);
    }
    double scaled{},raw{}; std::size_t active_mismatches{},global_free_columns{},block_free_columns{};
    std::vector<double> block_scaled_parameters; block_scaled_parameters.reserve(static_cast<std::size_t>(global.beta.size()));
    for(std::size_t k=0;k<layout.full_atoms.size();++k)
    {
        const auto atom=static_cast<Eigen::Index>(layout.full_atoms[k]);
        for(Eigen::Index kind=0;kind<2;++kind)
        {
            const auto local=2*static_cast<Eigen::Index>(k)+kind,global_column=2*atom+kind;
            const double difference=block.state.beta(global_column)-global.beta(local);
            raw=std::max(raw,std::abs(difference));
            const double parameter_scale=global.x.col(local).norm()/context.scale;
            scaled=std::max(scaled,std::abs(difference)*parameter_scale);
            block_scaled_parameters.push_back(block.state.beta(global_column)*parameter_scale);
            if(kind==0)
            {
                active_mismatches+=static_cast<std::size_t>((block.state.beta(global_column)==0)!=(global.beta(local)==0));
                global_free_columns+=static_cast<std::size_t>(global.beta(local)>0);
                block_free_columns+=static_cast<std::size_t>(block.state.beta(global_column)>0);
            }
            else {++global_free_columns; ++block_free_columns;}
        }
    }
    Outcome out; out.topology=topology; out.order=order; out.atoms=atoms;
    out.blocks=partition.cores.size();
    out.global_objective=global.certificate.objective/(context.scale*context.scale); out.block_objective=block.state.objective;
    out.prediction_inf=(block_prediction-global_prediction).lpNorm<Eigen::Infinity>();
    out.residual_inf=(block_residual-global_residual).lpNorm<Eigen::Infinity>();
    out.ac_scaled_inf=scaled; out.ac_raw_inf=raw; out.global_kkt=global.certificate.projected_kkt;
    out.block_kkt=block.sweeps.empty() ? std::numeric_limits<double>::infinity() : block.sweeps.back().global_ac_kkt;
    out.global_feasibility=global.certificate.feasible ? 0 : 1;
    out.block_feasibility=block.sweeps.empty() ? std::numeric_limits<double>::infinity() : block.sweeps.back().global_a_feasibility;
    out.search_seconds=search_seconds; out.sweeps=block.sweeps.size(); out.sweeps_to_global_kkt=block.sweeps_to_global_kkt;
    out.success=block.success; out.reason=block.reason; out.active_face_mismatches=active_mismatches;
    out.active_face_equal=active_mismatches==0; out.global_free_columns=global_free_columns;
    out.free_columns=block_free_columns; out.block_scaled_parameters=std::move(block_scaled_parameters);
    out.block_telemetry=BlockTelemetry(block.blocks);
    out.initial_objective=.5*profile_y.squaredNorm()/(context.scale*context.scale);
    out.objective_parity=std::abs(out.block_objective-out.global_objective)<=1e-12+2e-12*std::abs(out.global_objective);
    out.prediction_parity=out.prediction_inf<=2e-12+2e-13*std::max(1.0,global_prediction.cwiseAbs().maxCoeff());
    out.ac_parity=out.ac_scaled_inf<=1e-8;
    for(std::size_t k=0;k<block.sweeps.size();++k)
    {
        out.sweep_objectives.push_back(block.sweeps[k].objective_after);
        out.sweep_kkt.push_back(block.sweeps[k].global_ac_kkt);
        if(out.sweeps_to_objective_parity==0 && std::abs(block.sweeps[k].objective_after-out.global_objective)<=
            1e-12+2e-12*std::abs(out.global_objective)) out.sweeps_to_objective_parity=k+1;
    }
    out.individual={{"topology",topology},{"atoms",atoms},{"blocks",out.blocks},{"order",order},
        {"core_atoms",policy.core_atoms},{"maximum_sweeps",policy.maximum_sweeps},{"fixed_B",.5},
        {"status",block.success ? "passed" : "failed"},{"reason",block.reason},
        {"global",{{"objective",out.global_objective},{"projected_kkt",out.global_kkt},
            {"a_feasibility",out.global_feasibility},{"active_A",global.certificate.active_atoms.size()},
            {"free_columns",global_free_columns}}},
        {"block",{{"objective",out.block_objective},{"global_projected_kkt",out.block_kkt},
            {"a_feasibility",out.block_feasibility},{"sweeps",out.sweeps},
            {"sweeps_to_global_kkt",out.sweeps_to_global_kkt},
            {"sweeps_to_objective_parity",out.sweeps_to_objective_parity}}},
        {"comparison",{{"objective_difference",std::abs(out.block_objective-out.global_objective)},
            {"prediction_inf_difference",out.prediction_inf},{"residual_inf_difference",out.residual_inf},
            {"ac_scaled_inf_difference",out.ac_scaled_inf},{"ac_raw_inf_difference",out.ac_raw_inf},
            {"active_face_mismatches",active_mismatches},{"objective_parity",out.objective_parity},
            {"prediction_parity",out.prediction_parity},{"ac_parity",out.ac_parity}}},
        {"initial_objective",out.initial_objective},{"search_seconds",search_seconds},
        {"block_telemetry",out.block_telemetry},
        {"sweep_objective",[&] {j::array a; for(auto v:out.sweep_objectives) a.push_back(v); return a;}()},
        {"sweep_global_kkt",[&] {j::array a; for(auto v:out.sweep_kkt) a.push_back(v); return a;}()},
        {"diagnostics_captured",true}};
    return out;
}
void Write(const std::filesystem::path & path,const j::value & value)
{std::ofstream out(path); if(!out) throw std::runtime_error("Could not open campaign output."); out<<j::serialize(value)<<'\n';}
j::object Coupling(const std::string & topology,int atoms)
{
    auto input=std::make_shared<Input>(second_stage_test::OperatorWorkload(topology,atoms));
    const auto & layout=n::BuildParameterLayout(*input); const n::Domain parent_domain(input);
    const auto profile_domain=n::ProfileDomain(parent_domain,layout);
    const auto parent=n::CreateContext(input); const auto context=n::ProfileContext(parent,layout,parent_domain.rows);
    const n::Vector observations=Eigen::Map<const n::Vector>(input->observations.data(),static_cast<Eigen::Index>(input->observations.size()));
    const n::Vector eta=n::Vector::Constant(static_cast<Eigen::Index>(input->atom_ids.size()),std::log(.5));
    const auto profile_y=Select(observations,layout.informative_rows),profile_eta=Select(eta,layout.full_atoms);
    const auto global=n::EvaluateProfile(profile_domain,profile_y,profile_eta,false,&context);
    if(!global.valid) throw std::runtime_error("Global fixed-B profile failed: "+global.reason);
    const auto partition=n::BuildStructuralBlockPartition(*input,layout,128);
    if(partition.cores.size()!=2) throw std::runtime_error("Coupling diagnostic requires exactly two blocks.");

    std::vector<Eigen::Index> first_columns,second_columns;
    std::vector<Eigen::Index> atom_column(input->atom_ids.size(),-1);
    for(std::size_t k=0;k<layout.full_atoms.size();++k)
        atom_column.at(layout.full_atoms[k])=static_cast<Eigen::Index>(k);
    const std::array<const n::StructuralCore *,2> cores{&partition.cores[0],&partition.cores[1]};
    for(std::size_t block=0;block<cores.size();++block)
        for(auto atom:cores[block]->atoms)
        {
            const auto column=atom_column.at(static_cast<std::size_t>(atom));
            auto & columns=block==0 ? first_columns : second_columns;
            columns.push_back(2*column); columns.push_back(2*column+1);
        }

    n::Sparse gram_sparse=global.x.transpose()*global.x;
    n::Matrix gram(gram_sparse);
    n::Vector column_scale=gram.diagonal().cwiseSqrt();
    for(Eigen::Index row=0;row<gram.rows();++row)
        for(Eigen::Index column=0;column<gram.cols();++column)
            gram(row,column)/=column_scale(row)*column_scale(column);
    const auto n1=static_cast<Eigen::Index>(first_columns.size()),n2=static_cast<Eigen::Index>(second_columns.size());
    n::Matrix h11(n1,n1),h12(n1,n2),h21(n2,n1),h22(n2,n2);
    for(Eigen::Index row=0;row<n1;++row) for(Eigen::Index column=0;column<n1;++column)
        h11(row,column)=gram(first_columns[static_cast<std::size_t>(row)],first_columns[static_cast<std::size_t>(column)]);
    for(Eigen::Index row=0;row<n1;++row) for(Eigen::Index column=0;column<n2;++column)
        h12(row,column)=gram(first_columns[static_cast<std::size_t>(row)],second_columns[static_cast<std::size_t>(column)]);
    for(Eigen::Index row=0;row<n2;++row) for(Eigen::Index column=0;column<n1;++column)
        h21(row,column)=gram(second_columns[static_cast<std::size_t>(row)],first_columns[static_cast<std::size_t>(column)]);
    for(Eigen::Index row=0;row<n2;++row) for(Eigen::Index column=0;column<n2;++column)
        h22(row,column)=gram(second_columns[static_cast<std::size_t>(row)],second_columns[static_cast<std::size_t>(column)]);
    Eigen::LLT<n::Matrix> llt1(h11),llt2(h22);
    if(llt1.info()!=Eigen::Success || llt2.info()!=Eigen::Success)
        throw std::runtime_error("Block Gram factorization failed in coupling diagnostic.");
    const n::Matrix inverse_lower1=llt1.matrixL().solve(n::Matrix::Identity(n1,n1));
    const n::Matrix inverse_lower2=llt2.matrixL().solve(n::Matrix::Identity(n2,n2));
    const n::Matrix normalized_cross=(inverse_lower1*h12*inverse_lower2.transpose()).eval();
    Eigen::JacobiSVD<n::Matrix> coupling_svd(normalized_cross);
    const double canonical_coupling=coupling_svd.singularValues()(0);
    const n::Matrix gauss_seidel=llt2.solve(h21*llt1.solve(h12));
    Eigen::EigenSolver<n::Matrix> eigenvalues(gauss_seidel,false);
    if(eigenvalues.info()!=Eigen::Success) throw std::runtime_error("Gauss-Seidel spectrum failed.");
    double theoretical_radius{};
    for(Eigen::Index k=0;k<eigenvalues.eigenvalues().size();++k)
        theoretical_radius=std::max(theoretical_radius,std::abs(eigenvalues.eigenvalues()(k)));

    std::vector<std::uint8_t> row_owners(input->observations.size());
    for(std::size_t block=0;block<cores.size();++block)
        for(auto atom:cores[block]->atoms) for(const auto & support:input->support.at(static_cast<std::size_t>(atom)))
            if(std::binary_search(layout.informative_rows.begin(),layout.informative_rows.end(),support.row))
                row_owners[support.row]|=static_cast<std::uint8_t>(1u<<block);
    std::size_t touched_rows{},boundary_rows{};
    for(auto owner:row_owners) {touched_rows+=owner!=0; boundary_rows+=owner==3;}

    n::FixedBBlockPolicy policy; policy.core_atoms=128; policy.maximum_sweeps=7;
    const auto block=n::SearchFixedBBlocks(*input,layout,observations,eta,context,policy);
    std::vector<double> kkt_ratios;
    for(std::size_t k=1;k<block.sweeps.size();++k)
        if(block.sweeps[k-1].global_ac_kkt>0)
            kkt_ratios.push_back(block.sweeps[k].global_ac_kkt/block.sweeps[k-1].global_ac_kkt);
    if(kkt_ratios.empty()) throw std::runtime_error("No early KKT ratios available.");
    auto ordered_ratios=kkt_ratios; std::sort(ordered_ratios.begin(),ordered_ratios.end());
    const double median=ordered_ratios.size()%2 ? ordered_ratios[ordered_ratios.size()/2] :
        .5*(ordered_ratios[ordered_ratios.size()/2-1]+ordered_ratios[ordered_ratios.size()/2]);
    j::array ratios; for(auto ratio:kkt_ratios) ratios.push_back(ratio);
    const auto free_columns=static_cast<std::size_t>(global.beta.size())-global.certificate.active_atoms.size();
    return {{"topology",topology},{"atoms",atoms},{"blocks",partition.cores.size()},
        {"boundary_rows",boundary_rows},{"touched_rows",touched_rows},
        {"boundary_row_fraction",static_cast<double>(boundary_rows)/static_cast<double>(touched_rows)},
        {"canonical_block_coupling",canonical_coupling},{"theoretical_gs_radius",theoretical_radius},
        {"canonical_coupling_squared",canonical_coupling*canonical_coupling},
        {"radius_squared_difference",std::abs(theoretical_radius-canonical_coupling*canonical_coupling)},
        {"early_sweep_kkt_ratios",ratios},{"median_early_contraction",median},
        {"global_active_atoms",global.certificate.active_atoms.size()},{"global_free_columns",free_columns},
        {"global_kkt",global.certificate.projected_kkt},{"observed_sweeps",block.sweeps.size()}};
}
void RunCouplingCampaign(const std::filesystem::path & root)
{
    std::filesystem::create_directories(root);
    j::array cases; cases.push_back(Coupling("chain",256)); cases.push_back(Coupling("cube",256));
    Write(root/"coupling.json",j::object{{"cases",cases}});
    Write(root/"campaign-manifest.json",j::object{{"phase","F1.5 fixed-B block coupling diagnostics"},
        {"topologies",{"chain","cube"}},{"atoms",256},{"core_atoms",128},
        {"local_solver","existing constrained SolveLinear"},{"block_order","forward"},
        {"maximum_kkt_sweeps",7},{"global_kkt_tolerance",1e-10},
        {"normalization","column-scaled Gram; canonical metric invariant to per-column scaling"}});
}
}
int main(int argc,char ** argv)
{
    try {
        if(argc==3 && std::string(argv[1])=="--coupling")
        {RunCouplingCampaign(argv[2]); std::cout<<"fixed-B coupling diagnostics complete\n"; return 0;}
        if(argc!=2) throw std::invalid_argument("Usage: joint_fixed_b_block_experiment OUTPUT_DIR");
        const std::filesystem::path root(argv[1]); const auto individuals=root/"individual-results";
        std::filesystem::create_directories(individuals);
        j::array summary; std::vector<Outcome> outcomes;
        const std::vector<std::pair<std::string,int>> cases{{"chain",256},{"cube",256}};
        for(const auto & [topology,atoms]:cases) for(const std::string order:{"forward","reverse"})
        {
            auto result=Run(topology,atoms,order);
            const auto name=topology+"-"+std::to_string(atoms)+"-"+order+".json";
            Write(individuals/name,result.individual);
            j::object row{{"topology",topology},{"atoms",atoms},{"order",order},{"blocks",result.blocks},
                {"sweeps",result.sweeps},{"sweeps_to_global_kkt",result.sweeps_to_global_kkt},
                {"sweeps_to_objective_parity",result.sweeps_to_objective_parity},
                {"initial_objective",result.initial_objective},{"global_objective",result.global_objective},
                {"block_objective",result.block_objective},{"prediction_inf_difference",result.prediction_inf},
                {"residual_inf_difference",result.residual_inf},{"ac_scaled_inf_difference",result.ac_scaled_inf},
                {"ac_raw_inf_difference",result.ac_raw_inf},{"global_kkt",result.global_kkt},{"block_kkt",result.block_kkt},
                {"global_a_feasibility",result.global_feasibility},{"block_a_feasibility",result.block_feasibility},
                {"active_face_equal",result.active_face_equal},{"active_face_mismatches",result.active_face_mismatches},
                {"global_free_columns",result.global_free_columns},{"block_free_columns",result.free_columns},
                {"objective_parity",result.objective_parity},{"prediction_parity",result.prediction_parity},
                {"ac_parity",result.ac_parity},{"success",result.success},{"reason",result.reason},
                {"search_seconds",result.search_seconds}};
            summary.push_back(row); outcomes.push_back(std::move(result));
        }
        std::filesystem::create_directories(root);
        Write(root/"summary.json",j::object{{"cases",summary}});
        std::ofstream csv(root/"summary.csv");
        csv<<"topology,atoms,order,blocks,sweeps,sweeps_to_global_kkt,sweeps_to_objective_parity,global_objective,block_objective,ac_scaled_inf_difference,prediction_inf_difference,global_kkt,block_kkt,active_face_equal,global_free_columns,block_free_columns,success,reason,search_seconds\n";
        for(const auto & row:summary) csv<<j::value_to<std::string>(row.at("topology"))<<','<<j::value_to<int>(row.at("atoms"))<<','
            <<j::value_to<std::string>(row.at("order"))<<','<<j::value_to<std::size_t>(row.at("blocks"))<<','
            <<j::value_to<std::size_t>(row.at("sweeps"))<<','<<j::value_to<std::size_t>(row.at("sweeps_to_global_kkt"))<<','
            <<j::value_to<std::size_t>(row.at("sweeps_to_objective_parity"))<<','<<j::value_to<double>(row.at("global_objective"))<<','
            <<j::value_to<double>(row.at("block_objective"))<<','<<j::value_to<double>(row.at("ac_scaled_inf_difference"))<<','
            <<j::value_to<double>(row.at("prediction_inf_difference"))<<','<<j::value_to<double>(row.at("global_kkt"))<<','
            <<j::value_to<double>(row.at("block_kkt"))<<','<<j::value_to<bool>(row.at("active_face_equal"))<<','
            <<j::value_to<std::size_t>(row.at("global_free_columns"))<<','<<j::value_to<std::size_t>(row.at("block_free_columns"))<<','
            <<j::value_to<bool>(row.at("success"))<<','
            <<j::value_to<std::string>(row.at("reason"))<<','<<j::value_to<double>(row.at("search_seconds"))<<'\n';
        csv.close();
        bool parity=true,order_parity=true; j::array order_results,parity_failures,order_failures;
        for(const auto & item:outcomes) parity &= item.success && item.objective_parity &&
            item.ac_parity && item.active_face_equal && item.global_feasibility==0 && item.block_feasibility==0;
        for(const auto & item:outcomes) if(!(item.success && item.objective_parity && item.ac_parity &&
            item.active_face_equal && item.global_feasibility==0 && item.block_feasibility==0))
            parity_failures.push_back({{"topology",item.topology},{"atoms",item.atoms},{"order",item.order},
                {"reason",item.reason},{"global_kkt",item.global_kkt},{"block_kkt",item.block_kkt},
                {"ac_scaled_inf_difference",item.ac_scaled_inf},{"prediction_inf_difference",item.prediction_inf}});
        for(const auto & [topology,atoms]:cases)
        {
            const auto forward=std::find_if(outcomes.begin(),outcomes.end(),[&](const auto & r){return r.topology==topology && r.atoms==atoms && r.order=="forward";});
            const auto reverse=std::find_if(outcomes.begin(),outcomes.end(),[&](const auto & r){return r.topology==topology && r.atoms==atoms && r.order=="reverse";});
            const double objective_difference=std::abs(forward->block_objective-reverse->block_objective);
            double parameter_difference{};
            for(std::size_t k=0;k<std::min(forward->block_scaled_parameters.size(),reverse->block_scaled_parameters.size());++k)
                parameter_difference=std::max(parameter_difference,
                    std::abs(forward->block_scaled_parameters[k]-reverse->block_scaled_parameters[k]));
            const bool same_order_endpoint=forward->success && reverse->success &&
                objective_difference<=1e-12+2e-12*std::max(std::abs(forward->global_objective),std::abs(reverse->global_objective)) &&
                parameter_difference<=1e-8;
            order_parity &= same_order_endpoint;
            if(!same_order_endpoint) order_failures.push_back({{"topology",topology},{"atoms",atoms},
                {"objective_difference",objective_difference},{"scaled_ac_parameter_difference",parameter_difference},
                {"forward_success",forward->success},{"reverse_success",reverse->success}});
            order_results.push_back({{"topology",topology},{"atoms",atoms},
                {"forward_objective",forward->block_objective},{"reverse_objective",reverse->block_objective},
                {"objective_difference",objective_difference},{"ac_scaled_inf_difference",parameter_difference},
                {"sweep_count_difference",static_cast<long long>(forward->sweeps)-static_cast<long long>(reverse->sweeps)},
                {"forward_sweeps",forward->sweeps},{"reverse_sweeps",reverse->sweeps},
                {"convergence_parity",forward->success==reverse->success},{"passed",same_order_endpoint}});
        }
        Write(root/"analysis.json",j::object{{"fixed_b_convex_parity",parity ? "passed" : "failed"},
            {"forward_reverse_order_parity",order_parity ? "passed" : "failed"},
            {"parameter_scale_aware_tolerance",1e-8},{"global_kkt_tolerance",1e-10},
            {"convex_parity_failures",parity_failures},{"order_parity_failures",order_failures},
            {"prediction_comparison_role","diagnostic only; replay tolerance applies to cache-versus-full-replay checks"},
            {"objective_replay_tolerance","1e-12 + 2e-12 * abs(global_objective)"},
            {"prediction_replay_tolerance","2e-12 + 2e-13 * max(1, abs(global_prediction))"},
            {"order_comparisons",order_results}});
        Write(root/"campaign-manifest.json",j::object{{"phase","F1.5 fixed-B convergence-floor attribution"},
            {"topologies",{"chain","cube"}},{"atom_counts",{256}},
            {"fixed_B",.5},{"initial_A",0},{"initial_C",0},{"core_atoms",128},
            {"local_solver","existing constrained SolveLinear"},{"orderings",{"forward","reverse"}},
            {"maximum_sweeps",12},{"diagnostics","per-block local/global attribution; observer only"},
            {"normalization","parent scale max(1, norm(y))"},
            {"statistical_objective","unchanged FullABC informative-row profile objective"}});
        std::cout<<"fixed-B parity="<<(parity?"passed":"failed")
            <<" order parity="<<(order_parity?"passed":"failed")<<" results="<<outcomes.size()<<'\n';
        return parity && order_parity ? 0 : 1;
    } catch(const std::exception & error) {std::cerr<<error.what()<<'\n'; return 1;}
}
