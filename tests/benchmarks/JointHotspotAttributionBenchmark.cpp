#include "core/detail/joint_component/FixedNeighborBlockCoordinate.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/LinearSolve.hpp"
#include "core/detail/joint_component/TiledDerivative.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include "support/JointTestNumerics.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <boost/hash2/sha2.hpp>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <numbers>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <vector>

namespace {
namespace joint=rhbm_gem::core::joint_component;
using Clock=std::chrono::steady_clock;

struct Sample
{
    joint::FixedNeighborWork::CandidateReplayBreakdown replay;
    joint::ProfileSearchWork profile;
    std::vector<joint::TiledQrTileCensus> projected_width_census;
    double wall_seconds{},candidate_copy_seconds{},candidate_replay_seconds{},cache_update_seconds{},
        sweep_replay_seconds{},replay_audit_seconds{};
    std::size_t peak_rss_bytes{},unchanged_fast_paths{},full_replay_fallbacks{},sweeps{},accepted_updates{},
        active_width_bounds{},candidate_events{};
    double objective{},ac_kkt{},width_kkt{};
    std::string eta_digest,beta_digest,candidate_decision_digest,stop_reason;
    bool qualified{},independent_replay_passed{};
};

double Seconds(Clock::time_point start)
{return std::chrono::duration<double>(Clock::now()-start).count();}

std::size_t PeakRssBytes()
{
    rusage usage{};
    if(getrusage(RUSAGE_SELF,&usage)!=0) throw std::runtime_error("getrusage failed");
#ifdef __APPLE__
    return static_cast<std::size_t>(usage.ru_maxrss);
#else
    return static_cast<std::size_t>(usage.ru_maxrss)*1024;
#endif
}

std::string Digest(const Eigen::VectorXd & values)
{
    boost::hash2::sha2_256 hash;
    const auto add=[&hash](std::uint64_t value) {
        std::array<unsigned char,8> bytes{};
        for(unsigned k=0;k<8;++k) bytes[k]=static_cast<unsigned char>((value>>(8*k))&255);
        hash.update(bytes.data(),bytes.size());
    };
    add(static_cast<std::uint64_t>(values.size()));
    for(Eigen::Index i=0;i<values.size();++i) add(std::bit_cast<std::uint64_t>(values(i)));
    return boost::hash2::to_string(hash.result());
}

std::string Digest(const std::vector<std::uint64_t> & values)
{
    boost::hash2::sha2_256 hash;
    const auto add=[&hash](std::uint64_t value) {
        std::array<unsigned char,8> bytes{};
        for(unsigned k=0;k<8;++k) bytes[k]=static_cast<unsigned char>((value>>(8*k))&255);
        hash.update(bytes.data(),bytes.size());
    };
    add(static_cast<std::uint64_t>(values.size()));
    for(const auto value:values) add(value);
    return boost::hash2::to_string(hash.result());
}

double Median(std::vector<double> values)
{
    std::sort(values.begin(),values.end());
    const auto middle=values.size()/2;
    return values.size()%2 ? values[middle] : .5*(values[middle-1]+values[middle]);
}

template<class Getter>
void PrintTime(const char * name,const std::vector<Sample> & samples,Getter getter)
{
    std::vector<double> values;
    for(const auto & sample:samples) values.push_back(getter(sample));
    const auto [minimum,maximum]=std::minmax_element(values.begin(),values.end());
    std::cout<<' '<<name<<"_median_s="<<Median(values)<<' '<<name<<"_range_s=["<<*minimum<<','<<*maximum<<']';
}

double TileTime(const std::vector<joint::TiledQrTileCensus> & tiles,joint::TiledQrStage stage,
    double joint::TiledQrTileCensus::* member)
{
    double total{};
    for(const auto & tile:tiles) if(tile.stage==stage) total+=tile.*member;
    return total;
}

void PrintTileCensus(const char * label,const std::vector<joint::TiledQrTileCensus> & tiles,
    joint::TiledQrStage stage)
{
    std::map<std::pair<std::size_t,std::size_t>,std::size_t> shapes;
    std::size_t count{},rows{},active_rows{},structural{},numeric{},zero_rhs{},temporary_peak{},dense_cells{},
        design_cells{},derivative_cells{},active_columns{};
    std::size_t design_structural{},design_numeric{},derivative_structural{},derivative_numeric{},
        design_active{},derivative_active{},both_zero{},corrected_numeric{},corrected_active{},projection_active{};
    std::array<std::size_t,5> structural_bins{},exact_bins{},occupancy_bins{};
    for(const auto & tile:tiles) if(tile.stage==stage)
    {
        count+=tile.tile_count;
        rows+=tile.input_rows;
        active_rows+=tile.nonzero_design_rows;
        structural+=tile.structural_nonzeros;
        numeric+=tile.exact_nonzeros;
        zero_rhs+=tile.zero_design_nonzero_rhs_rows;
        active_columns+=tile.exact_nonzero_columns;
        for(std::size_t i=0;i<structural_bins.size();++i)
        {
            structural_bins[i]+=tile.structural_row_ratio_bins[i];
            exact_bins[i]+=tile.exact_row_ratio_bins[i];
            occupancy_bins[i]+=tile.dense_occupancy_bins[i];
        }
        dense_cells+=tile.input_rows*tile.columns;
        design_cells+=tile.input_rows*tile.columns;
        derivative_cells+=tile.input_rows*tile.rhs_columns;
        design_structural+=tile.design_structural_nonzeros;
        design_numeric+=tile.design_exact_nonzeros;
        derivative_structural+=tile.derivative_structural_nonzeros;
        derivative_numeric+=tile.derivative_exact_nonzeros;
        design_active+=tile.design_nonzero_rows;
        derivative_active+=tile.derivative_nonzero_rows;
        both_zero+=tile.both_zero_rows;
        corrected_numeric+=tile.corrected_jacobian_nonzeros;
        corrected_active+=tile.corrected_jacobian_nonzero_rows;
        projection_active+=tile.projection_nonzero_rows;
        temporary_peak=std::max(temporary_peak,tile.estimated_peak_temporary_bytes);
        shapes[{tile.factor_rows,tile.columns}]+=tile.tile_count;
    }
    if(!count) return;
    const auto ratio=[](std::size_t numerator,std::size_t denominator) {
        return denominator ? static_cast<double>(numerator)/static_cast<double>(denominator) : 0.0;
    };
    double occupancy=ratio(numeric,dense_cells);
    if(stage==joint::TiledQrStage::DerivativePrepare)
        occupancy=ratio(design_numeric+derivative_numeric,design_cells+derivative_cells);
    if(stage==joint::TiledQrStage::DerivativeReduce) occupancy=ratio(corrected_numeric,dense_cells);
    if(stage==joint::TiledQrStage::ProjectedWidth)
    {
        corrected_numeric=std::accumulate(tiles.begin(),tiles.end(),std::size_t{},[&](std::size_t total,const auto & tile) {
            return total+(tile.stage==stage ? tile.projection_nonzeros : 0);
        });
        occupancy=ratio(corrected_numeric,dense_cells);
    }
    std::cout<<' '<<label<<"_tile_count="<<count<<"_input_rows="<<rows
        <<"_active_rows="<<active_rows<<"_active_row_ratio="<<ratio(active_rows,rows)
        <<"_zero_rows_nonzero_rhs="<<zero_rhs<<"_structural_nnz="<<structural<<"_exact_nnz="<<numeric
        <<"_numeric_active_columns="<<active_columns
        <<"_dense_occupancy="<<occupancy
        <<"_peak_tile_temporary_bytes="<<temporary_peak<<"_factor_shapes=";
    bool first=true;
    for(const auto & [shape,frequency]:shapes)
    {
        std::cout<<(first?"":",")<<shape.first<<'x'<<shape.second<<':'<<frequency;
        first=false;
    }
    std::cout<<"_structural_row_ratio_bins=";
    for(std::size_t i=0;i<structural_bins.size();++i) std::cout<<(i?",":"")<<structural_bins[i];
    std::cout<<"_exact_row_ratio_bins=";
    for(std::size_t i=0;i<exact_bins.size();++i) std::cout<<(i?",":"")<<exact_bins[i];
    std::cout<<"_dense_occupancy_bins=";
    for(std::size_t i=0;i<occupancy_bins.size();++i) std::cout<<(i?",":"")<<occupancy_bins[i];
    if(stage==joint::TiledQrStage::DerivativePrepare)
        std::cout<<"_free_design_structural_nnz="<<design_structural
            <<"_free_design_exact_nnz="<<design_numeric<<"_raw_derivative_structural_nnz="<<derivative_structural
            <<"_raw_derivative_exact_nnz="<<derivative_numeric<<"_free_design_nonzero_rows="<<design_active
            <<"_free_design_dense_occupancy="<<ratio(design_numeric,design_cells)
            <<"_raw_derivative_nonzero_rows="<<derivative_active
            <<"_raw_derivative_dense_occupancy="<<ratio(derivative_numeric,derivative_cells)
            <<"_both_zero_rows="<<both_zero;
    if(stage==joint::TiledQrStage::DerivativeReduce)
        std::cout<<"_corrected_jacobian_exact_nnz="<<corrected_numeric
            <<"_corrected_jacobian_nonzero_rows="<<corrected_active;
    if(stage==joint::TiledQrStage::ProjectedWidth)
        std::cout<<"_projection_only_effective_rows="<<projection_active;
    std::cout<<'\n';
}

Sample Run(std::size_t targets,bool telemetry)
{
    joint::SetHotspotAttributionEnabled(telemetry);
    joint::ResetProjectedWidthTileCensus();
    rhbm_gem::core::JointProblem problem(second_stage_test::SyntheticJointWorkload(
        "chain",static_cast<int>(targets)));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const joint::Vector initial_eta=joint::Vector::Constant(
        static_cast<Eigen::Index>(targets),std::log(.55));
    joint::FixedNeighborPolicy policy;
    policy.core_atoms=12;
    policy.convergence_criteria=joint::PracticalJointConvergenceCriteria;
    policy.eta_upper_bound=std::log(2.2);
    policy.collect_telemetry=telemetry;
    policy.assess_local_conditioning=true;
    policy.assess_parameter_stability=true;
    std::vector<std::uint64_t> candidate_decisions;
    policy.candidate_observer=[&](const joint::FixedNeighborCandidateEvent & event) {
        candidate_decisions.push_back(static_cast<std::uint64_t>(event.sweep));
        candidate_decisions.push_back(static_cast<std::uint64_t>(event.block));
        candidate_decisions.push_back(event.accepted);
        candidate_decisions.push_back(event.unchanged);
        candidate_decisions.push_back(event.candidate_replay_ran);
    };

    const auto started=Clock::now();
    const auto result=joint::SearchFixedNeighborComponent(data.domain,data.layout,data.y,
        initial_eta,data.context,policy);
    Sample sample;
    sample.wall_seconds=Seconds(started);
    if(!result.final_sweep || !result.certification)
        throw std::runtime_error("missing completed or certified sweep: "+result.reason);
    const auto & certification=*result.certification;
    const auto & independent=certification.assessment.independent_replay;
    if(!result.search_converged || !certification.qualified() || !independent.passed)
        throw std::runtime_error("qualification failed: "+result.reason);
    sample.replay=result.work.candidate_replay_breakdown;
    sample.sweeps=result.sweep_count;
    sample.accepted_updates=result.total_accepted_local_updates;
    sample.active_width_bounds=result.final_sweep->active_width_bound_count;
    sample.candidate_events=candidate_decisions.size()/5;
    sample.candidate_decision_digest=Digest(candidate_decisions);
    sample.profile=result.work.local_profile_work;
    sample.candidate_copy_seconds=result.work.candidate_copy_seconds;
    sample.candidate_replay_seconds=result.work.candidate_replay_seconds;
    sample.cache_update_seconds=result.work.cache_update_seconds;
    sample.sweep_replay_seconds=result.work.sweep_replay_seconds;
    sample.replay_audit_seconds=result.work.sweep_replay_audit_seconds;
    sample.unchanged_fast_paths=result.work.candidate_replay_skipped_unchanged;
    sample.full_replay_fallbacks=result.work.full_candidate_replays;
    sample.peak_rss_bytes=PeakRssBytes();
    sample.objective=independent.independent_objective/(data.context.scale*data.context.scale);
    sample.ac_kkt=result.final_sweep->global_ac_kkt;
    sample.width_kkt=result.final_sweep->projected_width_gradient_inf_norm;
    sample.eta_digest=Digest(result.state.eta);
    sample.beta_digest=Digest(result.state.beta);
    sample.stop_reason=result.reason;
    sample.qualified=certification.qualified();
    sample.independent_replay_passed=independent.passed;
    sample.projected_width_census=joint::ProjectedWidthTileCensus();
    return sample;
}

void PrintHistogram(const char * label,const std::array<std::size_t,16> & histogram)
{
    std::cout<<' '<<label<<'=';
    for(std::size_t i=0;i<histogram.size();++i) std::cout<<(i ? "," : "")<<histogram[i];
}

void PrintSummary(std::size_t targets,bool telemetry,const std::vector<Sample> & samples)
{
    const auto & representative=samples[samples.size()/2];
    const auto & replay=representative.replay;
    const auto & work=representative.profile;
    const auto & linear=work.total.linear_solve;
    const auto & linear_detail=linear.breakdown;
    const auto & derivative=work.total.derivative_breakdown;
    const auto & lm=work.bounded_lm;
    const auto & lm_detail=lm.breakdown;
    const auto peak_rss=std::max_element(samples.begin(),samples.end(),[](const auto & a,const auto & b) {
        return a.peak_rss_bytes<b.peak_rss_bytes;
    })->peak_rss_bytes;
    std::cout<<std::setprecision(9)<<"targets="<<targets<<" telemetry="<<(telemetry ? "on" : "off")
        <<" warmups=1 measured_runs="<<samples.size();
    PrintTime("search_wall",samples,[](const auto & x){return x.wall_seconds;});
    PrintTime("profile_search",samples,[](const auto & x){return x.profile.total_seconds;});
    PrintTime("candidate_copy",samples,[](const auto & x){return x.candidate_copy_seconds;});
    PrintTime("candidate_replay",samples,[](const auto & x){return x.candidate_replay_seconds;});
    PrintTime("cache_update",samples,[](const auto & x){return x.cache_update_seconds;});
    PrintTime("sweep_full_replay",samples,[](const auto & x){return x.sweep_replay_seconds;});
    PrintTime("replay_audit",samples,[](const auto & x){return x.replay_audit_seconds;});
    PrintTime("ac_row_reduction",samples,[](const auto & x){return x.profile.total.linear_solve.breakdown.row_reduction_seconds;});
    PrintTime("ac_tile_dense_assembly",samples,[](const auto & x){return TileTime(
        x.profile.total.linear_solve.breakdown.tiled_qr_tiles,joint::TiledQrStage::SparseRows,
        &joint::TiledQrTileCensus::dense_tile_materialization_seconds);});
    PrintTime("ac_householder_factorization",samples,[](const auto & x){return TileTime(
        x.profile.total.linear_solve.breakdown.tiled_qr_tiles,joint::TiledQrStage::SparseRows,
        &joint::TiledQrTileCensus::householder_factorization_seconds);});
    PrintTime("ac_qtrhs_transform",samples,[](const auto & x){return TileTime(
        x.profile.total.linear_solve.breakdown.tiled_qr_tiles,joint::TiledQrStage::SparseRows,
        &joint::TiledQrTileCensus::qtrhs_transform_seconds);});
    PrintTime("ac_r_extraction",samples,[](const auto & x){return TileTime(
        x.profile.total.linear_solve.breakdown.tiled_qr_tiles,joint::TiledQrStage::SparseRows,
        &joint::TiledQrTileCensus::r_extraction_seconds);});
    PrintTime("ac_sparse_qr_analysis",samples,[](const auto & x){return x.profile.total.linear_solve.breakdown.sparse_qr_analysis_seconds;});
    PrintTime("ac_sparse_qr_factorization",samples,[](const auto & x){return x.profile.total.linear_solve.breakdown.sparse_qr_factorization_seconds;});
    PrintTime("ac_sparse_qr_solve",samples,[](const auto & x){return x.profile.total.linear_solve.breakdown.sparse_qr_solve_seconds;});
    PrintTime("derivative_prepare",samples,[](const auto & x){return x.profile.total.derivative_prepare_seconds;});
    PrintTime("derivative_prepare_qr",samples,[](const auto & x){return x.profile.total.derivative_breakdown.tiled_qr_seconds;});
    PrintTime("derivative_prepare_materialization",samples,[](const auto & x){return TileTime(
        x.profile.total.derivative_breakdown.tiled_qr_tiles,joint::TiledQrStage::DerivativePrepare,
        &joint::TiledQrTileCensus::dense_tile_materialization_seconds);});
    PrintTime("derivative_prepare_householder",samples,[](const auto & x){return TileTime(
        x.profile.total.derivative_breakdown.tiled_qr_tiles,joint::TiledQrStage::DerivativePrepare,
        &joint::TiledQrTileCensus::householder_factorization_seconds);});
    PrintTime("derivative_prepare_qtrhs",samples,[](const auto & x){return TileTime(
        x.profile.total.derivative_breakdown.tiled_qr_tiles,joint::TiledQrStage::DerivativePrepare,
        &joint::TiledQrTileCensus::qtrhs_transform_seconds);});
    PrintTime("derivative_reduce",samples,[](const auto & x){return x.profile.total.derivative_reduce_seconds;});
    PrintTime("derivative_corrected_jacobian",samples,[](const auto & x){return x.profile.total.derivative_breakdown.corrected_jacobian_seconds;});
    PrintTime("derivative_reduce_householder",samples,[](const auto & x){return TileTime(
        x.profile.total.derivative_breakdown.tiled_qr_tiles,joint::TiledQrStage::DerivativeReduce,
        &joint::TiledQrTileCensus::householder_factorization_seconds);});
    PrintTime("derivative_reduce_qtrhs",samples,[](const auto & x){return TileTime(
        x.profile.total.derivative_breakdown.tiled_qr_tiles,joint::TiledQrStage::DerivativeReduce,
        &joint::TiledQrTileCensus::qtrhs_transform_seconds);});
    std::cout<<" candidate_attempts="<<replay.attempts
        <<" affected_candidates="<<replay.affected_row_candidates
        <<" unchanged_fast_paths="<<representative.unchanged_fast_paths
        <<" full_replay_fallbacks="<<representative.full_replay_fallbacks
        <<" candidate_incidence_prepare_s="<<replay.incidence_preparation_seconds
        <<" candidate_incidence_ordering_s="<<replay.incidence_ordering_seconds
        <<" candidate_lookup_count_s="<<replay.incidence_lookup_count_seconds
        <<" candidate_incidence_construction_s="<<replay.incidence_construction_seconds
        <<" candidate_row_mapping_preparation_s="<<replay.row_mapping_preparation_seconds
        <<" candidate_prediction_residual_s="<<replay.prediction_residual_seconds
        <<" candidate_incremental_objective_s="<<replay.incremental_objective_seconds
        <<" candidate_scratch_peak_bytes="<<replay.peak_scratch_bytes
        <<" affected_rows_total="<<replay.affected_rows_total
        <<" affected_rows_min="<<(replay.affected_rows_min==std::numeric_limits<std::size_t>::max()?0:replay.affected_rows_min)
        <<" affected_rows_max="<<replay.affected_rows_max
        <<" contributor_entries_total="<<replay.contributor_entries_total
        <<" contributor_entries_min="<<(replay.contributor_entries_min==std::numeric_limits<std::size_t>::max()?0:replay.contributor_entries_min)
        <<" contributor_entries_max="<<replay.contributor_entries_max;
    PrintHistogram("affected_rows_log2_bins",replay.affected_row_count_bins);
    PrintHistogram("contributor_entries_log2_bins",replay.contributor_entry_count_bins);
    std::cout<<" profile_evaluations="<<work.total.evaluations
        <<" derivative_preparations="<<work.total.derivative_preparations
        <<" derivative_reductions="<<work.total.derivative_reductions
        <<" ac_solves="<<linear.solve_calls<<" ac_valid_solves="<<linear.valid_solves
        <<" ac_linear_solve_total_s="<<linear.solver_total_seconds
        <<" ac_sparse_qr_inclusive_s="<<linear.sparse_qr_seconds
        <<" ac_initial_projected_kkt_max="<<linear.initial_kkt_max
        <<" ac_final_projected_kkt_max="<<linear.final_kkt_max
        <<" ac_active_set_iterations="<<linear.active_set_iterations
        <<" ac_active_set_releases="<<linear.active_set_releases
        <<" sparse_structure_builds="<<linear.sparse_structure_build_count
        <<" sparse_qr_factorizations="<<linear_detail.sparse_qr_factorizations
        <<" free_face_changes="<<linear.free_face_changes
        <<" free_face_unchanged="<<linear.free_face_unchanged
        <<" free_face_repeat_transitions="<<linear_detail.free_face_repeats_across_builds
        <<" free_face_change_transitions="<<linear_detail.free_face_changes_across_builds
        <<" selected_pattern_repeat_transitions="<<linear_detail.selected_pattern_repeats
        <<" selected_pattern_change_transitions="<<linear_detail.selected_pattern_changes
        <<" reduced_pattern_repeat_transitions="<<linear_detail.reduced_pattern_repeats
        <<" reduced_pattern_change_transitions="<<linear_detail.reduced_pattern_changes
        <<" reduced_numeric_repeat_same_pattern="<<linear_detail.reduced_numeric_repeats_on_same_pattern
        <<" reduced_numeric_change_same_pattern="<<linear_detail.reduced_numeric_changes_on_same_pattern
        <<" selected_rows_range="<<(linear_detail.selected_matrix_builds?linear_detail.selected_rows_min:0)
        <<".."<<linear_detail.selected_rows_max
        <<" selected_columns_range="<<(linear_detail.selected_matrix_builds?linear_detail.selected_columns_min:0)
        <<".."<<linear_detail.selected_columns_max
        <<" selected_nonzeros_range="<<(linear_detail.selected_matrix_builds?linear_detail.selected_nonzeros_min:0)
        <<".."<<linear_detail.selected_nonzeros_max
        <<" reduced_rows_range="<<(linear_detail.selected_matrix_builds?linear_detail.reduced_rows_min:0)
        <<".."<<linear_detail.reduced_rows_max
        <<" reduced_nonzeros_range="<<(linear_detail.selected_matrix_builds?linear_detail.reduced_nonzeros_min:0)
        <<".."<<linear_detail.reduced_nonzeros_max
        <<" selected_nonzeros="<<linear.selected_nonzeros
        <<" reduced_nonzeros="<<linear.reduced_nonzeros
        <<" column_normalization_s="<<linear_detail.column_normalization_seconds
        <<" sparse_structure_build_s="<<linear.sparse_structure_build_seconds
        <<" row_reduction_s="<<linear_detail.row_reduction_seconds
        <<" row_major_conversion_s="<<linear.row_major_conversion_seconds
        <<" tile_mapping_build_s="<<linear.tile_mapping_build_seconds
        <<" tiled_sparse_row_qr_s="<<linear.tiled_qr_seconds
        <<" ac_dense_tile_materialization_s="<<TileTime(linear_detail.tiled_qr_tiles,
            joint::TiledQrStage::SparseRows,&joint::TiledQrTileCensus::dense_tile_materialization_seconds)
        <<" ac_householder_factorization_s="<<TileTime(linear_detail.tiled_qr_tiles,
            joint::TiledQrStage::SparseRows,&joint::TiledQrTileCensus::householder_factorization_seconds)
        <<" ac_qtrhs_transform_s="<<TileTime(linear_detail.tiled_qr_tiles,
            joint::TiledQrStage::SparseRows,&joint::TiledQrTileCensus::qtrhs_transform_seconds)
        <<" ac_r_extraction_s="<<TileTime(linear_detail.tiled_qr_tiles,
            joint::TiledQrStage::SparseRows,&joint::TiledQrTileCensus::r_extraction_seconds)
        <<" reduced_sparse_assembly_s="<<linear.reduced_sparse_assembly_seconds
        <<" sparse_qr_analysis_s="<<linear_detail.sparse_qr_analysis_seconds
        <<" sparse_qr_factorization_s="<<linear_detail.sparse_qr_factorization_seconds
        <<" sparse_qr_solve_s="<<linear_detail.sparse_qr_solve_seconds
        <<" pattern_fingerprint_s="<<linear_detail.pattern_fingerprint_seconds
        <<" reduced_pattern_fingerprint_s="<<linear_detail.reduced_pattern_fingerprint_seconds
        <<" refinement_eligible="<<linear.refinement_eligibility_count
        <<" refinement_triggered="<<linear.refinement_trigger_count
        <<" refinement_succeeded="<<linear.refinement_success_count
        <<" refinement_failed="<<linear.refinement_failure_count
        <<" refinement_iterations="<<linear.refinement_iterations
        <<" refinement_factorizations="<<linear.refinement_factorizations
        <<" refinement_s="<<linear.refinement_seconds
        <<" refinement_matrix_s="<<linear_detail.refinement_matrix_seconds
        <<" refinement_qr_s="<<linear_detail.refinement_qr_factorization_seconds
        <<" refinement_correction_solve_s="<<linear_detail.refinement_correction_solve_seconds
        <<" refinement_certification_s="<<linear_detail.refinement_certification_seconds
        <<" derivative_prepare_s="<<work.total.derivative_prepare_seconds
        <<" derivative_reduce_s="<<work.total.derivative_reduce_seconds
        <<" derivative_raw_construction_s="<<derivative.raw_construction_seconds
        <<" derivative_tiled_qr_s="<<derivative.tiled_qr_seconds
        <<" derivative_tiled_qr_calls="<<derivative.tiled_qr_calls
        <<" derivative_corrected_jacobian_s="<<derivative.corrected_jacobian_seconds
        <<" derivative_prepare_dense_tile_materialization_s="<<TileTime(derivative.tiled_qr_tiles,
            joint::TiledQrStage::DerivativePrepare,&joint::TiledQrTileCensus::dense_tile_materialization_seconds)
        <<" derivative_prepare_qr_assembly_s="<<TileTime(derivative.tiled_qr_tiles,
            joint::TiledQrStage::DerivativePrepare,&joint::TiledQrTileCensus::dense_assembly_seconds)
        <<" derivative_prepare_householder_s="<<TileTime(derivative.tiled_qr_tiles,
            joint::TiledQrStage::DerivativePrepare,&joint::TiledQrTileCensus::householder_factorization_seconds)
        <<" derivative_prepare_qtrhs_s="<<TileTime(derivative.tiled_qr_tiles,
            joint::TiledQrStage::DerivativePrepare,&joint::TiledQrTileCensus::qtrhs_transform_seconds)
        <<" derivative_prepare_r_extraction_s="<<TileTime(derivative.tiled_qr_tiles,
            joint::TiledQrStage::DerivativePrepare,&joint::TiledQrTileCensus::r_extraction_seconds)
        <<" derivative_reduce_qr_assembly_s="<<TileTime(derivative.tiled_qr_tiles,
            joint::TiledQrStage::DerivativeReduce,&joint::TiledQrTileCensus::dense_assembly_seconds)
        <<" derivative_reduce_householder_s="<<TileTime(derivative.tiled_qr_tiles,
            joint::TiledQrStage::DerivativeReduce,&joint::TiledQrTileCensus::householder_factorization_seconds)
        <<" derivative_reduce_qtrhs_s="<<TileTime(derivative.tiled_qr_tiles,
            joint::TiledQrStage::DerivativeReduce,&joint::TiledQrTileCensus::qtrhs_transform_seconds)
        <<" derivative_reduce_r_extraction_s="<<TileTime(derivative.tiled_qr_tiles,
            joint::TiledQrStage::DerivativeReduce,&joint::TiledQrTileCensus::r_extraction_seconds)
        <<" derivative_rank_svd_s="<<derivative.rank_svd_seconds
        <<" derivative_rank_svd_calls="<<derivative.rank_svd_calls
        <<" derivative_triangular_solves_s="<<derivative.triangular_solve_seconds
        <<" derivative_triangular_solve_calls="<<derivative.triangular_solve_calls
        <<" derivative_correction_matrix_s="<<derivative.correction_matrix_seconds
        <<" derivative_cancellation_fallbacks="<<derivative.cancellation_fallbacks
        <<" derivative_cancellation_fallback_s="<<derivative.cancellation_fallback_seconds
        <<" derivative_workspace_resize_events="<<derivative.workspace_resize_events
        <<" lm_subproblems="<<lm.subproblem_solves<<" lm_failed_subproblems="<<lm.failed_subproblems
        <<" lm_zero_damping_solves="<<lm.zero_damping_solves
        <<" lm_damping_solves="<<lm.damping_solves<<" lm_bracket_solves="<<lm.bracket_solves
        <<" lm_bisection_solves="<<lm.bisection_solves
        <<" lm_qr_factorizations="<<lm.damping_qr_factorizations
        <<" lm_active_set_iterations="<<lm.damping_active_set_iterations
        <<" lm_releases="<<lm.damping_releases
        <<" lm_active_selection_s="<<lm_detail.active_free_selection_seconds
        <<" lm_augmented_assembly_s="<<lm_detail.augmented_matrix_assembly_seconds
        <<" lm_dense_qr_s="<<lm_detail.dense_qr_factorization_seconds
        <<" lm_qr_solve_s="<<lm_detail.qr_solve_seconds
        <<" lm_kkt_validation_s="<<lm_detail.kkt_validation_seconds
        <<" lm_damping_search_s="<<lm_detail.damping_search_seconds
        <<" lm_updates_per_log2_bin=";
    for(std::size_t i=0;i<lm.damping_solves_per_update_bins.size();++i)
        std::cout<<(i ? "," : "")<<lm.damping_solves_per_update_bins[i];
    std::cout<<" lm_qr_per_update_log2_bin=";
    for(std::size_t i=0;i<lm.qr_factorizations_per_update_bins.size();++i)
        std::cout<<(i ? "," : "")<<lm.qr_factorizations_per_update_bins[i];
    std::cout<<" peak_rss_bytes="<<peak_rss
        <<" objective="<<representative.objective<<" ac_kkt="<<representative.ac_kkt
        <<" width_kkt="<<representative.width_kkt<<" stop_reason="<<representative.stop_reason
        <<" sweeps="<<representative.sweeps<<" accepted_local_updates="<<representative.accepted_updates
        <<" active_width_bounds="<<representative.active_width_bounds
        <<" candidate_events="<<representative.candidate_events
        <<" candidate_decisions_sha256="<<representative.candidate_decision_digest
        <<" eta_sha256="<<representative.eta_digest<<" beta_sha256="<<representative.beta_digest
        <<" qualified="<<representative.qualified
        <<" independent_replay="<<representative.independent_replay_passed<<'\n';
    PrintTileCensus("ac_sparse_rows",linear_detail.tiled_qr_tiles,joint::TiledQrStage::SparseRows);
    PrintTileCensus("derivative_prepare",derivative.tiled_qr_tiles,joint::TiledQrStage::DerivativePrepare);
    PrintTileCensus("derivative_reduce",derivative.tiled_qr_tiles,joint::TiledQrStage::DerivativeReduce);
    PrintTileCensus("projected_width",representative.projected_width_census,joint::TiledQrStage::ProjectedWidth);
}

void PrintCancellationSummary()
{
    std::vector<std::vector<joint::Support>> support(1);
    support[0]={{0,0.0},{1,0.25},{2,0.25}};
    joint::Vector y(3);
    for(Eigen::Index row=0;row<3;++row)
    {
        const double square=row==0 ? 0.0 : 0.25;
        const auto basis=joint::EvaluateKernel(square,.55,2.5);
        y(row)=3.0*basis.gaussian+.25*basis.charge+(row==1 ? 1e-4 : row==2 ? -1e-4 : 0.0);
    }
    const joint::Domain domain(3,std::move(support));
    const auto context=second_stage_test::matched::joint_abc::MakeContext(y,1);
    const auto eta=joint::Vector::Constant(1,std::log(.55));
    const auto evaluation=joint::EvaluateProfile(domain,y,eta,&context);
    if(!evaluation.valid) throw std::runtime_error("cancellation fixture profile failed: "+evaluation.reason);
    std::vector<double> total_seconds,fallback_seconds;
    std::size_t fallback_count{};
    for(std::size_t run=0;run<6;++run)
    {
        const auto started=Clock::now();
        const auto differential=joint::PrepareDerivative(evaluation,context.scale,&context,-1,3);
        const auto elapsed=Seconds(started);
        if(!differential.valid || !differential.reference_order || differential.breakdown.cancellation_fallbacks!=1)
            throw std::runtime_error("cancellation fixture did not use one valid reference-order fallback");
        fallback_count=differential.breakdown.cancellation_fallbacks;
        if(run)
        {
            total_seconds.push_back(elapsed);
            fallback_seconds.push_back(differential.breakdown.cancellation_fallback_seconds);
        }
    }
    const auto [minimum,maximum]=std::minmax_element(fallback_seconds.begin(),fallback_seconds.end());
    std::cout<<std::setprecision(9)<<"cancellation_fixture warmups=1 measured_runs=5 fallbacks="<<fallback_count
        <<" total_median_s="<<Median(total_seconds)<<" fallback_extra_median_s="<<Median(fallback_seconds)
        <<" fallback_extra_range_s=["<<*minimum<<','<<*maximum<<"] reference_order=1\n";
}

void PrintControlledLinearCensus()
{
    struct ControlledCase
    {
        const char * name;
        Eigen::Index columns;
        double active_row_ratio,support_ratio;
        bool active_a{},rank_deficient{},near_rank_deficient{};
    };
    std::vector<ControlledCase> cases;
    for(const auto columns:{Eigen::Index{4},Eigen::Index{24}})
        for(const double active_ratio:{.25,.5,.75,1.0})
            for(const double support_ratio:{.25,1.0})
                cases.push_back({columns==4 ? "narrow" : "wide",columns,active_ratio,support_ratio});
    cases.push_back({"active-a-signed-c",24,1.0,.25,true});
    cases.push_back({"rank-deficient",24,1.0,1.0,false,true});
    cases.push_back({"near-rank-deficient",24,1.0,1.0,false,false,true});

    constexpr Eigen::Index rows=2053;
    constexpr double pi=std::numbers::pi;
    for(const auto & test:cases)
    {
        const Eigen::Index active_rows=static_cast<Eigen::Index>(std::llround(
            test.active_row_ratio*static_cast<double>(rows)));
        joint::Matrix dense=joint::Matrix::Zero(rows,test.columns);
        for(Eigen::Index row=0;row<active_rows;++row)
        {
            const Eigen::Index supported=std::max<Eigen::Index>(1,
                static_cast<Eigen::Index>(std::llround(test.support_ratio*static_cast<double>(test.columns))));
            for(Eigen::Index col=0;col<test.columns;++col)
            {
                const bool present=test.support_ratio==1.0 ||
                    ((col-row%test.columns+test.columns)%test.columns)<supported;
                if(!present) continue;
                const double phase=2*pi*static_cast<double>((col+1)*(row+1))/static_cast<double>(active_rows);
                dense(row,col)=std::sin(phase)+.37*std::cos(phase)+.001*static_cast<double>(row+col);
            }
        }
        if(test.rank_deficient) dense.col(test.columns-1)=dense.col(test.columns-3);
        if(test.near_rank_deficient) dense.col(test.columns-1)=dense.col(test.columns-3)+
            2e-12*dense.col(test.columns-2);
        joint::Vector truth(test.columns);
        for(Eigen::Index col=0;col<test.columns;++col)
            truth(col)=col%2 ? -.12-.003*static_cast<double>(col) : .4+.02*static_cast<double>(col);
        if(test.active_a) truth(0)=-.8;
        joint::Vector rhs=dense*truth;
        for(Eigen::Index row=0;row<rows;++row)
            rhs(row)+=row<active_rows ? .001*std::sin(.13*static_cast<double>(row)) :
                .05*(row%2 ? 1.0 : -1.0);
        const joint::Sparse sparse=dense.sparseView(0.0,0.0);
        joint::SetHotspotAttributionEnabled(true);
        const auto solve=[&]() {
            joint::LinearSolveWorkspace workspace;
            const auto started=Clock::now();
            auto result=joint::SolveLinear(sparse,rhs,nullptr,0,&workspace);
            return std::pair{std::move(result),Seconds(started)};
        };
        (void)solve();
        std::vector<double> times;
        joint::LinearResult representative;
        for(std::size_t i=0;i<5;++i)
        {
            auto [result,elapsed]=solve();
            times.push_back(elapsed);
            if(i==2) representative=std::move(result);
        }
        if((test.rank_deficient && (representative.valid || representative.reason!="rank-deficient")) ||
            (!test.rank_deficient && !test.near_rank_deficient && !representative.valid))
            throw std::runtime_error(std::string("controlled QR case failed: ")+test.name);
        double reference_error=std::numeric_limits<double>::quiet_NaN();
        joint::Certificate certificate;
        if(representative.valid)
        {
            certificate=joint::CertifyLinear(sparse,rhs,representative.beta);
            if(!test.active_a && !test.near_rank_deficient)
            {
                const Eigen::ColPivHouseholderQR<joint::Matrix> reference(dense);
                const joint::Vector dense_beta=reference.solve(rhs);
                reference_error=(representative.beta-dense_beta).norm()/(1+dense_beta.norm());
                if(reference_error>2e-8) throw std::runtime_error("controlled dense QR reference mismatch");
            }
            if(!certificate.available || !certificate.kkt_passed)
                throw std::runtime_error("controlled QR KKT certification failed");
        }
        const auto [minimum,maximum]=std::minmax_element(times.begin(),times.end());
        std::cout<<std::setprecision(9)<<"controlled_case="<<test.name<<" columns="<<test.columns
            <<" rows="<<rows<<" nonzero_design_row_ratio="<<static_cast<double>(active_rows)/rows
            <<" requested_support_ratio="<<test.support_ratio<<" warmups=1 measured_runs=5"
            <<" solve_median_s="<<Median(times)<<" solve_range_s=["<<*minimum<<','<<*maximum<<']'
            <<" valid="<<representative.valid<<" reason="<<representative.reason<<" rank="<<representative.rank
            <<" active_set_iterations="<<representative.active_set_iterations
            <<" active_set_releases="<<representative.releases
            <<" objective="<<certificate.objective<<" projected_kkt="<<certificate.projected_kkt
            <<" dense_qr_relative_error="<<reference_error
            <<" zero_design_nonzero_rhs_rows="<<(rows-active_rows)
            <<" peak_rss_bytes="<<PeakRssBytes()<<'\n';
        PrintTileCensus("controlled_sparse_rows",representative.breakdown.tiled_qr_tiles,
            joint::TiledQrStage::SparseRows);
    }
}
}

int main(int argc,char ** argv)
{
    if(argc!=3)
    {
        std::cerr<<"usage: rhbm_joint_hotspot_attribution_benchmark <target-count> <on|off>\n";
        return 2;
    }
    try
    {
        const auto targets=static_cast<std::size_t>(std::stoul(argv[1]));
        if(targets==0 || targets>1000000) throw std::invalid_argument("target count is out of range");
        const std::string mode=argv[2];
        if(mode!="on" && mode!="off") throw std::invalid_argument("telemetry mode must be on or off");
        const bool telemetry=mode=="on";
        (void)Run(targets,telemetry);
        std::vector<Sample> samples;
        samples.reserve(5);
        for(std::size_t repetition=0;repetition<5;++repetition)
            samples.push_back(Run(targets,telemetry));
        PrintSummary(targets,telemetry,samples);
        if(telemetry)
        {
            PrintCancellationSummary();
            if(targets==300) PrintControlledLinearCensus();
        }
        return 0;
    }
    catch(const std::exception & error)
    {
        std::cerr<<"joint hotspot attribution benchmark failed: "<<error.what()<<'\n';
        return 1;
    }
}
