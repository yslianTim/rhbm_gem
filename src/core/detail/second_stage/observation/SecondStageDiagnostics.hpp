#pragma once

#include "core/detail/second_stage/CandidateEvidence.hpp"
#include "core/detail/second_stage/ConvergenceCertificate.hpp"

#include <array>
#include <string_view>

namespace rhbm_gem::core::detail {

constexpr std::size_t kAuditDetailLimit{ 5 };
enum class AuditStage { Proposal, LocalSearch, LocalPolish, BoundaryEndpoint, BoundaryCorrection,
    BoundaryBacktracking, RescueEndpoint, RescueCorrection, RescueBacktracking,
    SelectionOrdinary, SelectionRescue, Commit, Quarantine, Partition, FinalPolish, FinalCertification, Count };
enum class AuditCategory { None, Rejected, Unavailable, Invalid, Guard, Trust, Solver,
    Exhausted, Shrink, Enter, Retry, Release, Salvage, Rescue, Partition, Count };

// Fixed-size evidence only: no model, sample, history, or candidate ownership.
struct AuditEvent
{
    AuditStage stage{ AuditStage::LocalSearch };
    AuditCategory category{ AuditCategory::None };
    std::size_t first_atom{ 0 }, atom_count{ 0 }, trial{ 0 };
    std::string_view outcome{ "accepted" }, reason{}, scope{ "cluster" }, reference{ "iteration_previous" };
    std::optional<ObjectiveBreakdown> previous{}, candidate{}, best{};
    std::optional<double> factor{}, radius{};
    bool previous_checked{ false }, best_checked{ false };
};

struct AuditStageCount { std::size_t total{ 0 }, accepted{ 0 }, rejected{ 0 }, skipped{ 0 }; };
struct AuditBatch
{
    std::array<AuditStageCount, static_cast<std::size_t>(AuditStage::Count)> stages{};
    std::array<std::size_t, static_cast<std::size_t>(AuditCategory::Count)> categories{};
    std::array<AuditEvent, kAuditDetailLimit> details{};
    std::size_t detail_count{ 0 }, abnormal_count{ 0 };
    void Add(const AuditEvent &) noexcept;
    void Merge(const AuditBatch &) noexcept;
};

struct SelectionAuditDiagnostic
{
    bool executed{ false };
    std::string_view result{ "skipped" }, reason{ "no-ordinary-components" };
    std::size_t evaluations{ 0 }, removed_clusters{ 0 };
    std::optional<ObjectiveBreakdown> previous{}, candidate{}, best{};
};

struct IterationDiagnostics
{
    std::optional<double> accepted_maximum_transformed_change{};
    double proposal_maximum_transformed_change{ 0.0 };
};
struct IterationObservation { IterationDiagnostics diagnostics{}; };

// Basic final progress is independent of the optional audit payload.
struct FinalDependencyPolishDiagnostic
{
    std::size_t component_count{ 0 }, attempted_component_count{ 0 }, accepted_component_count{ 0 };
    std::size_t atom_count{ 0 }, parameter_count{ 0 }, round_count{ 0 }, suspicious_candidate_atom_count{ 0 };
    std::optional<double> objective_before{}, objective_after{};
    double elapsed_milliseconds{ 0.0 };
};

struct SecondStageSeedSummary { std::size_t local_mdpde{ 0 }, global_median{ 0 }; };

} // namespace rhbm_gem::core::detail
