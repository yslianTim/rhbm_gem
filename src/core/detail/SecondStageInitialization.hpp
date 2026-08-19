#pragma once

#include <array>
#include <cstddef>
#include <iomanip>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

#include "core/detail/FitStateView.hpp"
#include "core/detail/SecondStageContext.hpp"
#include "core/detail/TransformedGaussianModel.hpp"
#include "data/detail/AtomClassifier.hpp"

namespace rhbm_gem::core::detail {

enum class SecondStageSeedSource
{
    GroupPosterior,
    GroupPrior,
    GroupMedian,
    GlobalMedian
};

struct SecondStageSeedCandidates
{
    std::optional<GaussianModel3DWithUncertainty> group_posterior{};
    std::optional<GaussianModel3DWithUncertainty> group_prior{};
    std::optional<GaussianModel3DWithUncertainty> group_median{};
    std::optional<GaussianModel3DWithUncertainty> global_median{};
};

struct SecondStageSeedSelection
{
    SecondStageSeedSource source{ SecondStageSeedSource::GlobalMedian };
    GaussianModel3DWithUncertainty model{};
};

inline std::optional<SecondStageSeedSelection> SelectSecondStageSeed(
    const SecondStageSeedCandidates & candidates)
{
    const auto select = [](
        SecondStageSeedSource source,
        const std::optional<GaussianModel3DWithUncertainty> & candidate)
        -> std::optional<SecondStageSeedSelection>
    {
        if (!candidate.has_value() ||
            !IsValidSecondStageGaussianModel(candidate->GetModel()))
        {
            return std::nullopt;
        }
        return SecondStageSeedSelection{ source, *candidate };
    };

    if (const auto selected{
            select(SecondStageSeedSource::GroupPosterior, candidates.group_posterior) })
    {
        return selected;
    }
    if (const auto selected{
            select(SecondStageSeedSource::GroupPrior, candidates.group_prior) })
    {
        return selected;
    }
    if (const auto selected{
            select(SecondStageSeedSource::GroupMedian, candidates.group_median) })
    {
        return selected;
    }
    return select(SecondStageSeedSource::GlobalMedian, candidates.global_median);
}

constexpr double kNeighborContributionDistanceMax{ 2.5 };
constexpr double kNeighborAtomSearchRange{ 2.0 * kNeighborContributionDistanceMax };

struct SecondStageSeedSelectionRecord
{
    std::size_t atom_index{ 0 };
    detail::SecondStageSeedSource source{ detail::SecondStageSeedSource::GlobalMedian };
    GaussianModel3D original_model{};
    GaussianModel3D selected_model{};
};

struct UnselectedSecondStageSeedSelectionRecord
{
    int atom_serial_id{ 0 };
    detail::SecondStageSeedSource source{ detail::SecondStageSeedSource::GlobalMedian };
    GaussianModel3D selected_model{};
};

struct SecondStageInitialStateBuildResult
{
    FitState state{};
    std::vector<SecondStageSeedSelectionRecord> selection_record_list{};
    std::vector<UnselectedSecondStageSeedSelectionRecord> unselected_selection_record_list{};
    enum class Failure
    {
        None,
        SelectedSeedUnavailable,
        UnselectedSeedUnavailable
    } failure{ Failure::None };
};

inline SecondStageContext BuildSecondStageContext(
    const ModelObject & model_object,
    const FitOptions & options)
{
    SecondStageContext context;
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    context.selected_atom_list.reserve(atom_list.size());
    for (auto * atom : atom_list)
    {
        context.selected_atom_list.emplace_back(AtomContext{ atom });
    }
    std::unordered_map<GroupKey, std::size_t> selected_group_id_by_key;
    selected_group_id_by_key.reserve(context.size());
    const auto analysis_view{ model_object.GetAnalysisView() };
    std::unordered_map<const AtomObject *, std::size_t> atom_index_map;
    atom_index_map.reserve(context.size());
    for (std::size_t i = 0; i < context.size(); i++)
    {
        atom_index_map.emplace(context.at(i).atom, i);
    }
    std::unordered_map<const AtomObject *, std::size_t> unselected_atom_index_map;
    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        auto & atom_context{ context.at(atom_index) };
        const auto * atom{ atom_context.atom };
        atom_context.group_key = data_internal::GetGroupKey(atom);
        atom_context.residue_key = {
            atom->GetChainID(),
            atom->GetSequenceID()
        };
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
        atom_context.raw_sampling_entries = local_view.GetRawSamplingEntries(false);
        atom_context.initial_result = local_view.GetGaussianResult(FittingStage::Second);
        atom_context.group_prior = analysis_view.FindAtomGroupPriorWithUncertainty(FittingStage::Second, *atom);
        atom_context.alpha_r = local_view.GetAlphaR(FittingStage::Second);
        atom_context.refit_design_template = BuildLocalGaussianDesignTemplate(
            atom_context.raw_sampling_entries,
            options.distance_min,
            options.distance_max);
        auto [group_iter, inserted]{
            selected_group_id_by_key.emplace(
                atom_context.group_key,
                context.selected_atom_index_list_by_group.size())
        };
        if (inserted)
        {
            context.selected_atom_index_list_by_group.emplace_back();
        }
        atom_context.group_id = group_iter->second;
        context.selected_atom_index_list_by_group.at(atom_context.group_id).emplace_back(atom_index);
    }

    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        auto & atom_context{ context.at(atom_index) };
        const auto * atom{ atom_context.atom };
        const auto neighbor_atom_list{ atom->FindNeighborAtoms(kNeighborAtomSearchRange) };
        std::unordered_set<const AtomObject *> neighbor_atom_set;

        atom_context.neighbor_atom_sample_offset_list.reserve(atom_context.raw_sampling_entries.size() + 1);
        atom_context.neighbor_atom_sample_offset_list.emplace_back(0);
        atom_context.neighbor_atom_sample_list.reserve(
            atom_context.raw_sampling_entries.size() * neighbor_atom_list.size());
        for (std::size_t sample_index = 0; sample_index < atom_context.raw_sampling_entries.size(); sample_index++)
        {
            const auto & sample{ atom_context.raw_sampling_entries.at(sample_index) };
            for (auto * neighbor_atom : neighbor_atom_list)
            {
                if (options.exclude_hydrogen && neighbor_atom->GetElement() == Element::HYDROGEN)
                {
                    continue;
                }
                const auto distance{
                    static_cast<double>(
                        array_helper::ComputeNorm<float>(sample.point.position, neighbor_atom->GetPositionRef()))
                };
                if (distance > kNeighborContributionDistanceMax) continue;
                neighbor_atom_set.emplace(neighbor_atom);

                const auto selected_iter{ atom_index_map.find(neighbor_atom) };
                if (selected_iter != atom_index_map.end())
                {
                    atom_context.neighbor_atom_sample_list.emplace_back(
                        NeighborAtomSample{
                            true,
                            selected_iter->second,
                            distance
                        });
                    continue;
                }

                auto unselected_atom_contributor_iter{
                    unselected_atom_index_map.find(neighbor_atom)
                };
                if (unselected_atom_contributor_iter == unselected_atom_index_map.end())
                {
                    const auto unselected_atom_contributor_index{
                        context.unselected_atom_list.size()
                    };
                    const auto group_key{
                        data_internal::GetGroupKey(neighbor_atom)
                    };
                    const auto selected_group_iter{
                        selected_group_id_by_key.find(group_key)
                    };
                    context.unselected_atom_list.emplace_back(
                        UnselectedAtomContributor{
                            neighbor_atom,
                            group_key,
                            selected_group_iter ==
                                selected_group_id_by_key.end() ?
                                std::nullopt :
                                std::optional<std::size_t>{ selected_group_iter->second }
                        });
                    unselected_atom_contributor_iter = unselected_atom_index_map.emplace(
                        neighbor_atom,
                        unselected_atom_contributor_index).first;
                }
                atom_context.neighbor_atom_sample_list.emplace_back(
                    NeighborAtomSample{
                        false,
                        unselected_atom_contributor_iter->second,
                        distance
                    });
            }
            atom_context.neighbor_atom_sample_offset_list.emplace_back(
                atom_context.neighbor_atom_sample_list.size());
        }
        atom_context.neighbor_count_for_peeling = static_cast<int>(neighbor_atom_set.size());
    }

    return context;
}

inline void StoreSecondStageNeighborCounts(
    ModelObject & model_object,
    const SecondStageContext & context)
{
    auto analysis{ model_object.EditAnalysis() };
    for (const auto & atom_context : context)
    {
        analysis.SetAtomLocalNeighborCountForPeeling(
            *atom_context.atom,
            atom_context.neighbor_count_for_peeling);
    }
}


inline std::optional<GaussianModel3DWithUncertainty> BuildValidGaussianParameterMedian(
    const std::vector<GaussianModel3D> & model_list)
{
    const auto median_model{ BuildGaussianParameterMedian(model_list) };
    if (!median_model.has_value()) return std::nullopt;
    return GaussianModel3DWithUncertainty{
        *median_model,
        GaussianModel3DUncertainty{}
    };
}

inline SecondStageInitialStateBuildResult BuildInitialFitState(
    SecondStageContext & context)
{
    SecondStageInitialStateBuildResult build_result;
    auto & state{ build_result.state };
    state.resize(context.size());
    std::unordered_map<GroupKey, std::vector<GaussianModel3D>> models_by_group;
    std::vector<GaussianModel3D> global_models;
    global_models.reserve(context.size());

    for (std::size_t i = 0; i < context.size(); i++)
    {
        const auto & atom_context{ context.at(i) };
        state.at(i) = atom_context.initial_result;
        const auto group_key{ atom_context.group_key };

        const auto & result{ state.at(i) };
        const auto direct_selection{
            SelectSecondStageSeed(
                SecondStageSeedCandidates{
                    result.posterior,
                    atom_context.group_prior,
                    std::nullopt,
                    std::nullopt
                })
        };
        if (!direct_selection.has_value()) continue;

        models_by_group[group_key].emplace_back(direct_selection->model.GetModel());
        global_models.emplace_back(direct_selection->model.GetModel());
    }

    std::unordered_map<GroupKey, GaussianModel3DWithUncertainty> median_by_group;
    median_by_group.reserve(models_by_group.size());
    for (const auto & [group_key, models] : models_by_group)
    {
        const auto median_model{ BuildValidGaussianParameterMedian(models) };
        if (median_model.has_value())
        {
            median_by_group.emplace(group_key, *median_model);
        }
    }
    const auto global_median{ BuildValidGaussianParameterMedian(global_models) };

    for (std::size_t i = 0; i < context.size(); i++)
    {
        auto & result{ state.at(i) };
        const auto original_model{ result.mdpde.GetModel() };
        const auto & atom_context{ context.at(i) };
        const auto group_key{ atom_context.group_key };
        std::optional<GaussianModel3DWithUncertainty> group_median;
        const auto group_median_iter{ median_by_group.find(group_key) };
        if (group_median_iter != median_by_group.end())
        {
            group_median = group_median_iter->second;
        }
        const auto selection{
            SelectSecondStageSeed(
                SecondStageSeedCandidates{
                    result.posterior,
                    atom_context.group_prior,
                    group_median,
                    global_median
                })
        };
        if (!selection.has_value())
        {
            build_result.failure =
                SecondStageInitialStateBuildResult::Failure::SelectedSeedUnavailable;
            return build_result;
        }

        result.mdpde = selection->model;
        build_result.selection_record_list.emplace_back(
            SecondStageSeedSelectionRecord{
                i,
                selection->source,
                original_model,
                selection->model.GetModel()
            });
    }

    for (std::size_t i = 0; i < context.unselected_atom_list.size(); i++)
    {
        auto & unselected_atom_contributor{
            context.unselected_atom_list.at(i)
        };
        std::optional<GaussianModel3DWithUncertainty> group_median;
        const auto group_median_iter{
            median_by_group.find(unselected_atom_contributor.group_key)
        };
        if (group_median_iter != median_by_group.end())
        {
            group_median = group_median_iter->second;
        }
        const auto selection{
            SelectSecondStageSeed(
                SecondStageSeedCandidates{
                    std::nullopt,
                    std::nullopt,
                    group_median,
                    global_median
                })
        };
        if (!selection.has_value())
        {
            build_result.failure = SecondStageInitialStateBuildResult::Failure::UnselectedSeedUnavailable;
            return build_result;
        }

        unselected_atom_contributor.initial_seed = selection->model;
        build_result.unselected_selection_record_list.emplace_back(
            UnselectedSecondStageSeedSelectionRecord{
                unselected_atom_contributor.atom->GetSerialID(),
                selection->source,
                selection->model.GetModel()
            });
    }
    return build_result;
}

inline const char * GetSecondStageSeedSourceText(SecondStageSeedSource source)
{
    switch (source)
    {
    case SecondStageSeedSource::GroupPosterior:
        return "group-posterior";
    case SecondStageSeedSource::GroupPrior:
        return "group-prior";
    case SecondStageSeedSource::GroupMedian:
        return "group-median";
    case SecondStageSeedSource::GlobalMedian:
        return "global-median";
    }
    throw std::logic_error("Unknown second-stage seed source.");
}

inline void LogSecondStageSeedSelections(
    const std::vector<SecondStageSeedSelectionRecord> & selection_record_list,
    bool quiet_mode)
{
    if (quiet_mode || selection_record_list.empty()) return;

    constexpr std::array<SecondStageSeedSource, 4> source_list{
        SecondStageSeedSource::GroupPosterior,
        SecondStageSeedSource::GroupPrior,
        SecondStageSeedSource::GroupMedian,
        SecondStageSeedSource::GlobalMedian
    };
    std::array<std::size_t, source_list.size()> source_count{};
    for (const auto & record : selection_record_list)
    {
        source_count.at(static_cast<std::size_t>(record.source))++;
    }

    std::ostringstream summary;
    summary << "Selected second-stage initial seeds = "
        << selection_record_list.size() << ", sources = ";
    for (std::size_t i = 0; i < source_list.size(); i++)
    {
        if (i != 0) summary << ", ";
        summary << GetSecondStageSeedSourceText(source_list.at(i))
            << ":" << source_count.at(i);
    }
    summary << ".";
    Logger::Log(LogLevel::Info, summary.str());

    for (const auto & record : selection_record_list)
    {
        std::ostringstream detail_message;
        detail_message << "Second-stage seed selection: atom index = "
            << record.atom_index
            << ", source = " << GetSecondStageSeedSourceText(record.source)
            << std::scientific << std::setprecision(2)
            << ", original MDPDE A/B/C = "
            << record.original_model.GetAmplitude() << "/"
            << record.original_model.GetWidth() << "/"
            << record.original_model.GetOffset()
            << ", selected A/B/C = "
            << record.selected_model.GetAmplitude() << "/"
            << record.selected_model.GetWidth() << "/"
            << record.selected_model.GetOffset() << ".";
        Logger::Log(LogLevel::Debug, detail_message.str());
    }
}

inline void LogUnselectedSecondStageSeedSelections(
    const std::vector<UnselectedSecondStageSeedSelectionRecord> & selection_record_list,
    bool quiet_mode)
{
    if (quiet_mode || selection_record_list.empty()) return;

    std::size_t group_median_count{ 0 };
    std::size_t global_median_count{ 0 };
    for (const auto & record : selection_record_list)
    {
        if (record.source == SecondStageSeedSource::GroupMedian)
        {
            group_median_count++;
        }
        else if (record.source == SecondStageSeedSource::GlobalMedian)
        {
            global_median_count++;
        }
    }

    std::ostringstream summary;
    summary << "Unselected second-stage neighbor seeds = " << selection_record_list.size()
        << ", sources = group-median:" << group_median_count
        << ", global-median:" << global_median_count << ".";
    Logger::Log(LogLevel::Info, summary.str());

    for (const auto & record : selection_record_list)
    {
        std::ostringstream detail_message;
        detail_message
            << "Unselected second-stage neighbor seed selection: serial ID = "
            << record.atom_serial_id
            << ", source = " << GetSecondStageSeedSourceText(record.source)
            << std::scientific << std::setprecision(2)
            << ", seed A/B/C = "
            << record.selected_model.GetAmplitude() << "/"
            << record.selected_model.GetWidth() << "/"
            << record.selected_model.GetOffset() << ".";
        Logger::Log(LogLevel::Debug, detail_message.str());
    }
}

} // namespace rhbm_gem::core::detail
