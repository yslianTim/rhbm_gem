#pragma once

#include <cmath>
#include <cstddef>
#include <optional>
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
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

#include "core/detail/LocalFittingGroupMedian.hpp"
#include "core/detail/LocalFittingSeedRepair.hpp"
#include "core/detail/SecondStageLocalFittingContext.hpp"
#include "data/detail/AtomClassifier.hpp"

namespace rhbm_gem::core::detail {

constexpr double kSecondStageNeighborContributionDistanceMax{ 2.5 };
constexpr double kSecondStageNeighborAtomSearchRange{
    2.0 * kSecondStageNeighborContributionDistanceMax
};

struct SecondStageSeedSelectionRecord
{
    std::size_t atom_index{ 0 };
    detail::SecondStageSeedSource source{
        detail::SecondStageSeedSource::GlobalMedian
    };
    GaussianModel3D original_model{};
    GaussianModel3D selected_model{};
};

struct SecondStageInitialStateBuildResult
{
    LocalFittingState state{};
    std::vector<SecondStageSeedSelectionRecord> selection_record_list{};
    std::vector<SecondStageSeedSelectionRecord> unselected_selection_record_list{};
};


inline SecondStageLocalFittingContext BuildSecondStageLocalFittingContext(
    ModelObject & model_object,
    const FitOptions & options)
{
    SecondStageLocalFittingContext context;
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    context.selected_atom_list.reserve(atom_list.size());
    for (auto * atom : atom_list)
    {
        context.selected_atom_list.emplace_back(SecondStageAtomContext{ atom });
    }
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
        atom_context.group_prior = analysis_view.FindAtomGroupPriorWithUncertainty(
            FittingStage::Second,
            *atom);
        atom_context.alpha_r = local_view.GetAlphaR(FittingStage::Second);
        atom_context.refit_design_template = BuildLocalGaussianDesignTemplate(
            atom_context.raw_sampling_entries,
            options.distance_min,
            options.distance_max);
        auto [group_iter, inserted]{
            context.selected_group_id_by_key.emplace(
                atom_context.group_key,
                context.selected_atom_index_list_by_group.size())
        };
        if (inserted)
        {
            context.selected_atom_index_list_by_group.emplace_back();
        }
        atom_context.group_id = group_iter->second;
        context.selected_atom_index_list_by_group.at(atom_context.group_id)
            .emplace_back(atom_index);
    }

    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        auto & atom_context{ context.at(atom_index) };
        const auto * atom{ atom_context.atom };
        const auto neighbor_atom_list{ atom->FindNeighborAtoms(kSecondStageNeighborAtomSearchRange) };
        std::unordered_set<const AtomObject *> neighbor_atom_set;

        atom_context.sample_neighbor_offset_list.reserve(
            atom_context.raw_sampling_entries.size() + 1);
        atom_context.sample_neighbor_offset_list.emplace_back(0);
        atom_context.sample_neighbor_list.reserve(
            atom_context.raw_sampling_entries.size() *
            neighbor_atom_list.size());
        for (std::size_t sample_index = 0;
            sample_index < atom_context.raw_sampling_entries.size();
            sample_index++)
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
                if (distance > kSecondStageNeighborContributionDistanceMax) continue;
                neighbor_atom_set.emplace(neighbor_atom);

                const auto selected_iter{ atom_index_map.find(neighbor_atom) };
                if (selected_iter != atom_index_map.end())
                {
                    atom_context.sample_neighbor_list.emplace_back(
                        SecondStageNeighborSample{
                            true,
                            selected_iter->second,
                            distance
                        });
                    continue;
                }

                auto contributor_iter{ unselected_atom_index_map.find(neighbor_atom) };
                if (contributor_iter == unselected_atom_index_map.end())
                {
                    const auto contributor_index{
                        context.unselected_atom_list.size()
                    };
                    const auto group_key{
                        data_internal::GetGroupKey(neighbor_atom)
                    };
                    const auto selected_group_iter{
                        context.selected_group_id_by_key.find(group_key)
                    };
                    context.unselected_atom_list.emplace_back(
                        SecondStageUnselectedContributor{
                            neighbor_atom,
                            group_key,
                            selected_group_iter ==
                                context.selected_group_id_by_key.end() ?
                                std::nullopt :
                                std::optional<std::size_t>{
                                    selected_group_iter->second }
                        });
                    contributor_iter = unselected_atom_index_map.emplace(
                        neighbor_atom,
                        contributor_index).first;
                }
                atom_context.sample_neighbor_list.emplace_back(
                    SecondStageNeighborSample{
                        false,
                        contributor_iter->second,
                        distance
                    });
            }
            atom_context.sample_neighbor_offset_list.emplace_back(
                atom_context.sample_neighbor_list.size());
        }
        atom_context.neighbor_count_for_peeling =
            static_cast<int>(neighbor_atom_set.size());
    }

    return context;
}

inline void StoreSecondStageNeighborCounts(
    ModelObject & model_object,
    const SecondStageLocalFittingContext & context)
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
    const auto median_model{
        BuildLocalFittingGaussianParameterMedian(model_list)
    };
    if (!median_model.has_value()) return std::nullopt;
    return GaussianModel3DWithUncertainty{
        *median_model,
        GaussianModel3DUncertainty{}
    };
}

inline std::optional<SecondStageInitialStateBuildResult> BuildInitialLocalFittingState(
    SecondStageLocalFittingContext & context,
    bool & unselected_seed_failure)
{
    unselected_seed_failure = false;
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

        models_by_group[group_key].emplace_back(
            direct_selection->model.GetModel());
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
        if (!selection.has_value()) return std::nullopt;

        result.mdpde = selection->model;
        build_result.selection_record_list.emplace_back(
            SecondStageSeedSelectionRecord{
                i,
                selection->source,
                original_model,
                selection->model.GetModel()
            });
    }

    for (std::size_t i = 0;
        i < context.unselected_atom_list.size();
        i++)
    {
        auto & contributor{ context.unselected_atom_list.at(i) };
        std::optional<GaussianModel3DWithUncertainty> group_median;
        const auto group_median_iter{
            median_by_group.find(contributor.group_key)
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
            unselected_seed_failure = true;
            return std::nullopt;
        }

        contributor.initial_seed = selection->model;
        contributor.seed_source = selection->source;
        build_result.unselected_selection_record_list.emplace_back(
            SecondStageSeedSelectionRecord{
                i,
                selection->source,
                GaussianModel3D{},
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



} // namespace rhbm_gem::core::detail
