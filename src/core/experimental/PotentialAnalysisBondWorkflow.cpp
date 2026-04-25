#include "experimental/PotentialAnalysisBondWorkflow.hpp"

#include "command/detail/MapSampling.hpp"
#include "command/PotentialAnalysisCommand.hpp"
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/math/CylinderSampler.hpp>

#include <array>
#include <atomic>
#include <memory>
#include <tuple>
#include <unordered_map>
#include <vector>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace rhbm_gem::experimental {

namespace {
namespace ls = rhbm_gem::linearization_service;


struct PotentialAnalysisBondWorkflowContext
{
    ModelObject & model_object;
    MapObject & map_object;
    const PotentialAnalysisRequest & options;
    int thread_size;
};

RHBMExecutionOptions MakePotentialAnalysisExecutionOptions(
    int thread_size,
    bool quiet_mode)
{
    RHBMExecutionOptions execution_options;
    execution_options.quiet_mode = quiet_mode;
    execution_options.thread_size = thread_size;
    return execution_options;
}

const ls::LinearizationSpec & GaussianDatasetSpec()
{
    static const auto spec{ ls::LinearizationSpec::DefaultDataset() };
    return spec;
}

ls::LinearizationContext BuildLocalLinearizationContext(const LocalPotentialView & view)
{
    Eigen::VectorXd model_par_init{ Eigen::VectorXd::Zero(3) };
    model_par_init(0) = view.GetMomentZeroEstimate();
    model_par_init(1) = view.GetMomentTwoEstimate();
    return ls::LinearizationContext::FromModelParameters(model_par_init);
}

void RunBondSampling(
    ModelObject & model_object,
    const MapObject & map_object,
    const PotentialAnalysisRequest & options,
    int thread_size)
{
    ScopeTimer timer("PotentialAnalysisBondWorkflow::RunBondMapValueSampling");
    const auto & dataset_spec{ GaussianDatasetSpec() };
    CylinderSampler sampler;
    sampler.SetSampleCount(static_cast<unsigned int>(options.sampling_size));
    sampler.SetDistanceRange(options.sampling_range_min, options.sampling_range_max);
    sampler.SetHeight(options.sampling_height);
    sampler.Print();

    const auto & bond_list{ model_object.GetSelectedBonds() };
    const auto bond_size{ bond_list.size() };
    size_t bond_count{ 0 };
    auto analysis{ model_object.EditAnalysis() };
    std::vector<MutableLocalPotentialView> local_entry_list;
    local_entry_list.reserve(bond_size);
    for (auto * bond : bond_list)
    {
        local_entry_list.emplace_back(analysis.EnsureBondLocalPotential(*bond));
    }

#ifdef USE_OPENMP
    #pragma omp parallel num_threads(thread_size)
    {
        #pragma omp for
        for (size_t i = 0; i < bond_size; i++)
        {
            auto bond{ bond_list[i] };
            auto entry{ local_entry_list[i] };
            auto bond_vector{ bond->GetBondVector() };
            auto bond_position{ bond->GetPosition() };
            constexpr float adjusted_rate{ 0.0f };
            std::array<float, 3> adjusted_position{
                bond_position[0] + 0.5f * bond_vector[0] * adjusted_rate,
                bond_position[1] + 0.5f * bond_vector[1] * adjusted_rate,
                bond_position[2] + 0.5f * bond_vector[2] * adjusted_rate
            };
            auto sampling_entries{
                SampleMapValues(
                    map_object,
                    sampler,
                    adjusted_position,
                    bond_vector)
            };
            entry.SetSamplingEntries(sampling_entries);
            const auto local_view{ LocalPotentialView::RequireFor(*bond) };
            entry.SetDataset(
                ls::BuildDataset(
                    dataset_spec,
                    sampling_entries,
                    options.fit_range_min,
                    options.fit_range_max,
                    BuildLocalLinearizationContext(local_view))
            );
            entry.SetAlphaR(options.alpha_r);
            #pragma omp critical
            {
                bond_count++;
                Logger::ProgressPercent(bond_count, bond_size);
            }
        }
    }
#else
    for (size_t i = 0; i < bond_size; i++)
    {
        auto bond{ bond_list[i] };
        auto entry{ local_entry_list[i] };
        auto sampling_entries{
            SampleMapValues(
                map_object,
                sampler,
                bond->GetPosition(),
                bond->GetBondVector())
        };
        entry.SetSamplingEntries(sampling_entries);
        const auto local_view{ LocalPotentialView::RequireFor(*bond) };
        entry.SetDataset(
            ls::BuildDataset(
                dataset_spec,
                sampling_entries,
                options.fit_range_min,
                options.fit_range_max,
                BuildLocalLinearizationContext(local_view))
        );
        entry.SetAlphaR(options.alpha_r);
        bond_count++;
        Logger::ProgressPercent(bond_count, bond_size);
    }
#endif
}

void RunBondGrouping(ModelObject & model_object)
{
    ScopeTimer timer("PotentialAnalysisBondWorkflow::RunBondGroupClassification");
    auto analysis{ model_object.EditAnalysis() };
    analysis.RebuildBondGroupsFromSelection();
    const auto analysis_view{ model_object.GetAnalysisView() };
    Logger::Log(LogLevel::Info, analysis_view.DescribeBondGrouping());
}

void RunBondMapValueSampling(const PotentialAnalysisBondWorkflowContext & context)
{
    RunBondSampling(context.model_object, context.map_object, context.options, context.thread_size);
}

void RunBondGroupClassification(const PotentialAnalysisBondWorkflowContext & context)
{
    RunBondGrouping(context.model_object);
}

void RunLocalBondFitting(
    const PotentialAnalysisBondWorkflowContext & context,
    double universal_alpha_r)
{
    ScopeTimer timer("PotentialAnalysisBondWorkflow::RunLocalBondFitting");

    std::atomic<size_t> bond_count{ 0 };
    const auto & selected_bond_list{ context.model_object.GetSelectedBonds() };
    const auto selected_bond_size{ selected_bond_list.size() };
    auto analysis{ context.model_object.EditAnalysis() };
    std::vector<MutableLocalPotentialView> local_entry_list;
    local_entry_list.reserve(selected_bond_size);
    for (auto * bond : selected_bond_list)
    {
        local_entry_list.emplace_back(analysis.EnsureBondLocalPotential(*bond));
    }
    Logger::Log(
        LogLevel::Info,
        "Run Local bond fitting for " + std::to_string(selected_bond_size) + " bonds.");
#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(context.thread_size)
#endif
    for (size_t i = 0; i < selected_bond_size; i++)
    {
        auto local_entry{ local_entry_list[i] };
        const auto & dataset{ local_entry.GetDataset() };
        const auto result{
            rhbm_gem::rhbm_helper::EstimateBetaMDPDE(
                universal_alpha_r,
                dataset,
                MakePotentialAnalysisExecutionOptions(context.thread_size, true))
        };

        local_entry.SetFitResult(result);

#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            bond_count++;
            Logger::ProgressPercent(bond_count, selected_bond_size);
        }
    }
}

void RunBondPotentialFitting(const PotentialAnalysisBondWorkflowContext & context)
{
    ScopeTimer timer("PotentialAnalysisBondWorkflow::RunBondPotentialFitting");
    auto analysis{ context.model_object.EditAnalysis() };
    const auto analysis_view{ context.model_object.GetAnalysisView() };
    std::unordered_map<const BondObject *, MutableLocalPotentialView> local_entry_map;
    local_entry_map.reserve(context.model_object.GetSelectedBondCount());
    for (auto * bond : context.model_object.GetSelectedBonds())
    {
        local_entry_map.emplace(bond, analysis.EnsureBondLocalPotential(*bond));
    }
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupBondClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupBondClassKey(i) };
        Logger::Log(LogLevel::Info, "Class type: " + class_key);

        auto group_keys{ analysis_view.CollectBondGroupKeys(class_key) };
        const auto group_key_size{ group_keys.size() };
        std::atomic<size_t> key_count{ 0 };

#ifdef USE_OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(context.thread_size)
#endif
        for (size_t idx = 0; idx < group_key_size; idx++)
        {
            auto group_key{ group_keys[idx] };
            const auto & bond_list{ analysis_view.GetBondObjectList(group_key, class_key) };
            const auto group_size{ bond_list.size() };
            std::vector<RHBMMemberDataset> member_datasets;
            std::vector<RHBMBetaEstimateResult> member_fit_results;
            member_datasets.reserve(group_size);
            member_fit_results.reserve(group_size);
            for (const auto & bond : bond_list)
            {
                const auto local_entry{ local_entry_map.at(bond) };
                const auto & dataset{ local_entry.GetDataset() };
                member_datasets.emplace_back(dataset);
                member_fit_results.emplace_back(local_entry.GetFitResult());
            }
            const auto input{
                rhbm_gem::rhbm_helper::BuildGroupInput(member_datasets, member_fit_results)
            };
            const auto result{
                rhbm_gem::rhbm_helper::EstimateGroup(
                    context.options.alpha_g,
                    input,
                    MakePotentialAnalysisExecutionOptions(context.thread_size, true))
            };

#ifdef USE_OPENMP
            #pragma omp critical
#endif
            {
                analysis.ApplyBondGroupEstimateResult(
                    group_key,
                    class_key,
                    result,
                    context.options.alpha_g);
                key_count++;
                Logger::ProgressBar(key_count, group_key_size);
            }
        }
    }
}

void RunPotentialAnalysisBondWorkflowImpl(
    const PotentialAnalysisBondWorkflowContext & context)
{
    RunBondMapValueSampling(context);
    RunBondGroupClassification(context);
    RunLocalBondFitting(context, context.options.alpha_r);
    RunBondPotentialFitting(context);
}

} // namespace

void RunPotentialAnalysisBondWorkflow(
    ModelObject & model_object,
    MapObject & map_object,
    const PotentialAnalysisRequest & options,
    int thread_size)
{
    PotentialAnalysisBondWorkflowContext context{
        model_object,
        map_object,
        options,
        thread_size,
    };
    RunPotentialAnalysisBondWorkflowImpl(context);
}

} // namespace rhbm_gem::experimental
