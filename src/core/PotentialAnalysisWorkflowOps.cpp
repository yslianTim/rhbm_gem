#include "PotentialAnalysisWorkflowOps.hpp"

#include <array>
#include <atomic>
#include <memory>
#include <string>
#include <vector>

#include "AtomClassifier.hpp"
#include "AtomObject.hpp"
#include "BondClassifier.hpp"
#include "BondObject.hpp"
#include "ChemicalDataHelper.hpp"
#include "CylinderSampler.hpp"
#include "GausLinearTransformHelper.hpp"
#include "GroupPotentialEntry.hpp"
#include "HRLDataTransform.hpp"
#include "HRLModelAlgorithms.hpp"
#include "LocalPotentialEntry.hpp"
#include "Logger.hpp"
#include "MapObject.hpp"
#include "MapSampling.hpp"
#include "ModelObject.hpp"
#include "PotentialAnalysisCommand.hpp"
#include "PotentialAnalysisExecutionOptions.hpp"
#include "ScopeTimer.hpp"
#include "SphereSampler.hpp"

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace rhbm_gem {

void RunAtomSampling(
    ModelObject & model_object,
    const MapObject & map_object,
    const PotentialAnalysisCommandOptions & options)
{
    ScopeTimer timer("PotentialAnalysisCommand::RunAtomMapValueSampling");
    auto sampler{ std::make_unique<SphereSampler>() };
    sampler->SetSamplingSize(static_cast<unsigned int>(options.sampling_size));
    sampler->SetDistanceRangeMinimum(options.sampling_range_min);
    sampler->SetDistanceRangeMaximum(options.sampling_range_max);
    sampler->Print();

    const auto & atom_list{ model_object.GetSelectedAtomList() };
    const auto atom_size{ atom_list.size() };
    size_t atom_count{ 0 };

#ifdef USE_OPENMP
    #pragma omp parallel num_threads(options.thread_size)
    {
        #pragma omp for
        for (size_t i = 0; i < atom_size; i++)
        {
            auto atom{ atom_list[i] };
            auto entry{ atom->GetLocalPotentialEntry() };
            entry->AddDistanceAndMapValueList(
                SampleMapValues(map_object, *sampler, atom->GetPosition()));
            entry->AddBasisAndResponseEntryList(
                GausLinearTransformHelper::MapValueTransform(
                    entry->GetDistanceAndMapValueList(),
                    options.fit_range_min,
                    options.fit_range_max));
            if (!options.use_training_alpha)
            {
                entry->SetAlphaR(options.alpha_r);
            }
            #pragma omp critical
            {
                atom_count++;
                Logger::ProgressPercent(atom_count, atom_size);
            }
        }
    }
#else
    for (size_t i = 0; i < atom_size; i++)
    {
        auto atom{ atom_list[i] };
        auto entry{ atom->GetLocalPotentialEntry() };
        entry->AddDistanceAndMapValueList(
            SampleMapValues(map_object, *sampler, atom->GetPosition()));
        entry->AddBasisAndResponseEntryList(
            GausLinearTransformHelper::MapValueTransform(
                entry->GetDistanceAndMapValueList(),
                options.fit_range_min,
                options.fit_range_max));
        if (!options.use_training_alpha)
        {
            entry->SetAlphaR(options.alpha_r);
        }
        atom_count++;
        Logger::ProgressPercent(atom_count, atom_size);
    }
#endif
}

void RunAtomGrouping(ModelObject & model_object)
{
    ScopeTimer timer("RunAtomGroupClassification");
    Logger::Log(LogLevel::Info, "Atom Classification Summary:");
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        auto group_potential_entry{ std::make_unique<GroupPotentialEntry>() };
        for (auto atom : model_object.GetSelectedAtomList())
        {
            auto group_key{ AtomClassifier::GetGroupKeyInClass(atom, class_key) };
            group_potential_entry->AddAtomObjectPtr(group_key, atom);
            group_potential_entry->InsertGroupKey(group_key);
        }
        const auto group_size{ group_potential_entry->GetGroupKeySet().size() };
        model_object.AddAtomGroupPotentialEntry(class_key, group_potential_entry);
        Logger::Log(
            LogLevel::Info,
            " - Class type: " + class_key + " include " + std::to_string(group_size)
                + " groups.");
    }
}

void RunLocalAtomFitting(
    ModelObject & model_object,
    const PotentialAnalysisCommandOptions & options,
    double universal_alpha_r)
{
    ScopeTimer timer("PotentialAnalysisCommand::RunLocalAtomFitting");
    std::atomic<size_t> atom_count{ 0 };
    auto & selected_atom_list{ model_object.GetSelectedAtomList() };
    const auto selected_atom_size{ selected_atom_list.size() };
    Logger::Log(
        LogLevel::Info,
        "Run Local atom fitting for " + std::to_string(selected_atom_size) + " atoms.");
#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(options.thread_size)
#endif
    for (size_t i = 0; i < selected_atom_size; i++)
    {
        auto local_entry{ selected_atom_list[i]->GetLocalPotentialEntry() };
        auto & data_entry_list{ local_entry->GetBasisAndResponseEntryList() };
        const auto dataset{ HRLDataTransform::BuildMemberDataset(data_entry_list) };
        const auto result{
            HRLModelAlgorithms::EstimateBetaMDPDE(
                universal_alpha_r,
                dataset.X,
                dataset.y,
                detail::MakePotentialAnalysisExecutionOptions(options, true))
        };

        local_entry->SetBetaEstimateOLS(result.beta_ols);
        local_entry->SetBetaEstimateMDPDE(result.beta_mdpde);
        local_entry->SetSigmaSquare(result.sigma_square);
        local_entry->SetDataWeight(result.data_weight);
        local_entry->SetDataCovariance(result.data_covariance);

        Eigen::VectorXd model_par_init{ Eigen::VectorXd::Zero(3) };
        model_par_init(0) = local_entry->GetMomentZeroEstimate();
        model_par_init(1) = local_entry->GetMomentTwoEstimate();
        auto gaus_ols{
            GausLinearTransformHelper::BuildGaus3DModel(result.beta_ols, model_par_init)
        };
        auto gaus_mdpde{
            GausLinearTransformHelper::BuildGaus3DModel(result.beta_mdpde, model_par_init)
        };
        local_entry->AddGausEstimateOLS(gaus_ols(0), gaus_ols(1));
        local_entry->AddGausEstimateMDPDE(gaus_mdpde(0), gaus_mdpde(1));

#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            atom_count++;
            Logger::ProgressPercent(atom_count, selected_atom_size);
        }
    }
}

void RunBondSampling(
    ModelObject & model_object,
    const MapObject & map_object,
    const PotentialAnalysisCommandOptions & options)
{
    ScopeTimer timer("PotentialAnalysisBondWorkflow::RunBondMapValueSampling");
    auto sampler{ std::make_unique<CylinderSampler>() };
    sampler->SetSamplingSize(static_cast<unsigned int>(options.sampling_size));
    sampler->SetDistanceRangeMinimum(options.sampling_range_min);
    sampler->SetDistanceRangeMaximum(options.sampling_range_max);
    sampler->SetHeight(options.sampling_height);
    sampler->Print();

    const auto & bond_list{ model_object.GetSelectedBondList() };
    const auto bond_size{ bond_list.size() };
    size_t bond_count{ 0 };

#ifdef USE_OPENMP
    #pragma omp parallel num_threads(options.thread_size)
    {
        #pragma omp for
        for (size_t i = 0; i < bond_size; i++)
        {
            auto bond{ bond_list[i] };
            auto entry{ bond->GetLocalPotentialEntry() };
            auto bond_vector{ bond->GetBondVector() };
            auto bond_position{ bond->GetPosition() };
            constexpr float adjusted_rate{ 0.0f };
            std::array<float, 3> adjusted_position{
                bond_position[0] + 0.5f * bond_vector[0] * adjusted_rate,
                bond_position[1] + 0.5f * bond_vector[1] * adjusted_rate,
                bond_position[2] + 0.5f * bond_vector[2] * adjusted_rate
            };
            entry->AddDistanceAndMapValueList(
                SampleMapValues(
                    map_object,
                    *sampler,
                    adjusted_position,
                    bond_vector));
            entry->AddBasisAndResponseEntryList(
                GausLinearTransformHelper::MapValueTransform(
                    entry->GetDistanceAndMapValueList(),
                    options.fit_range_min,
                    options.fit_range_max));
            entry->SetAlphaR(options.alpha_r);
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
        auto entry{ bond->GetLocalPotentialEntry() };
        entry->AddDistanceAndMapValueList(
            SampleMapValues(
                map_object,
                *sampler,
                bond->GetPosition(),
                bond->GetBondVector()));
        entry->AddBasisAndResponseEntryList(
            GausLinearTransformHelper::MapValueTransform(
                entry->GetDistanceAndMapValueList(),
                options.fit_range_min,
                options.fit_range_max));
        entry->SetAlphaR(options.alpha_r);
        bond_count++;
        Logger::ProgressPercent(bond_count, bond_size);
    }
#endif
}

void RunBondGrouping(ModelObject & model_object)
{
    ScopeTimer timer("PotentialAnalysisBondWorkflow::RunBondGroupClassification");
    Logger::Log(LogLevel::Info, "Bond Classification Summary:");
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupBondClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupBondClassKey(i) };
        auto group_potential_entry{ std::make_unique<GroupPotentialEntry>() };
        for (auto bond : model_object.GetSelectedBondList())
        {
            auto group_key{ BondClassifier::GetGroupKeyInClass(bond, class_key) };
            group_potential_entry->AddBondObjectPtr(group_key, bond);
            group_potential_entry->InsertGroupKey(group_key);
        }
        const auto group_size{ group_potential_entry->GetGroupKeySet().size() };
        model_object.AddBondGroupPotentialEntry(class_key, group_potential_entry);
        Logger::Log(
            LogLevel::Info,
            " - Class type: " + class_key + " include "
                + std::to_string(group_size) + " groups.");
    }
}

} // namespace rhbm_gem
