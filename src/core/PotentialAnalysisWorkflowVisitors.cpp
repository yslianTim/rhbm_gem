#include "PotentialAnalysisWorkflowVisitors.hpp"

#include <atomic>
#include <memory>
#include <stdexcept>
#include <string>

#include "AtomClassifier.hpp"
#include "AtomObject.hpp"
#include "BondClassifier.hpp"
#include "BondObject.hpp"
#include "ChemicalDataHelper.hpp"
#include "CylinderSampler.hpp"
#include "DataObjectWorkflowVisitors.hpp"
#include "GausLinearTransformHelper.hpp"
#include "GroupPotentialEntry.hpp"
#include "HRLDataTransform.hpp"
#include "HRLModelAlgorithms.hpp"
#include "LocalPotentialEntry.hpp"
#include "Logger.hpp"
#include "MapObject.hpp"
#include "ModelObject.hpp"
#include "PotentialAnalysisCommand.hpp"
#include "PotentialAnalysisExecutionOptions.hpp"
#include "ScopeTimer.hpp"
#include "SphereSampler.hpp"

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace rhbm_gem {
namespace {

[[noreturn]] void ThrowUnsupportedType(const char * visitor_name, const char * supported_type)
{
    throw std::logic_error(
        std::string(visitor_name) + " supports " + supported_type + " only.");
}

} // namespace

AtomSamplingVisitor::AtomSamplingVisitor(
    const MapObject & map_object,
    const PotentialAnalysisCommandOptions & options) :
    m_map_object{ &map_object },
    m_options{ &options }
{
}

void AtomSamplingVisitor::VisitAtomObject(AtomObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("AtomSamplingVisitor", "ModelObject");
}

void AtomSamplingVisitor::VisitBondObject(BondObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("AtomSamplingVisitor", "ModelObject");
}

void AtomSamplingVisitor::VisitModelObject(ModelObject & data_object)
{
    if (m_map_object == nullptr || m_options == nullptr)
    {
        throw std::logic_error("AtomSamplingVisitor requires map and options.");
    }
    ScopeTimer timer("PotentialAnalysisCommand::RunAtomMapValueSampling");
    auto sampler{ std::make_unique<SphereSampler>() };
    sampler->SetSamplingSize(static_cast<unsigned int>(m_options->sampling_size));
    sampler->SetDistanceRangeMinimum(m_options->sampling_range_min);
    sampler->SetDistanceRangeMaximum(m_options->sampling_range_max);
    sampler->Print();

    const auto & atom_list{ data_object.GetSelectedAtomList() };
    const auto atom_size{ atom_list.size() };
    size_t atom_count{ 0 };

#ifdef USE_OPENMP
    #pragma omp parallel num_threads(m_options->thread_size)
    {
        MapSamplingWorkflow sampling_workflow{ sampler.get() };
        #pragma omp for
        for (size_t i = 0; i < atom_size; i++)
        {
            auto atom{ atom_list[i] };
            auto entry{ atom->GetLocalPotentialEntry() };
            entry->AddDistanceAndMapValueList(
                sampling_workflow.Sample(*m_map_object, atom->GetPosition()));
            entry->AddBasisAndResponseEntryList(
                GausLinearTransformHelper::MapValueTransform(
                    entry->GetDistanceAndMapValueList(),
                    m_options->fit_range_min,
                    m_options->fit_range_max));
            if (!m_options->use_training_alpha)
            {
                entry->SetAlphaR(m_options->alpha_r);
            }
            #pragma omp critical
            {
                atom_count++;
                Logger::ProgressPercent(atom_count, atom_size);
            }
        }
    }
#else
    MapSamplingWorkflow sampling_workflow{ sampler.get() };
    for (size_t i = 0; i < atom_size; i++)
    {
        auto atom{ atom_list[i] };
        auto entry{ atom->GetLocalPotentialEntry() };
        entry->AddDistanceAndMapValueList(
            sampling_workflow.Sample(*m_map_object, atom->GetPosition()));
        entry->AddBasisAndResponseEntryList(
            GausLinearTransformHelper::MapValueTransform(
                entry->GetDistanceAndMapValueList(),
                m_options->fit_range_min,
                m_options->fit_range_max));
        if (!m_options->use_training_alpha)
        {
            entry->SetAlphaR(m_options->alpha_r);
        }
        atom_count++;
        Logger::ProgressPercent(atom_count, atom_size);
    }
#endif
}

void AtomSamplingVisitor::VisitMapObject(MapObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("AtomSamplingVisitor", "ModelObject");
}

void AtomGroupingVisitor::VisitAtomObject(AtomObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("AtomGroupingVisitor", "ModelObject");
}

void AtomGroupingVisitor::VisitBondObject(BondObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("AtomGroupingVisitor", "ModelObject");
}

void AtomGroupingVisitor::VisitModelObject(ModelObject & data_object)
{
    ScopeTimer timer("RunAtomGroupClassification");
    Logger::Log(LogLevel::Info, "Atom Classification Summary:");
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        auto group_potential_entry{ std::make_unique<GroupPotentialEntry>() };
        for (auto atom : data_object.GetSelectedAtomList())
        {
            auto group_key{ AtomClassifier::GetGroupKeyInClass(atom, class_key) };
            group_potential_entry->AddAtomObjectPtr(group_key, atom);
            group_potential_entry->InsertGroupKey(group_key);
        }
        const auto group_size{ group_potential_entry->GetGroupKeySet().size() };
        data_object.AddAtomGroupPotentialEntry(class_key, group_potential_entry);
        Logger::Log(
            LogLevel::Info,
            " - Class type: " + class_key + " include " + std::to_string(group_size)
                + " groups.");
    }
}

void AtomGroupingVisitor::VisitMapObject(MapObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("AtomGroupingVisitor", "ModelObject");
}

LocalFittingVisitor::LocalFittingVisitor(
    const PotentialAnalysisCommandOptions & options,
    double universal_alpha_r) :
    m_options{ &options },
    m_universal_alpha_r{ universal_alpha_r }
{
}

void LocalFittingVisitor::VisitAtomObject(AtomObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("LocalFittingVisitor", "ModelObject");
}

void LocalFittingVisitor::VisitBondObject(BondObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("LocalFittingVisitor", "ModelObject");
}

void LocalFittingVisitor::VisitModelObject(ModelObject & data_object)
{
    if (m_options == nullptr)
    {
        throw std::logic_error("LocalFittingVisitor requires options.");
    }

    ScopeTimer timer("PotentialAnalysisCommand::RunLocalAtomFitting");
    std::atomic<size_t> atom_count{ 0 };
    auto & selected_atom_list{ data_object.GetSelectedAtomList() };
    const auto selected_atom_size{ selected_atom_list.size() };
    Logger::Log(
        LogLevel::Info,
        "Run Local atom fitting for " + std::to_string(selected_atom_size) + " atoms.");
#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(m_options->thread_size)
#endif
    for (size_t i = 0; i < selected_atom_size; i++)
    {
        auto local_entry{ selected_atom_list[i]->GetLocalPotentialEntry() };
        auto & data_entry_list{ local_entry->GetBasisAndResponseEntryList() };
        const auto dataset{ HRLDataTransform::BuildMemberDataset(data_entry_list) };
        const auto result{
            HRLModelAlgorithms::EstimateBetaMDPDE(
                m_universal_alpha_r,
                dataset.X,
                dataset.y,
                detail::MakePotentialAnalysisExecutionOptions(*m_options, true))
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

void LocalFittingVisitor::VisitMapObject(MapObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("LocalFittingVisitor", "ModelObject");
}

BondSamplingVisitor::BondSamplingVisitor(
    const MapObject & map_object,
    const PotentialAnalysisCommandOptions & options) :
    m_map_object{ &map_object },
    m_options{ &options }
{
}

void BondSamplingVisitor::VisitAtomObject(AtomObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("BondSamplingVisitor", "ModelObject");
}

void BondSamplingVisitor::VisitBondObject(BondObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("BondSamplingVisitor", "ModelObject");
}

void BondSamplingVisitor::VisitModelObject(ModelObject & data_object)
{
    if (m_map_object == nullptr || m_options == nullptr)
    {
        throw std::logic_error("BondSamplingVisitor requires map and options.");
    }
    ScopeTimer timer("PotentialAnalysisBondWorkflow::RunBondMapValueSampling");
    auto sampler{ std::make_unique<CylinderSampler>() };
    sampler->SetSamplingSize(static_cast<unsigned int>(m_options->sampling_size));
    sampler->SetDistanceRangeMinimum(m_options->sampling_range_min);
    sampler->SetDistanceRangeMaximum(m_options->sampling_range_max);
    sampler->SetHeight(m_options->sampling_height);
    sampler->Print();

    const auto & bond_list{ data_object.GetSelectedBondList() };
    const auto bond_size{ bond_list.size() };
    size_t bond_count{ 0 };

#ifdef USE_OPENMP
    #pragma omp parallel num_threads(m_options->thread_size)
    {
        MapSamplingWorkflow sampling_workflow{ sampler.get() };
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
                sampling_workflow.Sample(
                    *m_map_object,
                    adjusted_position,
                    bond_vector));
            entry->AddBasisAndResponseEntryList(
                GausLinearTransformHelper::MapValueTransform(
                    entry->GetDistanceAndMapValueList(),
                    m_options->fit_range_min,
                    m_options->fit_range_max));
            entry->SetAlphaR(m_options->alpha_r);
            #pragma omp critical
            {
                bond_count++;
                Logger::ProgressPercent(bond_count, bond_size);
            }
        }
    }
#else
    MapSamplingWorkflow sampling_workflow{ sampler.get() };
    for (size_t i = 0; i < bond_size; i++)
    {
        auto bond{ bond_list[i] };
        auto entry{ bond->GetLocalPotentialEntry() };
        entry->AddDistanceAndMapValueList(
            sampling_workflow.Sample(
                *m_map_object,
                bond->GetPosition(),
                bond->GetBondVector()));
        entry->AddBasisAndResponseEntryList(
            GausLinearTransformHelper::MapValueTransform(
                entry->GetDistanceAndMapValueList(),
                m_options->fit_range_min,
                m_options->fit_range_max));
        entry->SetAlphaR(m_options->alpha_r);
        bond_count++;
        Logger::ProgressPercent(bond_count, bond_size);
    }
#endif
}

void BondSamplingVisitor::VisitMapObject(MapObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("BondSamplingVisitor", "ModelObject");
}

void BondGroupingVisitor::VisitAtomObject(AtomObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("BondGroupingVisitor", "ModelObject");
}

void BondGroupingVisitor::VisitBondObject(BondObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("BondGroupingVisitor", "ModelObject");
}

void BondGroupingVisitor::VisitModelObject(ModelObject & data_object)
{
    ScopeTimer timer("PotentialAnalysisBondWorkflow::RunBondGroupClassification");
    Logger::Log(LogLevel::Info, "Bond Classification Summary:");
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupBondClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupBondClassKey(i) };
        auto group_potential_entry{ std::make_unique<GroupPotentialEntry>() };
        for (auto bond : data_object.GetSelectedBondList())
        {
            auto group_key{ BondClassifier::GetGroupKeyInClass(bond, class_key) };
            group_potential_entry->AddBondObjectPtr(group_key, bond);
            group_potential_entry->InsertGroupKey(group_key);
        }
        const auto group_size{ group_potential_entry->GetGroupKeySet().size() };
        data_object.AddBondGroupPotentialEntry(class_key, group_potential_entry);
        Logger::Log(
            LogLevel::Info,
            " - Class type: " + class_key + " include "
                + std::to_string(group_size) + " groups.");
    }
}

void BondGroupingVisitor::VisitMapObject(MapObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("BondGroupingVisitor", "ModelObject");
}

} // namespace rhbm_gem
