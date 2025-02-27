#include "PotentialAnalysisVisitor.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "MapObject.hpp"
#include "DataObjectManager.hpp"
#include "AtomSelector.hpp"
#include "SphereSampler.hpp"
#include "MapInterpolationVisitor.hpp"
#include "HRLModelHelper.hpp"
#include "GausLinearTransformHelper.hpp"
#include "ScopeTimer.hpp"
#include "AtomicPotentialEntry.hpp"
#include "GroupPotentialEntryBase.hpp"
#include "GroupPotentialEntry.hpp"
#include "GroupPotentialEntryFactory.hpp"

PotentialAnalysisVisitor::PotentialAnalysisVisitor(
    std::shared_ptr<AtomSelector> atom_selector,
    std::shared_ptr<SphereSampler> sphere_sampler) :
    m_alpha_r{ 0.0 }, m_alpha_g{ 0.0 },
    m_x_min{ 0.0 }, m_x_max{ 0.0 },
    m_atom_selector{ std::move(atom_selector) },
    m_sphere_sampler{ std::move(sphere_sampler) }
{

}

PotentialAnalysisVisitor::~PotentialAnalysisVisitor()
{

}

void PotentialAnalysisVisitor::VisitAtomObject(AtomObject * data_object)
{
    bool selected_flag
    {
        m_atom_selector->GetSelectionFlag(
            data_object->GetChainID(),
            data_object->GetIndicator(),
            data_object->GetResidue(),
            data_object->GetElement(),
            data_object->GetRemoteness(),
            data_object->GetBranch()
        )
    };
    data_object->SetSelectedFlag(selected_flag);
}

void PotentialAnalysisVisitor::VisitModelObject(ModelObject * data_object)
{
    ScopeTimer timer("PotentialAnalysisVisitor::VisitModelObject");
    std::cout <<"- Visiting ModelObject..." << std::endl;
    auto selected_atom_size{ data_object->GetNumberOfSelectedAtom() };
    const auto & atom_list{ data_object->GetComponentsList() };
    m_selected_atom_list.clear();
    m_selected_atom_list.reserve(selected_atom_size);
    for (auto & atom : atom_list)
    {
        if (atom->GetSelectedFlag() == false) continue;
        auto potential_entry{ std::make_unique<AtomicPotentialEntry>() };
        atom->AddAtomicPotentialEntry(std::move(potential_entry));
        m_selected_atom_list.emplace_back(atom.get());
    }
    std::cout <<" Number of selected atom = "<< selected_atom_size << std::endl;
}

void PotentialAnalysisVisitor::VisitMapObject(MapObject * data_object)
{
    ScopeTimer timer("PotentialAnalysisVisitor::VisitMapObject");
    std::cout <<"- Visiting MapObject..." << std::endl;
    MapInterpolationVisitor visitor{ m_sphere_sampler };
    for (auto & atom : m_selected_atom_list)
    {
        auto entry{ atom->GetAtomicPotentialEntry() };
        visitor.SetPosition(atom->GetPosition());
        data_object->Accept(&visitor);
        entry->AddDistanceAndMapValueList(visitor.GetSamplingDataList());
    }
}

void PotentialAnalysisVisitor::Analysis(DataObjectManager * data_manager)
{
    std::cout <<"- Analysis..." << std::endl;
    
    const auto & model_object{ data_manager->GetDataObjectRef(m_model_key_tag) };
    model_object->Accept(this);

    const auto & map_object{ data_manager->GetDataObjectRef(m_map_key_tag) };
    map_object->Accept(this);

    RunAtomClassification("residue_group", dynamic_cast<ModelObject*>(model_object.get()));
    RunPotentialFitting("residue_group", dynamic_cast<ModelObject*>(model_object.get()));
}

void PotentialAnalysisVisitor::RunAtomClassification(const std::string & key, ModelObject * model_object)
{
    ScopeTimer timer("PotentialAnalysisVisitor::RunAtomClassification");
    std::cout <<"- RunAtomClassification..." << std::endl;
    auto group_potential_entry( GroupPotentialEntryFactory::Create(key) );
    for (auto atom : m_selected_atom_list)
    {
        auto group_type{ std::make_tuple(atom->GetResidue(), atom->GetElement(), atom->GetRemoteness(), atom->GetBranch(), atom->GetSpecialAtomFlag()) };
        group_potential_entry->AddAtomObjectPtr(&group_type, atom);
        group_type_set.insert(group_type);
    }
    model_object->AddGroupPotentialEntry(key, std::unique_ptr<GroupPotentialEntryBase>(group_potential_entry.release()));
}

void PotentialAnalysisVisitor::RunPotentialFitting(const std::string & key, ModelObject * model_object)
{
    ScopeTimer timer("PotentialAnalysisVisitor::RunPotentialFitting");
    std::cout <<"- RunPotentialFitting..." << std::endl;
    auto group_potential_entry{ model_object->GetGroupPotentialEntry(key) };

    for (auto & group_type : group_type_set)
    {
        auto atom_list{ group_potential_entry->GetAtomObjectPtrList(&group_type) };
        auto group_size{ static_cast<int>(atom_list.size()) };
        std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> data_array;
        data_array.reserve(group_size);
        for (auto atom : atom_list)
        {
            auto entry{ atom->GetAtomicPotentialEntry() };
            std::vector<Eigen::VectorXd> sampling_entry_list;
            sampling_entry_list.reserve(entry->GetDistanceAndMapValueListSize());
            for (auto & data_entry : entry->GetDistanceAndMapValueList())
            {
                auto gaus_x{ static_cast<double>(std::get<0>(data_entry)) };
                auto gaus_y{ static_cast<double>(std::get<1>(data_entry)) };
                if (gaus_x < m_x_min || gaus_x > m_x_max) continue;
                if (gaus_y <= 0.0) continue;
                sampling_entry_list.emplace_back(GausLinearTransformHelper::BuildLinearModelDataVector(gaus_x, gaus_y));
            }
            data_array.emplace_back(std::make_tuple(sampling_entry_list, atom->GetInfo()));
        }
        auto model_estimator{ std::make_unique<HRLModelHelper>(2, group_size) };
        model_estimator->SetDataArray(data_array);
        model_estimator->RunEstimation(m_alpha_r, m_alpha_g);

        auto gaus_prior{
            GausLinearTransformHelper::BuildGausModelWithVariance(
                model_estimator->GetMuVectorPrior(), model_estimator->GetCapitalLambdaMatrix())
        };

        group_potential_entry->AddGausEstimateMean(&group_type, GausLinearTransformHelper::BuildGausModel(model_estimator->GetMuVectorMean()));
        group_potential_entry->AddGausEstimateMDPDE(&group_type, GausLinearTransformHelper::BuildGausModel(model_estimator->GetMuVectorMDPDE()));
        group_potential_entry->AddGausEstimatePrior(&group_type, std::get<0>(gaus_prior));
        group_potential_entry->AddGausVariancePrior(&group_type, std::get<1>(gaus_prior));

        auto count{ 0 };
        for (auto atom : atom_list)
        {
            auto atom_entry{ atom->GetAtomicPotentialEntry() };
            Eigen::VectorXd beta_vector_ols{ model_estimator->GetBetaMatrixOLS(count) };
            Eigen::VectorXd beta_vector_mdpde{ model_estimator->GetBetaMatrixMDPDE(count) };
            Eigen::VectorXd beta_vector_posterior{ model_estimator->GetBetaMatrixPosterior(count) };
            auto sigma_matrix_posterior{ model_estimator->GetCapitalSigmaMatrixPosterior(count) };
            auto gaus_posterior{ GausLinearTransformHelper::BuildGausModelWithVariance(beta_vector_posterior, sigma_matrix_posterior) };
            atom_entry->AddGausEstimateOLS(GausLinearTransformHelper::BuildGausModel(beta_vector_ols));
            atom_entry->AddGausEstimateMDPDE(GausLinearTransformHelper::BuildGausModel(beta_vector_mdpde));
            atom_entry->AddGausEstimatePosterior(key, std::get<0>(gaus_posterior));
            atom_entry->AddGausVariancePosterior(key, std::get<1>(gaus_posterior));
            atom_entry->AddOutlierTag(key, model_estimator->GetOutlierFlag(count));
            atom_entry->AddStatisticalDistance(key, model_estimator->GetStatisticalDistance(count));
            count++;
        }
    }
}

void PotentialAnalysisVisitor::SetFitRange(double x_min, double x_max)
{
    m_x_min = x_min;
    m_x_max = x_max;
}