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
#include "PotentialEntry.hpp"

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
        entry->distance_and_map_value_list = visitor.GetSamplingDataList();
    }
}

void PotentialAnalysisVisitor::Analysis(DataObjectManager * data_manager)
{
    std::cout <<"- Analysis..." << std::endl;
    
    const auto & model_object{ data_manager->GetDataObjectRef(m_model_key_tag) };
    model_object->Accept(this);

    const auto & map_object{ data_manager->GetDataObjectRef(m_map_key_tag) };
    map_object->Accept(this);

    RunAtomClassification();
    RunPotentialFitting();
}

void PotentialAnalysisVisitor::RunAtomClassification(void)
{
    ScopeTimer timer("PotentialAnalysisVisitor::RunAtomClassification");
    std::cout <<"- RunAtomClassification..." << std::endl;
    for (auto & atom : m_selected_atom_list)
    {
        auto group_type{ std::make_tuple(atom->GetResidue(), atom->GetElement(), atom->GetRemoteness(), atom->GetBranch()) };
        m_grouping_atom_map[group_type].emplace_back(atom);
    }
}

void PotentialAnalysisVisitor::RunPotentialFitting(void)
{
    ScopeTimer timer("PotentialAnalysisVisitor::RunPotentialFitting");
    std::cout <<"- RunPotentialFitting..." << std::endl;
    for (auto & [tuple_key, atom_list] : m_grouping_atom_map)
    {
        auto group_size{ static_cast<int>(atom_list.size()) };
        std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> data_array;
        data_array.reserve(group_size);
        for (auto atom : atom_list)
        {
            auto entry{ atom->GetAtomicPotentialEntry() };
            std::vector<Eigen::VectorXd> data_entry_list;
            data_entry_list.reserve(entry->distance_and_map_value_list.size());
            for (auto & data_entry : entry->distance_and_map_value_list)
            {
                auto gaus_x{ static_cast<double>(std::get<0>(data_entry)) };
                auto gaus_y{ static_cast<double>(std::get<1>(data_entry)) };
                if (gaus_x < m_x_min || gaus_x > m_x_max) continue;
                if (gaus_y <= 0.0) continue;
                data_entry_list.emplace_back(GausLinearTransformHelper::BuildLinearModelDataVector(gaus_x, gaus_y));
            }
            data_array.emplace_back(std::make_tuple(data_entry_list, atom->GetInfo()));
        }
        auto model_estimator{ std::make_unique<HRLModelHelper>(2, group_size) };
        model_estimator->SetDataArray(data_array);
        model_estimator->RunEstimation(m_alpha_r, m_alpha_g);

        auto gaus_prior{
            GausLinearTransformHelper::BuildGausModelWithVariance(
                model_estimator->GetMuVectorPrior(), model_estimator->GetCapitalLambdaMatrix())
        };
        auto gaus_estimate_prior{ std::get<0>(gaus_prior) };
        auto gaus_variance_prior{ std::get<1>(gaus_prior) };

        auto gaus_mdpde{ GausLinearTransformHelper::BuildGausModel(model_estimator->GetMuVectorMDPDE()) };
        auto gaus_mean{ GausLinearTransformHelper::BuildGausModel(model_estimator->GetMuVectorMean()) };

        auto group_gaus_entry{ std::make_unique<GroupPotentialEntry>() };
        group_gaus_entry->classification_type = 0;
        group_gaus_entry->gaus_estimate_prior = std::make_tuple(gaus_estimate_prior(0), gaus_estimate_prior(1));
        group_gaus_entry->gaus_variance_prior = std::make_tuple(gaus_variance_prior(0), gaus_variance_prior(1));
        group_gaus_entry->gaus_estimate_mdpde = std::make_tuple(gaus_mdpde(0), gaus_mdpde(1));
        group_gaus_entry->gaus_estimate_mean = std::make_tuple(gaus_mean(0), gaus_mean(1));
        m_grouping_gaus_entry_map[tuple_key] = std::move(group_gaus_entry);

        auto count{ 0 };
        for (auto atom : atom_list)
        {
            auto atom_entry{ atom->GetAtomicPotentialEntry() };
            Eigen::VectorXd beta_vector_posterior{ model_estimator->GetBetaMatrixPosterior(count) };
            auto sigma_matrix_posterior{ model_estimator->GetCapitalSigmaMatrixPosterior(count) };
            auto gaus_posterior{ GausLinearTransformHelper::BuildGausModelWithVariance(beta_vector_posterior, sigma_matrix_posterior) };
            auto gaus_estimate_posterior{ std::get<0>(gaus_posterior) };
            auto gaus_variance_posterior{ std::get<1>(gaus_posterior) };
            atom_entry->gaus_estimate_posterior = std::make_tuple(gaus_estimate_posterior(0), gaus_estimate_posterior(1));
            atom_entry->gaus_variance_posterior = std::make_tuple(gaus_variance_posterior(0), gaus_variance_posterior(1));

            Eigen::VectorXd beta_vector_mdpde{ model_estimator->GetBetaMatrixMDPDE(count) };
            auto gaus_mdpde{ GausLinearTransformHelper::BuildGausModel(beta_vector_mdpde) };
            atom_entry->gaus_estimate_mdpde = std::make_tuple(gaus_mdpde(0), gaus_mdpde(1));

            Eigen::VectorXd beta_vector_ols{ model_estimator->GetBetaMatrixOLS(count) };
            auto gaus_ols{ GausLinearTransformHelper::BuildGausModel(beta_vector_ols) };
            atom_entry->gaus_estimate_ols = std::make_tuple(gaus_ols(0), gaus_ols(1));

            atom_entry->outlier_tag = model_estimator->GetOutlierFlag(count);
            atom_entry->statistical_distance = model_estimator->GetStatisticalDistance(count);
            count++;
        }
    }
}

void PotentialAnalysisVisitor::SetFitRange(double x_min, double x_max)
{
    m_x_min = x_min;
    m_x_max = x_max;
}