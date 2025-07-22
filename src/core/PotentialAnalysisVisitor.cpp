#include "PotentialAnalysisVisitor.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "MapObject.hpp"
#include "DataObjectManager.hpp"
#include "SphereSampler.hpp"
#include "MapInterpolationVisitor.hpp"
#include "HRLModelHelper.hpp"
#include "ScopeTimer.hpp"
#include "AtomicPotentialEntry.hpp"
#include "GroupPotentialEntry.hpp"
#include "AtomicInfoHelper.hpp"
#include "KeyPacker.hpp"
#include "ArrayStats.hpp"
#include "AtomClassifier.hpp"
#include "GausLinearTransformHelper.hpp"
#include "Logger.hpp"

#include <tuple>
#include <fstream>
#include <unordered_set>

PotentialAnalysisVisitor::PotentialAnalysisVisitor(
    SphereSampler * sphere_sampler,
    const PotentialAnalysisCommand::Options & options) :
    m_is_asymmetry{ options.is_asymmetry },
    m_thread_size{ static_cast<unsigned int>(options.thread_size) },
    m_alpha_r{ options.alpha_r }, m_alpha_g{ options.alpha_g },
    m_x_min{ options.fit_range_min }, m_x_max{ options.fit_range_max },
    m_map_key_tag{ "map" }, m_model_key_tag{ "model" },
    m_sphere_sampler{ sphere_sampler }
{

}

void PotentialAnalysisVisitor::VisitAtomObject(AtomObject * data_object)
{
    if (data_object == nullptr) return;
    data_object->SetSelectedFlag(true);
}

void PotentialAnalysisVisitor::VisitModelObject(ModelObject * data_object)
{
    ScopeTimer timer("PotentialAnalysisVisitor::VisitModelObject");
    if (data_object == nullptr) return;
    data_object->FilterAtomFromSymmetry(m_is_asymmetry);
    data_object->Update();
    m_selected_atom_list = data_object->GetSelectedAtomList();
    for (auto & atom : m_selected_atom_list)
    {
        auto atom_potential_entry{ std::make_unique<AtomicPotentialEntry>() };
        atom->AddAtomicPotentialEntry(std::move(atom_potential_entry));
    }
    Logger::Log(LogLevel::Info, " Number of selected atom = "
                + std::to_string(m_selected_atom_list.size()));
}

void PotentialAnalysisVisitor::VisitMapObject(MapObject * data_object)
{
    ScopeTimer timer("PotentialAnalysisVisitor::VisitMapObject");
    if (data_object == nullptr) return;
    data_object->MapValueArrayNormalization();
    m_sphere_sampler->Print();
    MapInterpolationVisitor interpolation_visitor{ m_sphere_sampler };
    for (auto & atom : m_selected_atom_list)
    {
        auto entry{ atom->GetAtomicPotentialEntry() };
        interpolation_visitor.SetPosition(atom->GetPosition());
        data_object->Accept(&interpolation_visitor);
        entry->AddDistanceAndMapValueList(interpolation_visitor.GetSamplingDataList());
    }
}

void PotentialAnalysisVisitor::VisitDataObjectManager(DataObjectManager * data_manager)
{
    Logger::Log(LogLevel::Debug, "PotentialAnalysisVisitor::VisitDataObjectManager() called");
    try
    {
        auto model_object{ data_manager->GetTypedDataObjectPtr<ModelObject>(m_model_key_tag) };
        model_object->Accept(this);

        auto map_object{ data_manager->GetTypedDataObjectPtr<MapObject>(m_map_key_tag) };
        map_object->Accept(this);

        RunPotentialFitting(model_object);
    }
    catch(const std::exception & e)
    {
        Logger::Log(LogLevel::Error, e.what());
    }
}

void PotentialAnalysisVisitor::RunPotentialFitting(ModelObject * model_object)
{
    ScopeTimer timer("PotentialAnalysisVisitor::RunPotentialFitting");
    for (size_t i = 0; i < AtomicInfoHelper::GetGroupClassCount(); i++)
    {
        const auto & class_key{ AtomicInfoHelper::GetGroupClassKey(i) };

        // Atom Classification
        std::unordered_set<uint64_t> group_key_set;
        auto group_potential_entry( std::make_unique<GroupPotentialEntry>() );
        for (auto atom : m_selected_atom_list)
        {
            auto group_key{ AtomClassifier::GetGroupKeyInClass(atom, class_key) };
            group_potential_entry->AddAtomObjectPtr(group_key, atom);
            group_potential_entry->InsertGroupKey(group_key);
            group_key_set.insert(group_key);
        }

        // Group Potential Fitting
        for (const auto & group_key : group_key_set)
        {
            auto atom_list{ group_potential_entry->GetAtomObjectPtrList(group_key) };
            auto group_size{ atom_list.size() };
            std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> data_array;
            data_array.reserve(group_size);
            for (auto atom : atom_list)
            {
                auto entry{ atom->GetAtomicPotentialEntry() };
                std::vector<Eigen::VectorXd> sampling_entry_list;
                sampling_entry_list.reserve(static_cast<size_t>(entry->GetDistanceAndMapValueListSize()));
                for (auto & data_entry : entry->GetDistanceAndMapValueList())
                {
                    auto gaus_x{ static_cast<double>(std::get<0>(data_entry)) };
                    auto gaus_y{ static_cast<double>(std::get<1>(data_entry)) };
                    if (gaus_x < m_x_min || gaus_x > m_x_max) continue;
                    if (gaus_y <= 0.0) continue;
                    sampling_entry_list.emplace_back(
                        GausLinearTransformHelper::BuildLinearModelDataVector(gaus_x, gaus_y)
                    );
                }
                data_array.emplace_back(std::make_tuple(sampling_entry_list, atom->GetInfo()));
            }
            auto model_estimator{ std::make_unique<HRLModelHelper>(2, static_cast<int>(group_size)) };
            model_estimator->SetDataArray(data_array);
            model_estimator->RunEstimation(m_alpha_r, m_alpha_g);

            auto gaus_group_mean{
                GausLinearTransformHelper::BuildGausModel(model_estimator->GetMuVectorMean())
            };
            group_potential_entry->AddGausEstimateMean(
                group_key, gaus_group_mean(0), gaus_group_mean(1)
            );

            auto gaus_group_mdpde{
                GausLinearTransformHelper::BuildGausModel(model_estimator->GetMuVectorMDPDE())
            };
            group_potential_entry->AddGausEstimateMDPDE(
                group_key, gaus_group_mdpde(0), gaus_group_mdpde(1)
            );

            auto gaus_prior{
                GausLinearTransformHelper::BuildGausModelWithVariance(
                    model_estimator->GetMuVectorPrior(), model_estimator->GetCapitalLambdaMatrix())
            };
            auto prior_estimate{ std::get<0>(gaus_prior) };
            auto prior_variance{ std::get<1>(gaus_prior) };
            group_potential_entry->AddGausEstimatePrior(group_key, prior_estimate(0), prior_estimate(1));
            group_potential_entry->AddGausVariancePrior(group_key, prior_variance(0), prior_variance(1));

            auto count{ 0 };
            for (auto atom : atom_list)
            {
                auto atom_entry{ atom->GetAtomicPotentialEntry() };
                Eigen::VectorXd beta_vector_ols{ model_estimator->GetBetaMatrixOLS(count) };
                auto gaus_ols{ GausLinearTransformHelper::BuildGausModel(beta_vector_ols) };
                atom_entry->AddGausEstimateOLS(gaus_ols(0), gaus_ols(1));

                Eigen::VectorXd beta_vector_mdpde{ model_estimator->GetBetaMatrixMDPDE(count) };
                auto gaus_mdpde{ GausLinearTransformHelper::BuildGausModel(beta_vector_mdpde) };
                atom_entry->AddGausEstimateMDPDE(gaus_mdpde(0), gaus_mdpde(1));

                Eigen::VectorXd beta_vector_posterior{ model_estimator->GetBetaMatrixPosterior(count) };
                auto sigma_matrix_posterior{ model_estimator->GetCapitalSigmaMatrixPosterior(count) };
                auto gaus_posterior{
                    GausLinearTransformHelper::BuildGausModelWithVariance(
                        beta_vector_posterior, sigma_matrix_posterior)
                };
                auto posterior_estimate{ std::get<0>(gaus_posterior) };
                auto posterior_variance{ std::get<1>(gaus_posterior) };
                atom_entry->AddGausEstimatePosterior(class_key, posterior_estimate(0), posterior_estimate(1));
                atom_entry->AddGausVariancePosterior(class_key, posterior_variance(0), posterior_variance(1));
                atom_entry->AddOutlierTag(class_key, model_estimator->GetOutlierFlag(count));
                atom_entry->AddStatisticalDistance(class_key, model_estimator->GetStatisticalDistance(count));
                count++;
            }
        }
        model_object->AddGroupPotentialEntry(class_key, group_potential_entry);
    }
}
