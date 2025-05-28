#include "ChargeAnalysisVisitor.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "MapObject.hpp"
#include "DataObjectManager.hpp"
#include "SphereSampler.hpp"
#include "MapInterpolationVisitor.hpp"
#include "HRLModelHelper.hpp"
#include "ScopeTimer.hpp"
#include "AtomicChargeEntry.hpp"
#include "GroupChargeEntry.hpp"
#include "AtomicInfoHelper.hpp"
#include "KeyPacker.hpp"
#include "ArrayStats.hpp"
#include "AtomClassifier.hpp"

#include <iostream>
#include <tuple>
#include <fstream>
#include <unordered_set>

ChargeAnalysisVisitor::ChargeAnalysisVisitor(
    SphereSampler * sphere_sampler) :
    m_thread_size{ 1 },
    m_alpha_r{ 0.0 }, m_alpha_g{ 0.0 }, m_x_min{ 0.0 }, m_x_max{ 0.0 },
    m_sphere_sampler{ sphere_sampler }
{

}

ChargeAnalysisVisitor::~ChargeAnalysisVisitor()
{

}

void ChargeAnalysisVisitor::VisitAtomObject(AtomObject * data_object)
{
    if (data_object == nullptr) return;
    data_object->SetSelectedFlag(true);
}

void ChargeAnalysisVisitor::VisitModelObject(ModelObject * data_object)
{
    ScopeTimer timer("ChargeAnalysisVisitor::VisitModelObject");
    if (data_object == nullptr) return;
    data_object->FilterAtomFromSymmetry(m_is_asymmetry);
    data_object->Update();
    m_selected_atom_list = data_object->GetSelectedAtomList();
    for (auto & atom : m_selected_atom_list)
    {
        auto atom_charge_entry{ std::make_unique<AtomicChargeEntry>() };
        atom->AddAtomicChargeEntry(std::move(atom_charge_entry));
    }
    std::cout <<" Number of selected atom = "<< m_selected_atom_list.size() << std::endl;
}

void ChargeAnalysisVisitor::VisitMapObject(MapObject * data_object)
{
    ScopeTimer timer("ChargeAnalysisVisitor::VisitMapObject");
    if (data_object == nullptr) return;
}

void ChargeAnalysisVisitor::Analysis(DataObjectManager * data_manager)
{
    std::cout <<"- Analysis..." << std::endl;
    try
    {
        const auto & model_object{ data_manager->GetDataObjectRef(m_model_key_tag) };
        model_object->Accept(this);

        const auto & map_object{ data_manager->GetDataObjectRef(m_map_key_tag) };
        const auto & neutral_map_object{ data_manager->GetDataObjectRef(m_neutral_map_key_tag) };
        const auto & positive_map_object{ data_manager->GetDataObjectRef(m_positive_map_key_tag) };
        const auto & negative_map_object{ data_manager->GetDataObjectRef(m_negative_map_key_tag) };
        
        for (auto & atom : m_selected_atom_list)
        {
            auto sampling_points{ m_sphere_sampler->GenerateSamplingPoints(atom->GetPosition()) };
            MapInterpolationVisitor interpolation_visitor{ sampling_points };
            auto entry{ atom->GetAtomicChargeEntry() };
            map_object->Accept(&interpolation_visitor);
            entry->AddDistanceAndMapValueList(interpolation_visitor.GetSamplingDataList());
            neutral_map_object->Accept(&interpolation_visitor);
            entry->AddDistanceAndNeutralMapValueList(interpolation_visitor.GetSamplingDataList());
            positive_map_object->Accept(&interpolation_visitor);
            entry->AddDistanceAndPositiveMapValueList(interpolation_visitor.GetSamplingDataList());
            negative_map_object->Accept(&interpolation_visitor);
            entry->AddDistanceAndNegativeMapValueList(interpolation_visitor.GetSamplingDataList());
        }

        RunChargeFitting(dynamic_cast<ModelObject*>(model_object.get()));
    }
    catch(const std::exception & e)
    {
        std::cerr << e.what() << std::endl;
    }
}

void ChargeAnalysisVisitor::RunChargeFitting(ModelObject * model_object)
{
    ScopeTimer timer("ChargeAnalysisVisitor::RunChargeFitting");
    for (size_t i = 0; i < AtomicInfoHelper::GetGroupClassCount(); i++)
    {
        const auto & class_key{ AtomicInfoHelper::GetGroupClassKey(i) };

        // Atom Classification
        std::unordered_set<uint64_t> group_key_set;
        auto group_charge_entry( std::make_unique<GroupChargeEntry>() );
        for (auto atom : m_selected_atom_list)
        {
            auto group_key{ AtomClassifier::GetGroupKeyInClass(atom, class_key) };
            group_charge_entry->AddAtomObjectPtr(group_key, atom);
            group_charge_entry->InsertGroupKey(group_key);
            group_key_set.insert(group_key);
        }

        // Group Charge Fitting
        for (const auto & group_key : group_key_set)
        {
            auto atom_list{ group_charge_entry->GetAtomObjectPtrList(group_key) };
            auto group_size{ atom_list.size() };
            std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> data_array;
            data_array.reserve(group_size);
            for (auto atom : atom_list)
            {
                auto entry{ atom->GetAtomicChargeEntry() };
                if (entry == nullptr) continue;
                if (entry->IsDataListSizeValid() == false) continue;
                bool is_negative_charged{ false };
                std::vector<Eigen::VectorXd> sampling_entry_list;
                sampling_entry_list.reserve(entry->GetDistanceAndMapValueListSize());
                for (size_t p = 0; p < entry->GetDistanceAndMapValueListSize(); p++)
                {
                    auto data_entry{ entry->GetDistanceAndMapValueList().at(p) };
                    auto neutral_data_entry{ entry->GetDistanceAndNeutralMapValueList().at(p) };
                    auto positive_data_entry{ entry->GetDistanceAndPositiveMapValueList().at(p) };
                    auto negative_data_entry{ entry->GetDistanceAndNegativeMapValueList().at(p) };
                    auto distance{ std::get<0>(data_entry) };
                    if (distance < m_x_min || distance > m_x_max) continue;

                    auto y_data{ std::get<1>(data_entry) };
                    auto y_neutral{ std::get<1>(neutral_data_entry) };
                    auto y_positive{ std::get<1>(positive_data_entry) };
                    auto y_negative{ std::get<1>(negative_data_entry) };
                    auto model_x_positive{  y_positive - y_neutral };
                    auto model_x_negative{  y_neutral - y_negative };
                    auto model_y{ y_data - y_neutral };
                    
                    Eigen::VectorXd sampling_entry(2);
                    sampling_entry(0) = (is_negative_charged) ? model_x_negative : model_x_positive;
                    sampling_entry(1) = model_y;
                    sampling_entry_list.emplace_back(sampling_entry);
                }
                data_array.emplace_back(std::make_tuple(sampling_entry_list, atom->GetInfo()));
            }
            auto model_estimator{ std::make_unique<HRLModelHelper>(2, static_cast<int>(group_size)) };
            model_estimator->SetDataArray(data_array);
            model_estimator->RunEstimation(m_alpha_r, m_alpha_g);

            auto model_group_mean{ model_estimator->GetMuVectorMean() };
            group_charge_entry->AddModelEstimateMean(
                group_key, model_group_mean(0), model_group_mean(1)
            );

            auto model_group_mdpde{ model_estimator->GetMuVectorMDPDE() };
            group_charge_entry->AddModelEstimateMDPDE(
                group_key, model_group_mdpde(0), model_group_mdpde(1)
            );

            auto prior_estimate{ model_estimator->GetMuVectorPrior() };
            auto prior_variance{ model_estimator->GetCapitalLambdaMatrix() };
            group_charge_entry->AddModelEstimatePrior(group_key, prior_estimate(0), prior_estimate(1));
            group_charge_entry->AddModelVariancePrior(group_key, prior_variance(0, 0), prior_variance(1, 1));

            auto count{ 0 };
            for (auto atom : atom_list)
            {
                auto atom_entry{ atom->GetAtomicChargeEntry() };
                Eigen::VectorXd beta_vector_ols{ model_estimator->GetBetaMatrixOLS(count) };
                atom_entry->AddModelEstimateOLS(beta_vector_ols(0), beta_vector_ols(1));

                Eigen::VectorXd beta_vector_mdpde{ model_estimator->GetBetaMatrixMDPDE(count) };
                atom_entry->AddModelEstimateMDPDE(beta_vector_mdpde(0), beta_vector_mdpde(1));

                Eigen::VectorXd beta_vector_posterior{ model_estimator->GetBetaMatrixPosterior(count) };
                auto sigma_matrix_posterior{ model_estimator->GetCapitalSigmaMatrixPosterior(count) };
                atom_entry->AddModelEstimatePosterior(class_key, beta_vector_posterior(0), beta_vector_posterior(1));
                atom_entry->AddModelVariancePosterior(class_key, sigma_matrix_posterior(0, 0), sigma_matrix_posterior(1, 1));
                atom_entry->AddOutlierTag(class_key, model_estimator->GetOutlierFlag(count));
                atom_entry->AddStatisticalDistance(class_key, model_estimator->GetStatisticalDistance(count));
                count++;
            }
        }
        model_object->AddGroupChargeEntry(class_key, group_charge_entry);
    }
}

void ChargeAnalysisVisitor::SetAsymmetryFlag(bool value)
{
    m_is_asymmetry = value;
}

void ChargeAnalysisVisitor::SetThreadSize(unsigned int thread_size)
{
    m_thread_size = thread_size;
}

void ChargeAnalysisVisitor::SetFitRange(double x_min, double x_max)
{
    m_x_min = x_min;
    m_x_max = x_max;
}

void ChargeAnalysisVisitor::SetAlphaR(double alpha_r)
{
    m_alpha_r = alpha_r;
}

void ChargeAnalysisVisitor::SetAlphaG(double alpha_g)
{
    m_alpha_g = alpha_g;
}

void ChargeAnalysisVisitor::SetMapObjectKeyTag(const std::string & value)
{
    m_map_key_tag = value;
}

void ChargeAnalysisVisitor::SetModelObjectKeyTag(const std::string & value)
{
    m_model_key_tag = value;
}

void ChargeAnalysisVisitor::SetNeutralMapObjectKeyTag(const std::string & value)
{
    m_neutral_map_key_tag = value;
}

void ChargeAnalysisVisitor::SetPositiveMapObjectKeyTag(const std::string & value)
{
    m_positive_map_key_tag = value;
}

void ChargeAnalysisVisitor::SetNegativeMapObjectKeyTag(const std::string & value)
{
    m_negative_map_key_tag = value;
}