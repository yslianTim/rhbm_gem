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
#include "GroupPotentialEntry.hpp"
#include "AtomicInfoHelper.hpp"
#include "KeyPacker.hpp"
#include "ArrayStats.hpp"

#include <iostream>
#include <tuple>
#include <fstream>

PotentialAnalysisVisitor::PotentialAnalysisVisitor(
    AtomSelector * atom_selector,
    SphereSampler * sphere_sampler) :
    m_thread_size{ 1 },
    m_alpha_r{ 0.0 }, m_alpha_g{ 0.0 },
    m_x_min{ 0.0 }, m_x_max{ 0.0 },
    m_atom_selector{ atom_selector },
    m_sphere_sampler{ sphere_sampler }
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
    auto selected_atom_size{ static_cast<size_t>(data_object->GetNumberOfSelectedAtom()) };
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
    MapInterpolationVisitor interpolation_visitor{ m_sphere_sampler };
    for (auto & atom : m_selected_atom_list)
    {
        auto entry{ atom->GetAtomicPotentialEntry() };
        interpolation_visitor.SetPosition(atom->GetPosition());
        data_object->Accept(&interpolation_visitor);
        entry->AddDistanceAndMapValueList(interpolation_visitor.GetSamplingDataList());
    }
}

void PotentialAnalysisVisitor::Analysis(DataObjectManager * data_manager)
{
    std::cout <<"- Analysis..." << std::endl;
    try
    {
        const auto & model_object{ data_manager->GetDataObjectRef(m_model_key_tag) };
        model_object->Accept(this);

        const auto & map_object{ data_manager->GetDataObjectRef(m_map_key_tag) };
        map_object->Accept(this);

        //RunAtomPositionDumping(); // For test, to be move to other place
        //RunMapValueDumping(dynamic_cast<MapObject *>(map_object.get())); // For test, to be move to other place

        for (size_t i = 0; i < AtomicInfoHelper::GetGroupClassCount(); i++)
        {
            auto group_class_key{ AtomicInfoHelper::GetGroupClassKey(i) };
            RunAtomClassification(group_class_key, dynamic_cast<ModelObject*>(model_object.get()));
            RunPotentialFitting(group_class_key, dynamic_cast<ModelObject*>(model_object.get()));
        }
    }
    catch(const std::exception & e)
    {
        std::cerr << e.what() << std::endl;
    }
}

void PotentialAnalysisVisitor::RunAtomClassification(
    const std::string & class_key, ModelObject * model_object)
{
    ScopeTimer timer("PotentialAnalysisVisitor::RunAtomClassification");
    std::cout <<"- RunAtomClassification..." << std::endl;
    auto group_potential_entry( std::make_unique<GroupPotentialEntry>() );
    for (auto atom : m_selected_atom_list)
    {
        uint64_t group_key;
        if (class_key == AtomicInfoHelper::GetResidueClassKey())
        {
            group_key = KeyPackerResidueClass::Pack(
                atom->GetResidue(), atom->GetElement(),
                atom->GetRemoteness(), atom->GetBranch(),
                atom->GetSpecialAtomFlag());
        }
        else if (class_key == AtomicInfoHelper::GetElementClassKey())
        {
            group_key = KeyPackerElementClass::Pack(
                atom->GetElement(), atom->GetRemoteness(),
                atom->GetSpecialAtomFlag());
        }
        else if (class_key == AtomicInfoHelper::GetStructureClassKey())
        {
            group_key = KeyPackerStructureClass::Pack(
                atom->GetStructure(), atom->GetResidue(), atom->GetElement(),
                atom->GetRemoteness(), atom->GetBranch(),
                atom->GetSpecialAtomFlag());
        }
        else
        {
            throw std::runtime_error("PotentialAnalysisVisitor::RunAtomClassification()"
                                     " : Unsupported class key."+ class_key);
        }
        group_potential_entry->AddAtomObjectPtr(group_key, atom);
        group_potential_entry->InsertGroupKey(group_key);
        m_group_set_map[class_key].insert(group_key);
    }
    model_object->AddGroupPotentialEntry(class_key, group_potential_entry);
}

void PotentialAnalysisVisitor::RunPotentialFitting(
    const std::string & class_key, ModelObject * model_object)
{
    ScopeTimer timer("PotentialAnalysisVisitor::RunPotentialFitting");
    std::cout <<"- RunPotentialFitting..." << std::endl;
    auto group_potential_entry{ model_object->GetGroupPotentialEntry(class_key) };
    
    for (const auto & group_key : m_group_set_map.at(class_key))
    {
        auto atom_list{ group_potential_entry->GetAtomObjectPtrList(group_key) };
        auto group_size{ static_cast<int>(atom_list.size()) };
        std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> data_array;
        data_array.reserve(static_cast<size_t>(group_size));
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
                sampling_entry_list.emplace_back(GausLinearTransformHelper::BuildLinearModelDataVector(gaus_x, gaus_y));
            }
            data_array.emplace_back(std::make_tuple(sampling_entry_list, atom->GetInfo()));
        }
        auto model_estimator{ std::make_unique<HRLModelHelper>(2, group_size) };
        model_estimator->SetDataArray(data_array);
        model_estimator->RunEstimation(m_alpha_r, m_alpha_g);

        auto gaus_group_mean{ GausLinearTransformHelper::BuildGausModel(model_estimator->GetMuVectorMean()) };
        group_potential_entry->AddGausEstimateMean(group_key, gaus_group_mean(0), gaus_group_mean(1));

        auto gaus_group_mdpde{ GausLinearTransformHelper::BuildGausModel(model_estimator->GetMuVectorMDPDE()) };
        group_potential_entry->AddGausEstimateMDPDE(group_key, gaus_group_mdpde(0), gaus_group_mdpde(1));

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
            auto gaus_posterior{ GausLinearTransformHelper::BuildGausModelWithVariance(beta_vector_posterior, sigma_matrix_posterior) };
            auto posterior_estimate{ std::get<0>(gaus_posterior) };
            auto posterior_variance{ std::get<1>(gaus_posterior) };
            atom_entry->AddGausEstimatePosterior(class_key, posterior_estimate(0), posterior_estimate(1));
            atom_entry->AddGausVariancePosterior(class_key, posterior_variance(0), posterior_variance(1));
            atom_entry->AddOutlierTag(class_key, model_estimator->GetOutlierFlag(count));
            atom_entry->AddStatisticalDistance(class_key, model_estimator->GetStatisticalDistance(count));
            count++;
        }
    }
}

void PotentialAnalysisVisitor::RunMapValueDumping(MapObject * map_object)
{
    ScopeTimer timer("PotentialAnalysisVisitor::RunMapValueDumping");
    std::cout <<"- RunMapValueDumping..." << std::endl;
    auto selected_atom_size{ m_selected_atom_list.size() };
    std::array<float, 3> atom_range_min, atom_range_max;
    std::vector<float> x_list, y_list, z_list;
    x_list.reserve(selected_atom_size);
    y_list.reserve(selected_atom_size);
    z_list.reserve(selected_atom_size);
    for (auto & atom : m_selected_atom_list)
    {
        x_list.emplace_back(atom->GetPosition().at(0));
        y_list.emplace_back(atom->GetPosition().at(1));
        z_list.emplace_back(atom->GetPosition().at(2));
    }
    atom_range_min.at(0) = ArrayStats<float>::ComputeMin(x_list.data(), selected_atom_size) - 3.0f;
    atom_range_min.at(1) = ArrayStats<float>::ComputeMin(y_list.data(), selected_atom_size) - 3.0f;
    atom_range_min.at(2) = ArrayStats<float>::ComputeMin(z_list.data(), selected_atom_size) - 3.0f;
    atom_range_max.at(0) = ArrayStats<float>::ComputeMax(x_list.data(), selected_atom_size) + 3.0f;
    atom_range_max.at(1) = ArrayStats<float>::ComputeMax(y_list.data(), selected_atom_size) + 3.0f;
    atom_range_max.at(2) = ArrayStats<float>::ComputeMax(z_list.data(), selected_atom_size) + 3.0f;

    std::string output_path{ "/Users/yslian/Downloads/map_value_list_emd_11103_unsharpened.csv" };
    std::ofstream outfile(output_path);
    if (!outfile.is_open())
    {
        std::cerr << "Error: Could not open file " << output_path << " for writing.\n";
        return;
    }
    outfile << "GridID,X,Y,Z,MapValue\n";
    auto count{ 0 };
    for (size_t i = 0; i < map_object->GetMapValueArraySize(); i++)
    {
        auto grid_position{ map_object->GetGridPosition(i) };
        if (grid_position.at(0) < atom_range_min.at(0) || grid_position.at(0) > atom_range_max.at(0)) continue;
        if (grid_position.at(1) < atom_range_min.at(1) || grid_position.at(1) > atom_range_max.at(1)) continue;
        if (grid_position.at(2) < atom_range_min.at(2) || grid_position.at(2) > atom_range_max.at(2)) continue;
        auto map_value{ map_object->GetMapValue(i) };
        if (map_value <= 0.0f) continue;
        outfile << i <<','
                << grid_position.at(0) <<','<< grid_position.at(1) <<','<< grid_position.at(2) <<','
                << map_value <<'\n';
        count++;
    }
    std::cout <<"Selected map grid size = "<< count <<" / "<< map_object->GetMapValueArraySize() << std::endl;
}

void PotentialAnalysisVisitor::RunAtomPositionDumping(void)
{
    ScopeTimer timer("PotentialAnalysisVisitor::RunAtomPositionDumping");
    std::cout <<"- RunAtomPositionDumping..." << std::endl;

    std::string output_path{ "/Users/yslian/Downloads/atom_position_list_6z6u.csv" };
    std::ofstream outfile(output_path);
    if (!outfile.is_open())
    {
        std::cerr << "Error: Could not open file " << output_path << " for writing.\n";
        return;
    }
    outfile << "SerialID,X,Y,Z\n";
    for (auto & atom : m_selected_atom_list)
    {
        outfile << atom->GetSerialID() <<','
                << atom->GetPosition().at(0) <<','
                << atom->GetPosition().at(1) <<','
                << atom->GetPosition().at(2) <<'\n';
    }
    outfile.close();
    std::cout <<"Number of selected atom to be output = "<< m_selected_atom_list.size() << std::endl;
}

void PotentialAnalysisVisitor::SetThreadSize(unsigned int thread_size)
{
    m_thread_size = thread_size;
}

void PotentialAnalysisVisitor::SetFitRange(double x_min, double x_max)
{
    m_x_min = x_min;
    m_x_max = x_max;
}

void PotentialAnalysisVisitor::SetAlphaR(double alpha_r)
{
    m_alpha_r = alpha_r;
}

void PotentialAnalysisVisitor::SetAlphaG(double alpha_g)
{
    m_alpha_g = alpha_g;
}

void PotentialAnalysisVisitor::SetMapObjectKeyTag(const std::string & value)
{
    m_map_key_tag = value;
}

void PotentialAnalysisVisitor::SetModelObjectKeyTag(const std::string & value)
{
    m_model_key_tag = value;
}