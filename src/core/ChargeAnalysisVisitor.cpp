#include "ChargeAnalysisVisitor.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "MapObject.hpp"
#include "DataObjectManager.hpp"
#include "MapInterpolationVisitor.hpp"
#include "HRLModelHelper.hpp"
#include "ScopeTimer.hpp"
#include "AtomicChargeEntry.hpp"
#include "GroupChargeEntry.hpp"
#include "AtomicPotentialEntry.hpp"
#include "GroupPotentialEntry.hpp"
#include "AtomicInfoHelper.hpp"
#include "KeyPacker.hpp"
#include "ArrayStats.hpp"
#include "AtomClassifier.hpp"
#include "GlobalEnumClass.hpp"
#include "AminoAcidInfoHelper.hpp"
#include "ElectricPotential.hpp"

#include <iostream>
#include <tuple>
#include <fstream>
#include <unordered_set>

ChargeAnalysisVisitor::ChargeAnalysisVisitor(void) :
    m_thread_size{ 1 },
    m_alpha_r{ 0.0 }, m_alpha_g{ 0.0 }, m_x_min{ 0.0 }, m_x_max{ 0.0 }
{

}

ChargeAnalysisVisitor::~ChargeAnalysisVisitor()
{

}

void ChargeAnalysisVisitor::VisitAtomObject(AtomObject * data_object)
{
    if (data_object == nullptr) return;
    if (data_object->GetAtomicPotentialEntry() == nullptr)
    {
        data_object->SetSelectedFlag(false);
        return;
    }
    data_object->SetSelectedFlag(true);
}

void ChargeAnalysisVisitor::VisitModelObject(ModelObject * data_object)
{
    ScopeTimer timer("ChargeAnalysisVisitor::VisitModelObject");
    if (data_object == nullptr) return;
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
    if (data_object == nullptr) return;
}

void ChargeAnalysisVisitor::Analysis(DataObjectManager * data_manager)
{
    std::cout <<"- Analysis..." << std::endl;
    try
    {
        const auto & model_object{ data_manager->GetDataObjectRef(m_model_key_tag) };
        model_object->Accept(this);

        const auto & neutral_model_object{ data_manager->GetDataObjectRef(m_neutral_model_key_tag) };
        const auto & positive_model_object{ data_manager->GetDataObjectRef(m_positive_model_key_tag) };
        const auto & negative_model_object{ data_manager->GetDataObjectRef(m_negative_model_key_tag) };

        RunChargeFitting(dynamic_cast<ModelObject*>(neutral_model_object.get()),
                         dynamic_cast<ModelObject*>(positive_model_object.get()),
                         dynamic_cast<ModelObject*>(negative_model_object.get()));
    }
    catch(const std::exception & e)
    {
        std::cerr << e.what() << std::endl;
    }
}

void ChargeAnalysisVisitor::RunChargeFitting(
    ModelObject * model_neutral, ModelObject * model_pos, ModelObject * model_neg)
{
    ScopeTimer timer("ChargeAnalysisVisitor::RunChargeFitting");
    auto electric_potential{ std::make_unique<ElectricPotential>() };
    electric_potential->SetModelChoice(1);
    electric_potential->SetBlurringWidth(0.50);

    auto gaus_potential{ std::make_unique<ElectricPotential>() };
    gaus_potential->SetModelChoice(2);

    auto neutral_atom_map{ BuildSerialIDAtomObjectMap(model_neutral) };
    auto positive_atom_map{ BuildSerialIDAtomObjectMap(model_pos) };
    auto negative_atom_map{ BuildSerialIDAtomObjectMap(model_neg) };

    // Atom Classification
    std::unordered_map<int, std::vector<AtomObject *>> atom_list_map;
    for (auto atom : m_selected_atom_list)
    {
        if (atom->GetSpecialAtomFlag() == true) continue;
        if      (atom->GetElement() == Element::CARBON && atom->GetRemoteness() == Remoteness::ALPHA) atom_list_map[0].emplace_back(atom);
        else if (atom->GetElement() == Element::CARBON && atom->GetRemoteness() == Remoteness::NONE) atom_list_map[1].emplace_back(atom);
        else if (atom->GetElement() == Element::NITROGEN && atom->GetRemoteness() == Remoteness::NONE) atom_list_map[2].emplace_back(atom);
        else if (atom->GetElement() == Element::OXYGEN && atom->GetRemoteness() == Remoteness::NONE) atom_list_map[3].emplace_back(atom);
    }
    std::cout << "Size of group :" << atom_list_map.size() << std::endl;

    

    // Group Charge Fitting
    for (const auto & [id, atom_list] : atom_list_map)
    {
        auto group_size{ atom_list.size() };

        std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> data_array;
        std::vector<std::tuple<float, float, float>> regression_data_list;
        std::vector<Eigen::VectorXd> sampling_entry_list; // same as regression_data_list
        data_array.reserve(group_size);
        regression_data_list.reserve(group_size);
        sampling_entry_list.reserve(group_size);
        for (auto atom : atom_list)
        {
            auto potential_entry{ atom->GetAtomicPotentialEntry() };
            if (potential_entry == nullptr)
            {
                std::cout <<"Warning: Atomic potential entry is nullptr for atom "
                        << atom->GetInfo() << ". This atom will be skipped." << std::endl;
                continue;
            }
            auto serial_id{ atom->GetSerialID() };
            auto atom_neutral{ neutral_atom_map.at(serial_id) };
            auto atom_positive{ positive_atom_map.at(serial_id) };
            auto atom_negative{ negative_atom_map.at(serial_id) };

            auto potential_entry_neutral{ atom_neutral->GetAtomicPotentialEntry() };
            auto potential_entry_positive{ atom_positive->GetAtomicPotentialEntry() };
            auto potential_entry_negative{ atom_negative->GetAtomicPotentialEntry() };

            auto charge{ AminoAcidInfoHelper::GetPartialCharge(
                atom->GetResidue(),
                atom->GetElement(),
                atom->GetRemoteness(),
                atom->GetBranch(),
                atom->GetStructure()) };
            auto is_negative_charged{ charge < 0.0 };
            
            auto distance{ 0.5 };
            auto func_phi_0{
                electric_potential->GetPotentialValue(
                    atom->GetElement(), distance, 0.0,
                    potential_entry_neutral->GetAmplitudeEstimateMDPDE(),
                    potential_entry_neutral->GetWidthEstimateMDPDE()) };
            auto func_phi_pos{
                electric_potential->GetPotentialValue(
                    atom->GetElement(), distance, 1.0,
                    potential_entry_positive->GetAmplitudeEstimateMDPDE(),
                    potential_entry_positive->GetWidthEstimateMDPDE()) };
            auto func_phi_neg{
                electric_potential->GetPotentialValue(
                    atom->GetElement(), distance, -1.0,
                    potential_entry_negative->GetAmplitudeEstimateMDPDE(),
                    potential_entry_negative->GetWidthEstimateMDPDE()) };


            auto x0{ func_phi_0 };
            //auto x1_pos{ func_phi_pos - func_phi_0 };
            //auto x1_neg{ func_phi_0 - func_phi_neg };
            auto x1_pos{ func_phi_pos };
            auto x1_neg{ func_phi_neg };
            auto x1{ (is_negative_charged) ? -x1_neg : x1_pos };
            auto amplitude{ potential_entry->GetAmplitudeEstimateMDPDE() };
            auto width{ potential_entry->GetWidthEstimateMDPDE() };
            auto y{ gaus_potential->GetPotentialValue(atom->GetElement(), distance, 0.0, amplitude, width) };

            Eigen::VectorXd sampling_entry(4);
            sampling_entry(0) = 1.0;
            sampling_entry(1) = x0;
            sampling_entry(2) = x1;
            sampling_entry(3) = y;
            sampling_entry_list.emplace_back(sampling_entry);
            regression_data_list.emplace_back(std::make_tuple(x0, x1, y));
            data_array.emplace_back(std::make_tuple(sampling_entry_list, ""));
        }
        

        auto model_estimator{ std::make_unique<HRLModelHelper>(3, static_cast<int>(group_size)) };
        model_estimator->SetDataArray(data_array);
        model_estimator->RunEstimation(m_alpha_r, m_alpha_g);

        auto model_group_mean{ model_estimator->GetMuVectorMean() };
        auto model_group_mdpde{ model_estimator->GetMuVectorMDPDE() };
        auto prior_estimate{ model_estimator->GetMuVectorPrior() };
        auto prior_variance{ model_estimator->GetCapitalLambdaMatrix() };
        std::cout <<"[ID-"<< id <<"] "<< model_group_mdpde.transpose() <<"\t"<< model_group_mdpde(2)/(model_group_mdpde(1)+model_group_mdpde(2)) << std::endl;
        std::cout <<"[ID-"<< id <<"] "<< prior_estimate.transpose() <<"\t"<< prior_estimate(2)/(prior_estimate(1)+prior_estimate(2)) << std::endl;
    }
}

void ChargeAnalysisVisitor::SetFitRange(double x_min, double x_max)
{
    m_x_min = x_min;
    m_x_max = x_max;
}

std::unordered_map<int, AtomObject *>
ChargeAnalysisVisitor::BuildSerialIDAtomObjectMap(ModelObject * model_object)
{
    std::unordered_map<int, AtomObject *> map;
    for (auto & atom : model_object->GetSelectedAtomList())
    {
        if (atom->GetAtomicPotentialEntry() == nullptr) continue;
        auto serial_id{ atom->GetSerialID() };
        if (map.find(serial_id) == map.end())
        {
            map[serial_id] = atom;
        }
        else
        {
            std::cerr << "Warning: Duplicate serial ID found: " << serial_id << std::endl;
        }
    }
    std::cout << "Number of atoms in map: " << map.size() << std::endl;
    return map;
}