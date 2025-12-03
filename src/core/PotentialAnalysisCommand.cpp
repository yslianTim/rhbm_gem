#include "PotentialAnalysisCommand.hpp"
#include "DataObjectManager.hpp"
#include "AtomObject.hpp"
#include "BondObject.hpp"
#include "MapObject.hpp"
#include "ModelObject.hpp"
#include "MapInterpolationVisitor.hpp"
#include "HRLModelHelper.hpp"
#include "ScopeTimer.hpp"
#include "FilePathHelper.hpp"
#include "LocalPotentialEntry.hpp"
#include "GroupPotentialEntry.hpp"
#include "ChemicalDataHelper.hpp"
#include "AtomClassifier.hpp"
#include "BondClassifier.hpp"
#include "GausLinearTransformHelper.hpp"
#include "SphereSampler.hpp"
#include "CylinderSampler.hpp"
#include "Logger.hpp"
#include "CommandRegistry.hpp"

#include <unordered_set>
#include <tuple>
#include <vector>
#include <atomic>
#include <utility>
#include <iterator>
#include <stdexcept>
#include <limits>
#include <random>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace {
CommandRegistrar<PotentialAnalysisCommand> registrar_potential_analysis{
    "potential_analysis",
    "Run potential analysis"};
}

PotentialAnalysisCommand::PotentialAnalysisCommand(void) :
    CommandBase(), m_options{}, m_model_key_tag{"model"}, m_map_key_tag{"map"},
    m_map_object{ nullptr }, m_model_object{ nullptr }
{
    Logger::Log(LogLevel::Debug, "PotentialAnalysisCommand::PotentialAnalysisCommand() called.");
}

PotentialAnalysisCommand::~PotentialAnalysisCommand()
{
    Logger::Log(LogLevel::Debug, "PotentialAnalysisCommand::~PotentialAnalysisCommand() called.");
    m_map_object.reset();
    m_model_object.reset();
}

void PotentialAnalysisCommand::RegisterCLIOptionsExtend(CLI::App * cmd)
{
    Logger::Log(LogLevel::Debug, "PotentialAnalysisCommand::RegisterCLIOptionsExtend() called.");
    cmd->add_option_function<std::string>("-a,--model",
        [&](const std::string & value) { SetModelFilePath(value); },
        "Model file path")->required();
    cmd->add_option_function<std::string>("-m,--map",
        [&](const std::string & value) { SetMapFilePath(value); },
        "Map file path")->required();
    cmd->add_option_function<bool>("--simulation",
        [&](bool value) { SetSimulationFlag(value); },
        "Simulation flag")->default_val(m_options.is_simulation);
    cmd->add_option_function<double>("-r,--sim-resolution",
        [&](double value) { SetSimulatedMapResolution(value); },
        "Set simulated map's resolution (blurring width)")
        ->default_val(m_options.resolution_simulation);
    cmd->add_option_function<std::string>("-k,--save-key",
        [&](const std::string & value) { SetSavedKeyTag(value); },
        "New key tag for saving ModelObject results into database")
        ->default_val(m_options.saved_key_tag);
    cmd->add_option_function<bool>("--training-alpha",
        [&](bool value) { SetTrainingAlphaFlag(value); },
        "Turn On/Off alpha training flag")->default_val(m_options.use_training_alpha);
    cmd->add_option_function<bool>("--asymmetry",
        [&](bool value) { SetAsymmetryFlag(value); },
        "Turn On/Off asymmetry flag")->default_val(m_options.is_asymmetry);
    cmd->add_option_function<int>("-s,--sampling",
        [&](int value) { SetSamplingSize(value); },
        "Number of sampling points per atom")->default_val(m_options.sampling_size);
    cmd->add_option_function<double>("--sampling-min",
        [&](double value) { SetSamplingRangeMinimum(value); },
        "Minimum sampling range")->default_val(m_options.sampling_range_min);
    cmd->add_option_function<double>("--sampling-max",
        [&](double value) { SetSamplingRangeMaximum(value); },
        "Maximum sampling range")->default_val(m_options.sampling_range_max);
    cmd->add_option_function<double>("--sampling-height",
        [&](double value) { SetSamplingHeight(value); },
        "Maximum sampling height")->default_val(m_options.sampling_height);
    cmd->add_option_function<double>("--fit-min",
        [&](double value) { SetFitRangeMinimum(value); },
        "Minimum fitting range")->default_val(m_options.fit_range_min);
    cmd->add_option_function<double>("--fit-max",
        [&](double value) { SetFitRangeMaximum(value); },
        "Maximum fitting range")->default_val(m_options.fit_range_max);
    cmd->add_option_function<double>("--alpha-r",
        [&](double value) { SetAlphaR(value); },
        "Alpha value for R")->default_val(m_options.alpha_r);
    cmd->add_option_function<double>("--alpha-g",
        [&](double value) { SetAlphaG(value); },
        "Alpha value for G")->default_val(m_options.alpha_g);
}

bool PotentialAnalysisCommand::Execute(void)
{
    Logger::Log(LogLevel::Debug, "PotentialAnalysisCommand::Execute() called.");
    if (BuildDataObject() == false) return false;
    RunMapObjectPreprocessing();
    RunModelObjectPreprocessing();

    RunAtomMapValueSampling();
    RunAtomGroupClassification();
    if (m_options.use_training_alpha) RunAtomAlphaTraining();
    else RunLocalAtomFitting();
    RunAtomPotentialFitting();

    //RunBondMapValueSampling();
    //RunBondGroupClassification();
    //RunLocalBondFitting();
    //RunBondPotentialFitting();
    SaveDataObject();
    return true;
}

void PotentialAnalysisCommand::SetTrainingAlphaFlag(bool value)
{
    m_options.use_training_alpha = value;
}

void PotentialAnalysisCommand::SetAsymmetryFlag(bool value)
{
    m_options.is_asymmetry = value; 
}
void PotentialAnalysisCommand::SetSimulationFlag(bool value)
{
    m_options.is_simulation = value;
}

void PotentialAnalysisCommand::SetSimulatedMapResolution(double value)
{
    m_options.resolution_simulation = value;
}

void PotentialAnalysisCommand::SetFitRangeMinimum(double value)
{
    m_options.fit_range_min = value;
}

void PotentialAnalysisCommand::SetFitRangeMaximum(double value)
{
    m_options.fit_range_max = value;
}

void PotentialAnalysisCommand::SetAlphaR(double value)
{
    m_options.alpha_r = value;
}

void PotentialAnalysisCommand::SetAlphaG(double value)
{
    m_options.alpha_g = value;
}

void PotentialAnalysisCommand::SetModelFilePath(const std::filesystem::path & path)
{
    m_options.model_file_path = path;
    if (!FilePathHelper::EnsureFileExists(m_options.model_file_path, "Model file"))
    {
        Logger::Log(LogLevel::Error,
            "Model file does not exist: " + m_options.model_file_path.string());
        m_valiate_options = false;
    }
}

void PotentialAnalysisCommand::SetMapFilePath(const std::filesystem::path & path)
{
    m_options.map_file_path = path;
    if (!FilePathHelper::EnsureFileExists(m_options.map_file_path, "Map file"))
    {
        Logger::Log(LogLevel::Error,
            "Map file does not exist: " + m_options.map_file_path.string());
        m_valiate_options = false;
    }
}

void PotentialAnalysisCommand::SetSavedKeyTag(const std::string & tag)
{
    m_options.saved_key_tag = tag;
}

void PotentialAnalysisCommand::SetSamplingSize(int value)
{
    m_options.sampling_size = value;
    if (m_options.sampling_size <= 0)
    {
        Logger::Log(LogLevel::Warning,
            "Sampling size must be positive, reset to default value = 1500");
        m_options.sampling_size = 1500;
    }
}

void PotentialAnalysisCommand::SetSamplingRangeMinimum(double value)
{
    m_options.sampling_range_min = value;
}

void PotentialAnalysisCommand::SetSamplingRangeMaximum(double value)
{
    m_options.sampling_range_max = value;
}

void PotentialAnalysisCommand::SetSamplingHeight(double value)
{
    m_options.sampling_height = value;
}

bool PotentialAnalysisCommand::BuildDataObject(void)
{
    Logger::Log(LogLevel::Debug, "PotentialAnalysisCommand::BuildDataObject() called");
    ScopeTimer timer("PotentialAnalysisCommand::BuildDataObject");
    auto data_manager{ GetDataManagerPtr() };
    data_manager->SetDatabaseManager(m_options.database_path);
    try
    {
        data_manager->ProcessFile(m_options.model_file_path, m_model_key_tag);
        data_manager->ProcessFile(m_options.map_file_path, m_map_key_tag);
        if (m_options.is_simulation == true)
        {
            auto model_object{ data_manager->GetTypedDataObject<ModelObject>(m_model_key_tag) };
            UpdateModelObjectForSimulation(model_object.get());
        }
    }
    catch (const std::exception & e)
    {
        Logger::Log(LogLevel::Error,
            "PotentialAnalysisCommand::Execute() : " + std::string(e.what()));
        return false;
    }
    return true;
}

void PotentialAnalysisCommand::UpdateModelObjectForSimulation(ModelObject * model_object)
{
    Logger::Log(LogLevel::Debug, "PotentialAnalysisCommand::UpdateModelObjectForSimulation() called");
    if (model_object == nullptr) return;
    if (m_options.resolution_simulation == 0.0)
    {
        Logger::Log(LogLevel::Warning,
            "[Warning] The resolution of input simulated map hasn't been set.\n"
            "          Please give the corresponding resolution value for this map.\n"
            "          (-r, --sim-resolution)");
    }
    model_object->SetEmdID("Simulation");
    model_object->SetResolution(m_options.resolution_simulation);
    model_object->SetResolutionMethod("Blurring Width");
}

void PotentialAnalysisCommand::RunMapObjectPreprocessing(void)
{
    Logger::Log(LogLevel::Debug, "PotentialAnalysisCommand::RunMapObjectPreprocessing() called");
    ScopeTimer timer("PotentialAnalysisCommand::RunMapObjectPreprocessing");
    auto data_manager{ GetDataManagerPtr() };
    m_map_object = data_manager->GetTypedDataObject<MapObject>(m_map_key_tag);
    m_map_object->MapValueArrayNormalization();
}

void PotentialAnalysisCommand::RunModelObjectPreprocessing(void)
{
    Logger::Log(LogLevel::Debug, "PotentialAnalysisCommand::RunModelObjectPreprocessing() called");
    ScopeTimer timer("PotentialAnalysisCommand::RunModelObjectPreprocessing");
    auto data_manager{ GetDataManagerPtr() };
    m_model_object = data_manager->GetTypedDataObject<ModelObject>(m_model_key_tag);
    for (auto & atom : m_model_object->GetAtomList()) atom->SetSelectedFlag(true);
    for (auto & bond : m_model_object->GetBondList()) bond->SetSelectedFlag(true);
    m_model_object->FilterAtomFromSymmetry(m_options.is_asymmetry);
    m_model_object->FilterBondFromSymmetry(m_options.is_asymmetry);
    m_model_object->Update();
    for (auto & atom : m_model_object->GetSelectedAtomList())
    {
        auto local_potential_entry{ std::make_unique<LocalPotentialEntry>() };
        atom->AddLocalPotentialEntry(std::move(local_potential_entry));
    }
    for (auto & bond : m_model_object->GetSelectedBondList())
    {
        auto local_potential_entry{ std::make_unique<LocalPotentialEntry>() };
        bond->AddLocalPotentialEntry(std::move(local_potential_entry));
    }
    Logger::Log(LogLevel::Info,
        "Number of selected atom = " + std::to_string(m_model_object->GetNumberOfSelectedAtom()));
    Logger::Log(LogLevel::Info,
        "Number of selected bond = " + std::to_string(m_model_object->GetNumberOfSelectedBond()));
}

void PotentialAnalysisCommand::RunAtomMapValueSampling(void)
{
    Logger::Log(LogLevel::Debug, "PotentialAnalysisCommand::RunAtomMapValueSampling() called");
    ScopeTimer timer("PotentialAnalysisCommand::RunAtomMapValueSampling");
    if (m_map_object == nullptr) return;
    auto sampler{ std::make_unique<SphereSampler>() };
    sampler->SetSamplingSize(static_cast<unsigned int>(m_options.sampling_size));
    sampler->SetDistanceRangeMinimum(m_options.sampling_range_min);
    sampler->SetDistanceRangeMaximum(m_options.sampling_range_max);
    sampler->Print();
    
    const auto & atom_list{ m_model_object->GetSelectedAtomList() };
    auto atom_size{ atom_list.size() };
    size_t atom_count{ 0 };

#ifdef USE_OPENMP
    #pragma omp parallel num_threads(m_options.thread_size)
    {
        MapInterpolationVisitor interpolation_visitor{ sampler.get() };
        #pragma omp for
        for (size_t i = 0; i < atom_size; i++)
        {
            auto atom{ atom_list[i] };
            auto entry{ atom->GetLocalPotentialEntry() };
            interpolation_visitor.SetPosition(atom->GetPosition());
            m_map_object->Accept(&interpolation_visitor);
            entry->AddDistanceAndMapValueList(interpolation_visitor.MoveSamplingDataList());
            entry->AddBasisAndResponseEntryList(
                GausLinearTransformHelper::MapValueTransform(
                    entry->GetDistanceAndMapValueList(),
                    m_options.fit_range_min, m_options.fit_range_max)
            );
            if (m_options.use_training_alpha == false) entry->SetAlphaR(m_options.alpha_r);
            #pragma omp critical
            {
                atom_count++;
                Logger::ProgressPercent(atom_count, atom_size);
            }
        }
    }
#else
    MapInterpolationVisitor interpolation_visitor{ sampler.get() };
    for (size_t i = 0; i < atom_size; i++)
    {
        auto atom{ atom_list[i] };
        auto entry{ atom->GetLocalPotentialEntry() };
        interpolation_visitor.SetPosition(atom->GetPosition());
        m_map_object->Accept(&interpolation_visitor);
        entry->AddDistanceAndMapValueList(interpolation_visitor.GetSamplingDataList());
        entry->AddBasisAndResponseEntryList(
            GausLinearTransformHelper::MapValueTransform(
                entry->GetDistanceAndMapValueList(),
                m_options.fit_range_min, m_options.fit_range_max)
        );
        if (m_options.use_training_alpha == false) entry->SetAlphaR(m_options.alpha_r);
        atom_count++;
        Logger::ProgressPercent(atom_count, atom_size);
    }
#endif
}

void PotentialAnalysisCommand::RunBondMapValueSampling(void)
{
    Logger::Log(LogLevel::Debug, "PotentialAnalysisCommand::RunBondMapValueSampling() called");
    ScopeTimer timer("PotentialAnalysisCommand::RunBondMapValueSampling");
    if (m_map_object == nullptr) return;
    auto sampler{ std::make_unique<CylinderSampler>() };
    sampler->SetSamplingSize(static_cast<unsigned int>(m_options.sampling_size));
    sampler->SetDistanceRangeMinimum(m_options.sampling_range_min);
    sampler->SetDistanceRangeMaximum(m_options.sampling_range_max);
    sampler->SetHeight(m_options.sampling_height);
    sampler->Print();
    
    const auto & bond_list{ m_model_object->GetSelectedBondList() };
    auto bond_size{ bond_list.size() };
    size_t bond_count{ 0 };

#ifdef USE_OPENMP
    #pragma omp parallel num_threads(m_options.thread_size)
    {
        MapInterpolationVisitor interpolation_visitor{ sampler.get() };
        #pragma omp for
        for (size_t i = 0; i < bond_size; i++)
        {
            auto bond{ bond_list[i] };
            auto entry{ bond->GetLocalPotentialEntry() };
            auto bond_vector{ bond->GetBondVector() };
            auto bond_position{ bond->GetPosition() };
            auto adjusted_rate{ 0.0f };
            std::array<float, 3> adjusted_position{
                bond_position[0] + 0.5f * bond_vector[0] * adjusted_rate,
                bond_position[1] + 0.5f * bond_vector[1] * adjusted_rate,
                bond_position[2] + 0.5f * bond_vector[2] * adjusted_rate
            };
            interpolation_visitor.SetPosition(adjusted_position);
            interpolation_visitor.SetAxisVector(bond_vector);
            m_map_object->Accept(&interpolation_visitor);
            entry->AddDistanceAndMapValueList(interpolation_visitor.MoveSamplingDataList());
            entry->AddBasisAndResponseEntryList(
                GausLinearTransformHelper::MapValueTransform(
                    entry->GetDistanceAndMapValueList(),
                    m_options.fit_range_min, m_options.fit_range_max)
            );
            entry->SetAlphaR(m_options.alpha_r);
            #pragma omp critical
            {
                bond_count++;
                Logger::ProgressPercent(bond_count, bond_size);
            }
        }
    }
#else
    MapInterpolationVisitor interpolation_visitor{ sampler.get() };
    for (size_t i = 0; i < bond_size; i++)
    {
        auto bond{ bond_list[i] };
        auto entry{ bond->GetLocalPotentialEntry() };
        interpolation_visitor.SetPosition(bond->GetPosition());
        interpolation_visitor.SetAxisVector(bond->GetBondVector());
        m_map_object->Accept(&interpolation_visitor);
        entry->AddDistanceAndMapValueList(interpolation_visitor.GetSamplingDataList());
        entry->AddBasisAndResponseEntryList(
            GausLinearTransformHelper::MapValueTransform(
                entry->GetDistanceAndMapValueList(),
                m_options.fit_range_min, m_options.fit_range_max)
        );
        entry->SetAlphaR(m_options.alpha_r);
        bond_count++;
        Logger::ProgressPercent(bond_count, bond_size);
    }
#endif
}

void PotentialAnalysisCommand::RunAtomGroupClassification(void)
{
    Logger::Log(LogLevel::Debug, "RunAtomGroupClassification() called");
    ScopeTimer timer("RunAtomGroupClassification");
    if (m_map_object == nullptr) return;

    Logger::Log(LogLevel::Info, "Atom Classification Summary:");
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        auto group_potential_entry( std::make_unique<GroupPotentialEntry>() );
        for (auto atom : m_model_object->GetSelectedAtomList())
        {
            auto group_key{ AtomClassifier::GetGroupKeyInClass(atom, class_key) };
            group_potential_entry->AddAtomObjectPtr(group_key, atom);
            group_potential_entry->InsertGroupKey(group_key);
        }
        auto group_size{ group_potential_entry->GetGroupKeySet().size() };
        m_model_object->AddAtomGroupPotentialEntry(class_key, group_potential_entry);
        Logger::Log(LogLevel::Info,
            " - Class type: " + class_key + " include " + std::to_string(group_size) + " groups.");
    }
}

void PotentialAnalysisCommand::RunBondGroupClassification(void)
{
    Logger::Log(LogLevel::Debug, "RunBondGroupClassification() called");
    ScopeTimer timer("RunBondGroupClassification");
    if (m_map_object == nullptr) return;

    Logger::Log(LogLevel::Info, "Bond Classification Summary:");
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupBondClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupBondClassKey(i) };
        auto group_potential_entry( std::make_unique<GroupPotentialEntry>() );
        for (auto bond : m_model_object->GetSelectedBondList())
        {
            auto group_key{ BondClassifier::GetGroupKeyInClass(bond, class_key) };
            group_potential_entry->AddBondObjectPtr(group_key, bond);
            group_potential_entry->InsertGroupKey(group_key);
        }
        auto group_size{ group_potential_entry->GetGroupKeySet().size() };
        m_model_object->AddBondGroupPotentialEntry(class_key, group_potential_entry);
        Logger::Log(LogLevel::Info,
            " - Class type: " + class_key + " include " + std::to_string(group_size) + " groups.");
    }
}

void PotentialAnalysisCommand::RunAtomAlphaTraining(void)
{
    Logger::Log(LogLevel::Debug, "PotentialAnalysisCommand::RunAtomAlphaTraining() called");
    ScopeTimer timer("PotentialAnalysisCommand::RunAtomAlphaTraining");
    if (m_map_object == nullptr) return;

    const size_t subset_size_alpha_r{ 5 };
    std::vector<double> alpha_r_list{ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5 };
    auto ordered_alpha_r_list{ alpha_r_list };
    std::sort(ordered_alpha_r_list.begin(), ordered_alpha_r_list.end());
    
    // Alpha_R Training
    auto & atom_list{ m_model_object->GetSelectedAtomList() };
    std::vector<AtomObject *> selected_atom_list;
    selected_atom_list.reserve(atom_list.size());
    for (auto & atom : atom_list)
    {
        if (atom->IsMainChainAtom() == false) continue;
        selected_atom_list.emplace_back(atom);
    }
    selected_atom_list.shrink_to_fit();
    
    auto selected_atom_size{ selected_atom_list.size() };
    Logger::Log(LogLevel::Info,
        "Run Alpha_R Training with "+ std::to_string(selected_atom_size) +" atoms.");
    auto alpha_r{ TrainUniversalAlphaR(
        selected_atom_list, subset_size_alpha_r, ordered_alpha_r_list)
    };
    
    for (auto & atom : atom_list)
    {
        atom->GetLocalPotentialEntry()->SetAlphaR(alpha_r);
    }
    
    // Alpha_G Training
    RunLocalAtomFitting();
    const size_t subset_size_alpha_g{ 10 };
    std::vector<double> alpha_g_list{ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5 };
    auto ordered_alpha_g_list{ alpha_g_list };
    std::sort(ordered_alpha_g_list.begin(), ordered_alpha_g_list.end());

    std::vector<std::vector<AtomObject *>> atom_list_set;
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        auto group_potential_entry{ m_model_object->GetAtomGroupPotentialEntry(class_key) };
        const auto & group_key_set{ group_potential_entry->GetGroupKeySet() };
        for (auto group_key : group_key_set)
        {
            auto & group_atom_list{ group_potential_entry->GetAtomObjectPtrList(group_key) };
            if (group_atom_list.size() < 10) continue;
            atom_list_set.emplace_back(group_atom_list);
        }
    }

    auto selected_group_size{ atom_list_set.size() };
    Logger::Log(LogLevel::Info,
        "Run Alpha_G Training with "+ std::to_string(selected_group_size) +" groups.");
    auto alpha_g{ TrainUniversalAlphaG(
        atom_list_set, subset_size_alpha_g, ordered_alpha_g_list)
    };

    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        auto group_potential_entry{ m_model_object->GetAtomGroupPotentialEntry(class_key) };
        const auto & group_key_set{ group_potential_entry->GetGroupKeySet() };
        for (auto group_key : group_key_set)
        {
            group_potential_entry->AddAlphaG(group_key, alpha_g);
        }
    }
}

std::vector<double> PotentialAnalysisCommand::TrainAlphaR(
    const AtomObject * atom, const size_t group_size, const std::vector<double> & alpha_list)
{
    auto local_entry{ atom->GetLocalPotentialEntry() };
    const int basis_size{ 2 };
    std::vector<std::vector<Eigen::VectorXd>> grouped_data_list(group_size);
    auto total_entry_size{ static_cast<size_t>(local_entry->GetBasisAndResponseEntryListSize()) };
    size_t entries_in_atom_size{ total_entry_size / group_size + 1};
    for (size_t i = 0; i < group_size; i++)
    {
        grouped_data_list[i].reserve(entries_in_atom_size);
    }

    size_t count{ 0 };
    for (auto & data_entry : local_entry->GetBasisAndResponseEntryList())
    {
        auto group_index{ count % group_size };
        grouped_data_list[group_index].emplace_back(data_entry);
        count++;
    }

    std::vector<std::vector<Eigen::VectorXd>> data_test(group_size);
    std::vector<std::vector<Eigen::VectorXd>> data_training(group_size);
    for (size_t i = 0; i < group_size; i++)
    {
        data_test[i].reserve(entries_in_atom_size);
        data_training[i].reserve(total_entry_size - entries_in_atom_size);
        data_test[i].insert(
            data_test[i].end(), grouped_data_list[i].begin(), grouped_data_list[i].end());
        for (size_t j = 0; j < group_size; j++)
        {
            if (i == j) continue;
            data_training[i].insert(
                data_training[i].end(), grouped_data_list[j].begin(), grouped_data_list[j].end());
        }
    }
    
    auto alpha_size{ alpha_list.size() };
    std::vector<double> error_sum_list;
    error_sum_list.resize(alpha_size, 0.0);
    for (size_t p = 0; p < alpha_size; p++)
    {
        auto alpha{ alpha_list.at(static_cast<size_t>(p)) };
        auto beta_error_sum{ 0.0 };
        for (size_t i = 0; i < group_size; i++)
        {
            auto estimator{ std::make_unique<HRLModelHelper>(basis_size, 1) };
            estimator->SetQuietMode();
            estimator->SetThreadSize(1);
            Eigen::VectorXd beta_ols_test;
            Eigen::VectorXd beta_mdpde_test;
            double sigma_square_test;
            Eigen::DiagonalMatrix<double, Eigen::Dynamic> W_test;
            Eigen::DiagonalMatrix<double, Eigen::Dynamic> capital_sigma_test;
            estimator->RunBetaMDPDE(
                data_test.at(i), alpha, beta_ols_test, beta_mdpde_test,
                sigma_square_test, W_test, capital_sigma_test);

            Eigen::VectorXd beta_ols_training;
            Eigen::VectorXd beta_mdpde_training;
            double sigma_square_training;
            Eigen::DiagonalMatrix<double, Eigen::Dynamic> W_training;
            Eigen::DiagonalMatrix<double, Eigen::Dynamic> capital_sigma_training;
            estimator->RunBetaMDPDE(
                data_training.at(i), alpha, beta_ols_training,
                beta_mdpde_training, sigma_square_training, W_training, capital_sigma_training);

            beta_error_sum += (beta_mdpde_test - beta_mdpde_training).norm();
        }
        error_sum_list[p] = beta_error_sum;
    }

    return error_sum_list;
}

std::vector<double> PotentialAnalysisCommand::TrainAlphaG(
    const std::vector<AtomObject *> & atom_list,
    const size_t group_size, const std::vector<double> & alpha_list)
{
    if (atom_list.size() < 10)
    {
        std::vector<double> empty_result(alpha_list.size(), 0.0);
        return empty_result;
    }
    auto atom_in_half_size{ atom_list.size() / 2 };
    const int basis_size{ 2 };
    
    std::vector<std::vector<Eigen::VectorXd>> data_test(group_size);
    std::vector<std::vector<Eigen::VectorXd>> data_training(group_size);
    for (size_t i = 0; i < group_size; i++)
    {
        // Randomly pick the half of atoms into test set and training set for each group
        std::vector<AtomObject *> data1, data2;
        std::vector<AtomObject *> shuffled{ atom_list };
        std::shuffle(shuffled.begin(), shuffled.end(), std::mt19937{std::random_device{}()});
        auto diff{ static_cast<std::vector<AtomObject *>::difference_type>(atom_in_half_size) };
        data1.assign(shuffled.begin(), shuffled.begin() + diff);
        data2.assign(shuffled.begin() + diff, shuffled.end());
        data_test[i].reserve(data1.size());
        data_training[i].reserve(data2.size());
        for (auto atom : data1)
        {
            data_test[i].emplace_back(atom->GetLocalPotentialEntry()->GetBetaEstimateMDPDE());
        }
        for (auto atom : data2)
        {
            data_training[i].emplace_back(atom->GetLocalPotentialEntry()->GetBetaEstimateMDPDE());
        }
    }

    auto alpha_size{ alpha_list.size() };
    std::vector<double> error_sum_list;
    error_sum_list.resize(alpha_size, 0.0);
    for (size_t p = 0; p < alpha_size; p++)
    {
        auto alpha{ alpha_list.at(p) };
        auto mu_error_sum{ 0.0 };
        for (size_t i = 0; i < group_size; i++)
        {
            auto estimator{ std::make_unique<HRLModelHelper>(basis_size, 1) };
            estimator->SetQuietMode();
            estimator->SetThreadSize(1);
            Eigen::VectorXd mu_mdpde_test;
            Eigen::ArrayXd omega_array_test;
            double omega_sum_test;
            Eigen::MatrixXd capital_lambda_test;
            std::vector<Eigen::MatrixXd> member_capital_lambda_list_test;
            estimator->RunMuMDPDE(
                data_test.at(i), alpha, mu_mdpde_test,
                omega_array_test, omega_sum_test, capital_lambda_test,
                member_capital_lambda_list_test);

            Eigen::VectorXd mu_mdpde_training;
            Eigen::ArrayXd omega_array_training;
            double omega_sum_training;
            Eigen::MatrixXd capital_lambda_training;
            std::vector<Eigen::MatrixXd> member_capital_lambda_list_training;
            estimator->RunMuMDPDE(
                data_training.at(i), alpha, mu_mdpde_training,
                omega_array_training, omega_sum_training, capital_lambda_training,
                member_capital_lambda_list_training);

            mu_error_sum += (mu_mdpde_test - mu_mdpde_training).norm();
        }
        error_sum_list[p] = mu_error_sum;
    }

    return error_sum_list;
}

double PotentialAnalysisCommand::TrainUniversalAlphaR(
    const std::vector<AtomObject *> & atom_list,
    const size_t subset_size,
    const std::vector<double> & alpha_list)
{
    auto atom_size{ atom_list.size() };
    auto alpha_size{ static_cast<int>(alpha_list.size()) };

    std::atomic<size_t> atom_count{ 0 };
    Eigen::ArrayXd beta_error_sum_array{ Eigen::ArrayXd::Zero(alpha_size) };

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(m_options.thread_size)
#endif
    for (size_t i = 0; i < atom_size; i++)
    {
        auto error_list{ TrainAlphaR(atom_list[i], subset_size, alpha_list) };
        
#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            for (int p = 0; p < alpha_size; p++)
            {
                auto error{ error_list.at(static_cast<size_t>(p)) };
                beta_error_sum_array(p) += error;
            }
            atom_count++;
            Logger::ProgressPercent(atom_count, atom_size);
        }
    }

    int error_min_id;
    beta_error_sum_array.minCoeff(&error_min_id);
    auto alpha_r_error{ alpha_list.at(static_cast<size_t>(error_min_id)) };

    std::cout << "\nAlpha_R Training Results Summary:";
    std::cout << " - Minimum Beta Error Sum Alpha_R: " << alpha_r_error << "\n";

    return alpha_r_error;
}

double PotentialAnalysisCommand::TrainUniversalAlphaG(
    const std::vector<std::vector<AtomObject *>> & atom_list_set,
    const size_t subset_size,
    const std::vector<double> & alpha_list)
{
    auto alpha_size{ static_cast<int>(alpha_list.size()) };
    auto group_size{ atom_list_set.size() };
    std::atomic<size_t> group_count{ 0 };
    Eigen::ArrayXd mu_error_sum_array{ Eigen::ArrayXd::Zero(alpha_size) };

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(m_options.thread_size)
#endif
    for (size_t i = 0; i < group_size; i++)
    {
        auto & group_atom_list{ atom_list_set[i] };
        auto error_list{ TrainAlphaG(group_atom_list, subset_size, alpha_list) };
        
#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            for (int p = 0; p < alpha_size; p++)
            {
                auto error{ error_list.at(static_cast<size_t>(p)) };
                mu_error_sum_array(p) += error;
            }
            group_count++;
            Logger::ProgressPercent(group_count, group_size);
        }
    }

    int error_min_id;
    mu_error_sum_array.minCoeff(&error_min_id);
    auto alpha_g_error{ alpha_list.at(static_cast<size_t>(error_min_id)) };

    std::cout << "\nAlpha_G Training Results Summary:";
    std::cout << " - Minimum Mu Error Sum Alpha_G: " << alpha_g_error << "\n";

    return alpha_g_error;
}

void PotentialAnalysisCommand::RunLocalAtomFitting(void)
{
    Logger::Log(LogLevel::Debug, "PotentialAnalysisCommand::RunLocalAtomFitting() called");
    ScopeTimer timer("PotentialAnalysisCommand::RunLocalAtomFitting");
    if (m_model_object == nullptr) return;

    const int basis_size{ 2 };

    std::atomic<size_t> atom_count{ 0 };
    auto & selected_atom_list{ m_model_object->GetSelectedAtomList() };
    auto selected_atom_size{ selected_atom_list.size() };
    Logger::Log(LogLevel::Info,
        "Run Local atom fitting for "+ std::to_string(selected_atom_size) +" atoms.");
#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(m_options.thread_size)
#endif
    for (size_t i = 0; i < selected_atom_size; i++)
    {
        auto local_entry{ selected_atom_list[i]->GetLocalPotentialEntry() };
        auto & data_entry_list{ local_entry->GetBasisAndResponseEntryList() };
        auto estimator{ std::make_unique<HRLModelHelper>(basis_size, 1) };
        estimator->SetQuietMode();
        estimator->SetThreadSize(1);
        Eigen::VectorXd beta_ols;
        Eigen::VectorXd beta_mdpde;
        double sigma_square;
        Eigen::DiagonalMatrix<double, Eigen::Dynamic> W;
        Eigen::DiagonalMatrix<double, Eigen::Dynamic> capital_sigma;
        estimator->RunBetaMDPDE(
            data_entry_list, local_entry->GetAlphaR(), beta_ols, beta_mdpde,
            sigma_square, W, capital_sigma);

        local_entry->SetBetaEstimateOLS(beta_ols);
        local_entry->SetBetaEstimateMDPDE(beta_mdpde);
        local_entry->SetSigmaSquare(sigma_square);
        local_entry->SetDataWeight(W);
        local_entry->SetDataCovariance(capital_sigma);

        Eigen::VectorXd model_par_init{ Eigen::VectorXd::Zero(3) };
        model_par_init(0) = local_entry->GetMomentZeroEstimate();
        model_par_init(1) = local_entry->GetMomentTwoEstimate();
        auto gaus_ols{ GausLinearTransformHelper::BuildGaus3DModel(beta_ols, model_par_init) };
        auto gaus_mdpde{ GausLinearTransformHelper::BuildGaus3DModel(beta_mdpde, model_par_init) };
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

void PotentialAnalysisCommand::RunLocalBondFitting(void)
{
    Logger::Log(LogLevel::Debug, "PotentialAnalysisCommand::RunLocalBondFitting() called");
    ScopeTimer timer("PotentialAnalysisCommand::RunLocalBondFitting");
    if (m_model_object == nullptr) return;

    const int basis_size{ 2 };

    std::atomic<size_t> atom_count{ 0 };
    auto & selected_bond_list{ m_model_object->GetSelectedBondList() };
    auto selected_bond_size{ selected_bond_list.size() };
    Logger::Log(LogLevel::Info,
        "Run Local bond fitting for "+ std::to_string(selected_bond_size) +" bonds.");
#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(m_options.thread_size)
#endif
    for (size_t i = 0; i < selected_bond_size; i++)
    {
        auto local_entry{ selected_bond_list[i]->GetLocalPotentialEntry() };
        auto & data_entry_list{ local_entry->GetBasisAndResponseEntryList() };
        auto estimator{ std::make_unique<HRLModelHelper>(basis_size, 1) };
        estimator->SetQuietMode();
        estimator->SetThreadSize(1);
        Eigen::VectorXd beta_ols;
        Eigen::VectorXd beta_mdpde;
        double sigma_square;
        Eigen::DiagonalMatrix<double, Eigen::Dynamic> W;
        Eigen::DiagonalMatrix<double, Eigen::Dynamic> capital_sigma;
        estimator->RunBetaMDPDE(
            data_entry_list, local_entry->GetAlphaR(), beta_ols, beta_mdpde,
            sigma_square, W, capital_sigma);

        local_entry->SetBetaEstimateOLS(beta_ols);
        local_entry->SetBetaEstimateMDPDE(beta_mdpde);
        local_entry->SetSigmaSquare(sigma_square);
        local_entry->SetDataWeight(W);
        local_entry->SetDataCovariance(capital_sigma);

        Eigen::VectorXd model_par_init{ Eigen::VectorXd::Zero(3) };
        model_par_init(0) = local_entry->GetMomentZeroEstimate();
        model_par_init(1) = local_entry->GetMomentTwoEstimate();
        auto gaus_ols{ GausLinearTransformHelper::BuildGaus3DModel(beta_ols, model_par_init) };
        auto gaus_mdpde{ GausLinearTransformHelper::BuildGaus3DModel(beta_mdpde, model_par_init) };
        local_entry->AddGausEstimateOLS(gaus_ols(0), gaus_ols(1));
        local_entry->AddGausEstimateMDPDE(gaus_mdpde(0), gaus_mdpde(1));

        #ifdef USE_OPENMP
            #pragma omp critical
        #endif
        {
            atom_count++;
            Logger::ProgressPercent(atom_count, selected_bond_size);
        }
    }
}

void PotentialAnalysisCommand::RunAtomPotentialFitting(void)
{
    Logger::Log(LogLevel::Debug, "PotentialAnalysisCommand::RunAtomPotentialFitting() called");
    ScopeTimer timer("PotentialAnalysisCommand::RunAtomPotentialFitting");
    if (m_model_object == nullptr) return;
    const int basis_size{ 2 };
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        Logger::Log(LogLevel::Info, "Class type: " + class_key);

        // Group Atom Potential Fitting
        auto group_potential_entry{ m_model_object->GetAtomGroupPotentialEntry(class_key) };
        const auto & key_set{ group_potential_entry->GetGroupKeySet() };
        std::vector<GroupKey> group_keys(key_set.begin(), key_set.end());
        auto group_key_size{ group_keys.size() };
        std::atomic<size_t> key_count{ 0 };

#ifdef USE_OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(m_options.thread_size)
#endif
        for (size_t idx = 0; idx < group_key_size; idx++)
        {
            auto group_key{ group_keys[idx] };
            const auto & atom_list{ group_potential_entry->GetAtomObjectPtrList(group_key) };
            auto group_size{ atom_list.size() };
            std::vector<std::vector<Eigen::VectorXd>> data_entry_list;
            std::vector<Eigen::VectorXd> beta_mdpde_list;
            std::vector<double> sigma_square_list;
            std::vector<Eigen::DiagonalMatrix<double, Eigen::Dynamic>> data_weight_list;
            std::vector<Eigen::DiagonalMatrix<double, Eigen::Dynamic>> data_covariance_list;
            data_entry_list.reserve(group_size);
            beta_mdpde_list.reserve(group_size);
            sigma_square_list.reserve(group_size);
            data_weight_list.reserve(group_size);
            data_covariance_list.reserve(group_size);
            for (const auto & atom : atom_list)
            {
                auto entry{ atom->GetLocalPotentialEntry() };
                data_entry_list.emplace_back(entry->GetBasisAndResponseEntryList());
                beta_mdpde_list.emplace_back(entry->GetBetaEstimateMDPDE());
                sigma_square_list.emplace_back(entry->GetSigmaSquare());
                data_weight_list.emplace_back(entry->GetDataWeight());
                data_covariance_list.emplace_back(entry->GetDataCovariance());
            }
            auto model_estimator{ std::make_unique<HRLModelHelper>(basis_size, static_cast<int>(group_size)) };
            model_estimator->SetThreadSize(1);
            model_estimator->SetMemberDataEntriesList(data_entry_list);
            model_estimator->SetMemberBetaMDPDEList(
                beta_mdpde_list, sigma_square_list, data_weight_list, data_covariance_list
            );
            auto alpha_g{ (m_options.use_training_alpha) ?
                group_potential_entry->GetAlphaG(group_key) : m_options.alpha_g
            };
            model_estimator->RunGroupEstimation(alpha_g);

            auto gaus_group_mean{
                GausLinearTransformHelper::BuildGaus3DModel(model_estimator->GetMuVectorMean())
            };

            auto gaus_group_mdpde{
                GausLinearTransformHelper::BuildGaus3DModel(model_estimator->GetMuVectorMDPDE())
            };

            auto gaus_prior{
                GausLinearTransformHelper::BuildGaus3DModelWithVariance(
                    model_estimator->GetMuVectorPrior(), model_estimator->GetCapitalLambdaMatrix())
            };
            auto prior_estimate{ std::get<0>(gaus_prior) };
            auto prior_variance{ std::get<1>(gaus_prior) };

            auto count{ 0 };
            for (const auto & atom : atom_list)
            {
                auto atom_entry{ atom->GetLocalPotentialEntry() };
                const auto & beta_vector_posterior{ model_estimator->GetBetaPosterior(count) };
                const auto & sigma_matrix_posterior{ model_estimator->GetCapitalSigmaMatrixPosterior(count) };
                auto gaus_posterior{
                    GausLinearTransformHelper::BuildGaus3DModelWithVariance(
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
            model_estimator.reset();
            
#ifdef USE_OPENMP
            #pragma omp critical
#endif
            {
                group_potential_entry->AddGausEstimateMean(
                    group_key, gaus_group_mean(0), gaus_group_mean(1)
                );
                group_potential_entry->AddGausEstimateMDPDE(
                    group_key, gaus_group_mdpde(0), gaus_group_mdpde(1)
                );
                group_potential_entry->AddGausEstimatePrior(
                    group_key, prior_estimate(0), prior_estimate(1)
                );
                group_potential_entry->AddGausVariancePrior(
                    group_key, prior_variance(0), prior_variance(1)
                );
                group_potential_entry->AddAlphaG(group_key, alpha_g);
                key_count++;
                Logger::ProgressBar(key_count, group_key_size);
            }
        }
    }
}

void PotentialAnalysisCommand::RunBondPotentialFitting(void)
{
    Logger::Log(LogLevel::Debug, "PotentialAnalysisCommand::RunBondPotentialFitting() called");
    ScopeTimer timer("PotentialAnalysisCommand::RunBondPotentialFitting");
    if (m_model_object == nullptr) return;
    const int basis_size{ 2 };
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupBondClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupBondClassKey(i) };
        Logger::Log(LogLevel::Info, "Class type: " + class_key);

        // Group Bond Potential Fitting
        auto group_potential_entry{ m_model_object->GetBondGroupPotentialEntry(class_key) };
        const auto & key_set{ group_potential_entry->GetGroupKeySet() };
        std::vector<GroupKey> group_keys(key_set.begin(), key_set.end());
        auto group_key_size{ group_keys.size() };
        std::atomic<size_t> key_count{ 0 };

#ifdef USE_OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(m_options.thread_size)
#endif
        for (size_t idx = 0; idx < group_key_size; idx++)
        {
            auto group_key{ group_keys[idx] };
            const auto & bond_list{ group_potential_entry->GetBondObjectPtrList(group_key) };
            auto group_size{ bond_list.size() };
            std::vector<std::vector<Eigen::VectorXd>> data_entry_list;
            std::vector<Eigen::VectorXd> beta_mdpde_list;
            std::vector<double> sigma_square_list;
            std::vector<Eigen::DiagonalMatrix<double, Eigen::Dynamic>> data_weight_list;
            std::vector<Eigen::DiagonalMatrix<double, Eigen::Dynamic>> data_covariance_list;
            data_entry_list.reserve(group_size);
            beta_mdpde_list.reserve(group_size);
            sigma_square_list.reserve(group_size);
            data_weight_list.reserve(group_size);
            data_covariance_list.reserve(group_size);
            for (const auto & bond : bond_list)
            {
                auto entry{ bond->GetLocalPotentialEntry() };
                data_entry_list.emplace_back(entry->GetBasisAndResponseEntryList());
                beta_mdpde_list.emplace_back(entry->GetBetaEstimateMDPDE());
                sigma_square_list.emplace_back(entry->GetSigmaSquare());
                data_weight_list.emplace_back(entry->GetDataWeight());
                data_covariance_list.emplace_back(entry->GetDataCovariance());
            }
            auto model_estimator{ std::make_unique<HRLModelHelper>(basis_size, static_cast<int>(group_size)) };
            model_estimator->SetThreadSize(1);
            model_estimator->SetMemberDataEntriesList(data_entry_list);
            model_estimator->SetMemberBetaMDPDEList(
                beta_mdpde_list, sigma_square_list, data_weight_list, data_covariance_list
            );
            auto alpha_g{ m_options.alpha_g };
            model_estimator->RunGroupEstimation(alpha_g);

            auto gaus_group_mean{
                GausLinearTransformHelper::BuildGaus2DModel(model_estimator->GetMuVectorMean())
            };

            auto gaus_group_mdpde{
                GausLinearTransformHelper::BuildGaus2DModel(model_estimator->GetMuVectorMDPDE())
            };

            auto gaus_prior{
                GausLinearTransformHelper::BuildGaus2DModelWithVariance(
                    model_estimator->GetMuVectorPrior(), model_estimator->GetCapitalLambdaMatrix())
            };
            auto prior_estimate{ std::get<0>(gaus_prior) };
            auto prior_variance{ std::get<1>(gaus_prior) };

            auto count{ 0 };
            for (const auto & bond : bond_list)
            {
                auto bond_entry{ bond->GetLocalPotentialEntry() };
                const auto & beta_vector_posterior{ model_estimator->GetBetaPosterior(count) };
                const auto & sigma_matrix_posterior{ model_estimator->GetCapitalSigmaMatrixPosterior(count) };
                auto gaus_posterior{
                    GausLinearTransformHelper::BuildGaus2DModelWithVariance(
                        beta_vector_posterior, sigma_matrix_posterior)
                };
                auto posterior_estimate{ std::get<0>(gaus_posterior) };
                auto posterior_variance{ std::get<1>(gaus_posterior) };
                bond_entry->AddGausEstimatePosterior(class_key, posterior_estimate(0), posterior_estimate(1));
                bond_entry->AddGausVariancePosterior(class_key, posterior_variance(0), posterior_variance(1));
                bond_entry->AddOutlierTag(class_key, model_estimator->GetOutlierFlag(count));
                bond_entry->AddStatisticalDistance(class_key, model_estimator->GetStatisticalDistance(count));
                count++;
            }
            model_estimator.reset();
            
#ifdef USE_OPENMP
            #pragma omp critical
#endif
            {
                group_potential_entry->AddGausEstimateMean(
                    group_key, gaus_group_mean(0), gaus_group_mean(1)
                );
                group_potential_entry->AddGausEstimateMDPDE(
                    group_key, gaus_group_mdpde(0), gaus_group_mdpde(1)
                );
                group_potential_entry->AddGausEstimatePrior(
                    group_key, prior_estimate(0), prior_estimate(1)
                );
                group_potential_entry->AddGausVariancePrior(
                    group_key, prior_variance(0), prior_variance(1)
                );
                group_potential_entry->AddAlphaG(group_key, m_options.alpha_g);
                key_count++;
                Logger::ProgressBar(key_count, group_key_size);
            }
        }
    }
}

void PotentialAnalysisCommand::SaveDataObject(void)
{
    Logger::Log(LogLevel::Debug, "PotentialAnalysisCommand::SaveDataObject() called");
    ScopeTimer timer("PotentialAnalysisCommand::SaveDataObject");
    if (m_model_object == nullptr) return;

    auto data_manager{ GetDataManagerPtr() };
    data_manager->SaveDataObject(m_model_key_tag, m_options.saved_key_tag);

    for (auto atom : m_model_object->GetSelectedAtomList())
    {
        auto entry{ atom->GetLocalPotentialEntry() };
        if (entry != nullptr)
        {
            entry->ClearDistanceAndMapValueList();
        }
    }
}
