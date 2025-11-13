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
    RunAtomAlphaTraining();
    RunAtomPotentialFitting();

    RunBondMapValueSampling();
    RunBondGroupClassification();
    RunBondPotentialFitting();
    SaveDataObject();
    return true;
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
            "Class type: " + class_key + " include " + std::to_string(group_size) + " groups.");
    }
}

void PotentialAnalysisCommand::RunBondGroupClassification(void)
{
    Logger::Log(LogLevel::Debug, "RunBondGroupClassification() called");
    ScopeTimer timer("RunBondGroupClassification");
    if (m_map_object == nullptr) return;

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
            "Class type: " + class_key + " include " + std::to_string(group_size) + " groups.");
    }
}

void PotentialAnalysisCommand::RunAtomAlphaTraining(void)
{
    Logger::Log(LogLevel::Debug, "PotentialAnalysisCommand::RunAtomAlphaTraining() called");
    ScopeTimer timer("PotentialAnalysisCommand::RunAtomAlphaTraining");
    if (m_map_object == nullptr) return;

    const size_t group_size{ 5 };
    std::vector<double> alpha_list{ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
    auto ordered_alpha_list{ alpha_list };
    std::sort(ordered_alpha_list.begin(), ordered_alpha_list.end());
    
    // Alpha_R Training
    auto & selected_atom_list{ m_model_object->GetSelectedAtomList() };
#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(m_options.thread_size)
#endif
    for (size_t i = 0; i < selected_atom_list.size(); i++)
    {
        auto atom{ selected_atom_list[i] };
        auto entry{ atom->GetLocalPotentialEntry() };
        auto alpha_r{ TrainAlphaR(atom, group_size, ordered_alpha_list) };
        entry->SetAlphaR(alpha_r);
        #ifdef USE_OPENMP
            #pragma omp critical
        #endif
        {
            Logger::ProgressPercent(i+1, m_model_object->GetNumberOfSelectedAtom());
        }
    }

    // Alpha_G Training
    const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(0) };
    auto group_potential_entry{ m_model_object->GetAtomGroupPotentialEntry(class_key) };
    auto classifier{ std::make_unique<AtomClassifier>() };
    auto group_key{ classifier->GetMainChainSimpleAtomClassGroupKey(0) };
    Logger::Log(LogLevel::Info,
        "Using group_key: " + std::to_string(group_key) + " for alpha_g training...");

    const auto & atom_list{ group_potential_entry->GetAtomObjectPtrList(group_key) };
    auto total_atom_size{ atom_list.size() };
    Logger::Log(LogLevel::Info, " - Include " + std::to_string(total_atom_size) + " atoms.");
    auto alpha_g{ TrainAlphaG(atom_list, 5, ordered_alpha_list) };
    Logger::Log(LogLevel::Info, "Alpha G = " + std::to_string(alpha_g));
}

double PotentialAnalysisCommand::TrainAlphaR(
    const AtomObject * atom, const size_t group_size, const std::vector<double> & alpha_list)
{
    auto local_entry{ atom->GetLocalPotentialEntry() };
    const int basis_size{ 2 };
    Eigen::VectorXd model_par_init{ Eigen::VectorXd::Zero(3) };
    std::map<size_t, std::vector<Eigen::VectorXd>> sampling_entry_list_map;
    auto total_entry_size{ static_cast<size_t>(local_entry->GetDistanceAndMapValueListSize()) };
    size_t entries_in_atom_size{ total_entry_size / group_size + 1};
    for (size_t i = 0; i < group_size; i++)
    {
        sampling_entry_list_map[i].reserve(entries_in_atom_size);
    }

    size_t count{ 0 };
    for (auto & data_entry : local_entry->GetDistanceAndMapValueList())
    {
        auto gaus_x{ static_cast<double>(std::get<0>(data_entry)) };
        auto gaus_y{ static_cast<double>(std::get<1>(data_entry)) };
        if (gaus_x < m_options.fit_range_min || gaus_x > m_options.fit_range_max) continue;
        if (gaus_y <= 0.0) continue;
        auto group_index{ count % group_size };
        sampling_entry_list_map.at(group_index).emplace_back(
            GausLinearTransformHelper::BuildLinearModelDataVector(
                gaus_x, gaus_y, model_par_init, basis_size)
        );
        count++;
    }
    size_t total_selected_entry_size{ count };

    std::map<size_t, std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>>> data_array_test_map;
    std::map<size_t, std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>>> data_array_training_map;
    for (size_t i = 0; i < group_size; i++)
    {
        auto test_set_entry_size{ sampling_entry_list_map.at(i).size() };
        data_array_test_map[i].emplace_back(
            sampling_entry_list_map.at(i), "test set - group:" + std::to_string(i));
        std::vector<Eigen::VectorXd> training_set;
        training_set.reserve(total_selected_entry_size - test_set_entry_size);
        for (size_t j = 0; j < group_size; j++)
        {
            if (j == i) continue;
            training_set.insert(
                training_set.end(),
                sampling_entry_list_map.at(j).begin(), sampling_entry_list_map.at(j).end());
        }
        data_array_training_map[i].emplace_back(
            std::move(training_set), "training set - group:" + std::to_string(i));
    }

    std::map<size_t, std::unique_ptr<HRLModelHelper>> estimator_test_map;
    std::map<size_t, std::unique_ptr<HRLModelHelper>> estimator_training_map;
    for (size_t i = 0; i < group_size; i++)
    {
        estimator_test_map[i] = std::make_unique<HRLModelHelper>(
            basis_size, static_cast<int>(data_array_test_map.at(i).size()));
        estimator_training_map[i] = std::make_unique<HRLModelHelper>(
            basis_size, static_cast<int>(data_array_training_map.at(i).size()));
        estimator_test_map.at(i)->SetQuietMode();
        estimator_test_map.at(i)->SetThreadSize(1);
        estimator_test_map.at(i)->SetDataArray(std::move(data_array_test_map.at(i)));
        estimator_training_map.at(i)->SetQuietMode();
        estimator_training_map.at(i)->SetThreadSize(1);
        estimator_training_map.at(i)->SetDataArray(std::move(data_array_training_map.at(i)));
    }
    
    auto alpha_size{ static_cast<int>(alpha_list.size()) };
    Eigen::ArrayXd beta_error_sum_array{ Eigen::ArrayXd::Zero(alpha_size) };
    Eigen::ArrayXd beta_variance_trace_array{ Eigen::ArrayXd::Zero(alpha_size) };
    for (int p = 0; p < alpha_size; p++)
    {
        auto alpha{ alpha_list.at(static_cast<size_t>(p)) };
        Eigen::VectorXd beta_mean_training{ Eigen::VectorXd::Zero(basis_size) };
        std::map<size_t, Eigen::VectorXd> beta_mdpde_training_map;
        auto beta_error_sum{ 0.0 };
        for (size_t i = 0; i < group_size; i++)
        {
            auto estimator_test{ estimator_test_map.at(i).get() };
            auto estimator_training{ estimator_training_map.at(i).get() };
            estimator_test->SetUniversalAlphaR(alpha);
            estimator_training->SetUniversalAlphaR(alpha);
            estimator_test->RunEstimation(0.0);
            estimator_training->RunEstimation(0.0);
            auto beta_mdpde_test{ estimator_test->GetBetaMatrixMDPDE(0) };
            auto beta_mdpde_training{ estimator_training->GetBetaMatrixMDPDE(0) };
            beta_mdpde_training_map.emplace(i, beta_mdpde_training);
            beta_error_sum += (beta_mdpde_test - beta_mdpde_training).norm();
            beta_mean_training += beta_mdpde_training;
        }
        beta_mean_training /= static_cast<double>(group_size);
        beta_error_sum_array(p) = beta_error_sum;
        
        // Calculate variance of beta
        Eigen::MatrixXd beta_variance_training{ Eigen::MatrixXd::Zero(basis_size, basis_size) };
        for (size_t i = 0; i < group_size; i++)
        {
            auto beta_mdpde_training{ beta_mdpde_training_map.at(i) };
            auto beta_deviation{ beta_mdpde_training - beta_mean_training };
            beta_variance_training += beta_deviation * beta_deviation.transpose();
        }
        beta_variance_training /= static_cast<double>(group_size - 1);
        auto trace{ beta_variance_training.trace() };
        beta_variance_trace_array(p) = trace;
        //std::cout << "Training : Variance trace = " << trace << std::endl;
        //std::cout << "Test : Error sum = " << beta_error_sum << std::endl;
    }

    int trace_min_id, error_min_id;
    beta_variance_trace_array.minCoeff(&trace_min_id);
    beta_error_sum_array.minCoeff(&error_min_id);
    //auto trace_min{ beta_variance_trace_array.minCoeff(&trace_min_id) };
    //auto error_min{ beta_error_sum_array.minCoeff(&error_min_id) };
    //std::cout << "Min trace = " << trace_min << " @ Id: "<< trace_min_id << std::endl;
    //std::cout << "Min error = " << error_min << " @ Id: " << error_min_id << std::endl;
    auto result_alpha_id{ std::max(trace_min_id, error_min_id) };

    return alpha_list.at(static_cast<size_t>(result_alpha_id));
}

double PotentialAnalysisCommand::TrainAlphaG(
    const std::vector<AtomObject *> & atom_list,
    const size_t group_size, const std::vector<double> & alpha_list)
{
    size_t atom_in_group_size{ atom_list.size() / group_size + 1};
    const int basis_size{ 2 };
    std::map<size_t, std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>>> total_data_array_map;
    std::map<size_t, std::vector<double>> total_alpha_r_map;
    for (size_t i = 0; i < group_size; i++)
    {
        total_data_array_map[i].reserve(atom_in_group_size);
        total_alpha_r_map[i].reserve(atom_in_group_size);
    }

    size_t count{ 0 };
    for (const auto & atom : atom_list)
    {
        auto entry{ atom->GetLocalPotentialEntry() };
        Eigen::VectorXd model_par_init{ Eigen::VectorXd::Zero(3) };
        std::vector<Eigen::VectorXd> sampling_entry_list;
        sampling_entry_list.reserve(static_cast<size_t>(entry->GetDistanceAndMapValueListSize()));
        for (auto & data_entry : entry->GetDistanceAndMapValueList())
        {
            auto gaus_x{ static_cast<double>(std::get<0>(data_entry)) };
            auto gaus_y{ static_cast<double>(std::get<1>(data_entry)) };
            if (gaus_x < m_options.fit_range_min || gaus_x > m_options.fit_range_max) continue;
            if (gaus_y <= 0.0) continue;
            sampling_entry_list.emplace_back(
                GausLinearTransformHelper::BuildLinearModelDataVector(
                    gaus_x, gaus_y, model_par_init, basis_size)
            );
        }
        auto group_index{ count % group_size };
        total_data_array_map.at(group_index).emplace_back(std::move(sampling_entry_list), atom->GetInfo());
        total_alpha_r_map.at(group_index).emplace_back(entry->GetAlphaR());
        count++;
    }
    
    std::map<size_t, std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>>> data_array_test_map;
    std::map<size_t, std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>>> data_array_training_map;
    std::map<size_t, std::vector<double>> alpha_r_list_test_map;
    std::map<size_t, std::vector<double>> alpha_r_list_training_map;
    for (size_t i = 0; i < group_size; i++)
    {
        auto test_set_atom_size{ total_data_array_map.at(i).size() };
        data_array_test_map.emplace(i, total_data_array_map.at(i));
        alpha_r_list_test_map.emplace(i, total_alpha_r_map.at(i));
        std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> data_array_training_set;
        std::vector<double> alpha_r_list_training_set;
        data_array_training_set.reserve(atom_list.size() - test_set_atom_size);
        alpha_r_list_training_set.reserve(atom_list.size() - test_set_atom_size);
        for (size_t j = 0; j < group_size; j++)
        {
            if (j == i) continue;
            data_array_training_set.insert(
                data_array_training_set.end(),
                total_data_array_map.at(j).begin(), total_data_array_map.at(j).end());
            alpha_r_list_training_set.insert(
                alpha_r_list_training_set.end(),
                total_alpha_r_map.at(j).begin(), total_alpha_r_map.at(j).end());
        }
        data_array_training_map.emplace(i, std::move(data_array_training_set));
        alpha_r_list_training_map.emplace(i, std::move(alpha_r_list_training_set));
    }

    std::map<size_t, std::unique_ptr<HRLModelHelper>> estimator_test_map;
    std::map<size_t, std::unique_ptr<HRLModelHelper>> estimator_training_map;
    for (size_t i = 0; i < group_size; i++)
    {
        estimator_test_map[i] = std::make_unique<HRLModelHelper>(
            basis_size, static_cast<int>(data_array_test_map.at(i).size()));
        estimator_training_map[i] = std::make_unique<HRLModelHelper>(
            basis_size, static_cast<int>(data_array_training_map.at(i).size()));
        estimator_test_map.at(i)->SetThreadSize(1);
        estimator_test_map.at(i)->SetDataArray(std::move(data_array_test_map.at(i)));
        estimator_test_map.at(i)->SetDedicateAlphaRList(std::move(alpha_r_list_test_map.at(i)));
        estimator_training_map.at(i)->SetThreadSize(1);
        estimator_training_map.at(i)->SetDataArray(std::move(data_array_training_map.at(i)));
        estimator_training_map.at(i)->SetDedicateAlphaRList(std::move(alpha_r_list_training_map.at(i)));
    }

    auto alpha_size{ static_cast<int>(alpha_list.size()) };
    Eigen::ArrayXd mu_error_sum_array{ Eigen::ArrayXd::Zero(alpha_size) };
    Eigen::ArrayXd mu_variance_trace_array{ Eigen::ArrayXd::Zero(alpha_size) };
    for (int p = 0; p < alpha_size; p++)
    {
        auto alpha{ alpha_list.at(static_cast<size_t>(p)) };
        Eigen::VectorXd mu_mean_training{ Eigen::VectorXd::Zero(basis_size) };
        std::map<size_t, Eigen::VectorXd> mu_mdpde_training_map;
        auto mu_error_sum{ 0.0 };
        for (size_t i = 0; i < group_size; i++)
        {
            auto estimator_test{ estimator_test_map.at(i).get() };
            auto estimator_training{ estimator_training_map.at(i).get() };
            estimator_test->RunEstimation(alpha);
            estimator_training->RunEstimation(alpha);
            auto mu_mdpde_test{ estimator_test->GetMuVectorMDPDE() };
            auto mu_mdpde_training{ estimator_training->GetMuVectorMDPDE() };
            mu_mdpde_training_map.emplace(i, mu_mdpde_training);
            mu_error_sum += (mu_mdpde_test - mu_mdpde_training).norm();
            mu_mean_training += mu_mdpde_training;
        }
        mu_mean_training /= static_cast<double>(group_size);
        mu_error_sum_array(p) = mu_error_sum;
        
        // Calculate variance of mu
        Eigen::MatrixXd mu_variance_training{ Eigen::MatrixXd::Zero(basis_size, basis_size) };
        for (size_t i = 0; i < group_size; i++)
        {
            auto mu_mdpde_training{ mu_mdpde_training_map.at(i) };
            auto mu_deviation{ mu_mdpde_training - mu_mean_training };
            mu_variance_training += mu_deviation * mu_deviation.transpose();
        }
        mu_variance_training /= static_cast<double>(group_size - 1);
        auto trace{ mu_variance_training.trace() };
        mu_variance_trace_array(p) = trace;
    }

    int trace_min_id, error_min_id;
    auto trace_min{ mu_variance_trace_array.minCoeff(&trace_min_id) };
    auto error_min{ mu_error_sum_array.minCoeff(&error_min_id) };
    std::cout << "Min trace = " << trace_min << " @ Id: "<< trace_min_id << std::endl;
    std::cout << "Min error = " << error_min << " @ Id: " << error_min_id << std::endl;
    auto result_alpha_id{ std::max(trace_min_id, error_min_id) };

    return alpha_list.at(static_cast<size_t>(result_alpha_id));
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
            std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> data_array;
            std::vector<double> alpha_r_list;
            data_array.reserve(group_size);
            alpha_r_list.reserve(group_size);
            for (const auto & atom : atom_list)
            {
                auto entry{ atom->GetLocalPotentialEntry() };
                Eigen::VectorXd model_par_init{ Eigen::VectorXd::Zero(3) };
                model_par_init(0) = entry->GetMomentZeroEstimate();
                model_par_init(1) = entry->GetMomentTwoEstimate();
                std::vector<Eigen::VectorXd> sampling_entry_list;
                sampling_entry_list.reserve(static_cast<size_t>(entry->GetDistanceAndMapValueListSize()));
                for (auto & data_entry : entry->GetDistanceAndMapValueList())
                {
                    auto gaus_x{ static_cast<double>(std::get<0>(data_entry)) };
                    auto gaus_y{ static_cast<double>(std::get<1>(data_entry)) };
                    if (gaus_x < m_options.fit_range_min || gaus_x > m_options.fit_range_max) continue;
                    if (gaus_y <= 0.0) continue;
                    sampling_entry_list.emplace_back(
                        GausLinearTransformHelper::BuildLinearModelDataVector(
                            gaus_x, gaus_y, model_par_init, basis_size)
                    );
                }
                data_array.emplace_back(std::move(sampling_entry_list), atom->GetInfo());
                alpha_r_list.emplace_back(entry->GetAlphaR());
            }
            auto model_estimator{ std::make_unique<HRLModelHelper>(basis_size, static_cast<int>(group_size)) };
            model_estimator->SetThreadSize(1);
            model_estimator->SetDataArray(std::move(data_array));
            //model_estimator->SetUniversalAlphaR(m_options.alpha_r);
            model_estimator->SetDedicateAlphaRList(alpha_r_list);
            model_estimator->RunEstimation(m_options.alpha_g);

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
                Eigen::VectorXd model_par_init{ Eigen::VectorXd::Zero(3) };
                model_par_init(0) = atom_entry->GetMomentZeroEstimate();
                model_par_init(1) = atom_entry->GetMomentTwoEstimate();
                const auto & beta_vector_ols{ model_estimator->GetBetaMatrixOLS(count) };
                auto gaus_ols{ GausLinearTransformHelper::BuildGaus3DModel(beta_vector_ols, model_par_init) };
                atom_entry->AddGausEstimateOLS(gaus_ols(0), gaus_ols(1));

                const auto & beta_vector_mdpde{ model_estimator->GetBetaMatrixMDPDE(count) };
                auto gaus_mdpde{ GausLinearTransformHelper::BuildGaus3DModel(beta_vector_mdpde, model_par_init) };
                atom_entry->AddGausEstimateMDPDE(gaus_mdpde(0), gaus_mdpde(1));

                const auto & beta_vector_posterior{ model_estimator->GetBetaMatrixPosterior(count) };
                auto sigma_matrix_posterior{ model_estimator->GetCapitalSigmaMatrixPosterior(count) };
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
            std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> data_array;
            data_array.reserve(group_size);
            for (const auto & bond : bond_list)
            {
                auto entry{ bond->GetLocalPotentialEntry() };
                Eigen::VectorXd model_par_init{ Eigen::VectorXd::Zero(3) };
                std::vector<Eigen::VectorXd> sampling_entry_list;
                sampling_entry_list.reserve(static_cast<size_t>(entry->GetDistanceAndMapValueListSize()));
                for (auto & data_entry : entry->GetDistanceAndMapValueList())
                {
                    auto gaus_x{ static_cast<double>(std::get<0>(data_entry)) };
                    auto gaus_y{ static_cast<double>(std::get<1>(data_entry)) };
                    if (gaus_x < m_options.fit_range_min || gaus_x > m_options.fit_range_max) continue;
                    if (gaus_y <= 0.0) continue;
                    sampling_entry_list.emplace_back(
                        GausLinearTransformHelper::BuildLinearModelDataVector(gaus_x, gaus_y, model_par_init)
                    );
                }
                data_array.emplace_back(std::move(sampling_entry_list), bond->GetInfo());
            }
            auto model_estimator{ std::make_unique<HRLModelHelper>(2, static_cast<int>(group_size)) };
            model_estimator->SetThreadSize(1);
            model_estimator->SetDataArray(std::move(data_array));
            model_estimator->SetUniversalAlphaR(m_options.alpha_r);
            model_estimator->RunEstimation(m_options.alpha_g);

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
                const auto & beta_vector_ols{ model_estimator->GetBetaMatrixOLS(count) };
                auto gaus_ols{ GausLinearTransformHelper::BuildGaus2DModel(beta_vector_ols) };
                bond_entry->AddGausEstimateOLS(gaus_ols(0), gaus_ols(1));

                const auto & beta_vector_mdpde{ model_estimator->GetBetaMatrixMDPDE(count) };
                auto gaus_mdpde{ GausLinearTransformHelper::BuildGaus2DModel(beta_vector_mdpde) };
                bond_entry->AddGausEstimateMDPDE(gaus_mdpde(0), gaus_mdpde(1));

                const auto & beta_vector_posterior{ model_estimator->GetBetaMatrixPosterior(count) };
                auto sigma_matrix_posterior{ model_estimator->GetCapitalSigmaMatrixPosterior(count) };
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
