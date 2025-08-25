#include "PotentialAnalysisCommand.hpp"
#include "DataObjectManager.hpp"
#include "AtomObject.hpp"
#include "MapObject.hpp"
#include "ModelObject.hpp"
#include "MapInterpolationVisitor.hpp"
#include "HRLModelHelper.hpp"
#include "ScopeTimer.hpp"
#include "FilePathHelper.hpp"
#include "AtomicPotentialEntry.hpp"
#include "GroupPotentialEntry.hpp"
#include "AtomicInfoHelper.hpp"
#include "AtomClassifier.hpp"
#include "GausLinearTransformHelper.hpp"
#include "SphereSampler.hpp"
#include "Logger.hpp"
#include "CommandRegistry.hpp"

#include <unordered_set>
#include <tuple>

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
    m_sphere_sampler{ std::make_unique<SphereSampler>() },
    m_map_object{ nullptr }, m_model_object{ nullptr }
{
    Logger::Log(LogLevel::Debug, "PotentialAnalysisCommand::PotentialAnalysisCommand() called.");
}

PotentialAnalysisCommand::~PotentialAnalysisCommand()
{
    Logger::Log(LogLevel::Debug, "PotentialAnalysisCommand::~PotentialAnalysisCommand() called.");
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
    RunMapValueSampling();
    RunPotentialFitting();
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
    for (auto & atom : m_model_object->GetComponentsList())
    {
        atom->SetSelectedFlag(true);
    }
    m_model_object->FilterAtomFromSymmetry(m_options.is_asymmetry);
    m_model_object->Update();
    for (auto & atom : m_model_object->GetSelectedAtomList())
    {
        auto atom_potential_entry{ std::make_unique<AtomicPotentialEntry>() };
        atom->AddAtomicPotentialEntry(std::move(atom_potential_entry));
    }
    Logger::Log(LogLevel::Info,
        "Number of selected atom = " + std::to_string(m_model_object->GetNumberOfSelectedAtom()));
}

void PotentialAnalysisCommand::RunMapValueSampling(void)
{
    Logger::Log(LogLevel::Debug, "PotentialAnalysisCommand::RunMapValueSampling() called");
    ScopeTimer timer("PotentialAnalysisCommand::RunMapValueSampling");
    if (m_map_object == nullptr || m_sphere_sampler == nullptr) return;
    m_sphere_sampler->SetSamplingSize(static_cast<unsigned int>(m_options.sampling_size));
    m_sphere_sampler->SetDistanceRangeMinimum(m_options.sampling_range_min);
    m_sphere_sampler->SetDistanceRangeMaximum(m_options.sampling_range_max);
    m_sphere_sampler->Print();
    
    const auto & atom_list{ m_model_object->GetSelectedAtomList() };
    auto atom_size{ atom_list.size() };
    size_t atom_count{ 0 };

#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(m_options.thread_size)
#endif
    for (size_t i = 0; i < atom_size; i++)
    {
        MapInterpolationVisitor interpolation_visitor{ m_sphere_sampler.get() };
        auto atom{ atom_list[i] };
        auto entry{ atom->GetAtomicPotentialEntry() };
        interpolation_visitor.SetPosition(atom->GetPosition());
        m_map_object->Accept(&interpolation_visitor);
        entry->AddDistanceAndMapValueList(interpolation_visitor.GetSamplingDataList());
#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            atom_count++;
            Logger::ProgressBar(atom_count, atom_size);
        }
    }
}

void PotentialAnalysisCommand::RunPotentialFitting(void)
{
    Logger::Log(LogLevel::Debug, "PotentialAnalysisCommand::RunPotentialFitting() called");
    ScopeTimer timer("PotentialAnalysisCommand::RunPotentialFitting");
    if (m_model_object == nullptr) return;
    for (size_t i = 0; i < AtomicInfoHelper::GetGroupClassCount(); i++)
    {
        const auto & class_key{ AtomicInfoHelper::GetGroupClassKey(i) };
        Logger::Log(LogLevel::Info, "Class type: " + class_key);

        // Atom Classification
        std::unordered_set<uint64_t> group_key_set;
        auto group_potential_entry( std::make_unique<GroupPotentialEntry>() );
        for (auto atom : m_model_object->GetSelectedAtomList())
        {
            auto group_key{ AtomClassifier::GetGroupKeyInClass(atom, class_key) };
            group_potential_entry->AddAtomObjectPtr(group_key, atom);
            group_potential_entry->InsertGroupKey(group_key);
            group_key_set.insert(group_key);
        }

        // Group Potential Fitting
        auto group_key_size{ group_key_set.size() };
        size_t key_count{ 1 };
        for (const auto & group_key : group_key_set)
        {
            Logger::ProgressBar(key_count, group_key_size);
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
                    if (gaus_x < m_options.fit_range_min || gaus_x > m_options.fit_range_max) continue;
                    if (gaus_y <= 0.0) continue;
                    sampling_entry_list.emplace_back(
                        GausLinearTransformHelper::BuildLinearModelDataVector(gaus_x, gaus_y)
                    );
                }
                data_array.emplace_back(std::make_tuple(sampling_entry_list, atom->GetInfo()));
            }
            auto model_estimator{ std::make_unique<HRLModelHelper>(2, static_cast<int>(group_size)) };
            model_estimator->SetDataArray(data_array);
            model_estimator->RunEstimation(m_options.alpha_r, m_options.alpha_g);

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
                const auto & beta_vector_ols{ model_estimator->GetBetaMatrixOLS(count) };
                auto gaus_ols{ GausLinearTransformHelper::BuildGausModel(beta_vector_ols) };
                atom_entry->AddGausEstimateOLS(gaus_ols(0), gaus_ols(1));

                const auto & beta_vector_mdpde{ model_estimator->GetBetaMatrixMDPDE(count) };
                auto gaus_mdpde{ GausLinearTransformHelper::BuildGausModel(beta_vector_mdpde) };
                atom_entry->AddGausEstimateMDPDE(gaus_mdpde(0), gaus_mdpde(1));

                const auto & beta_vector_posterior{ model_estimator->GetBetaMatrixPosterior(count) };
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
            key_count++;
        }
        m_model_object->AddGroupPotentialEntry(class_key, group_potential_entry);
    }

    auto data_manager{ GetDataManagerPtr() };
    data_manager->SaveDataObject(m_model_key_tag, m_options.saved_key_tag);
}
