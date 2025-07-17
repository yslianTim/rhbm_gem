#include "PotentialAnalysisCommand.hpp"
#include "DataObjectManager.hpp"
#include "PotentialAnalysisVisitor.hpp"
#include "SphereSampler.hpp"
#include "ModelObject.hpp"
#include "Logger.hpp"

PotentialAnalysisCommand::PotentialAnalysisCommand(void) :
    m_sphere_sampler{ std::make_unique<SphereSampler>() }
{
    m_sphere_sampler->SetThreadSize(static_cast<unsigned int>(m_globals.thread_size));
    m_sphere_sampler->SetSamplingSize(static_cast<unsigned int>(m_options.sampling_size));
    m_sphere_sampler->SetDistanceRangeMinimum(m_options.sampling_range_min);
    m_sphere_sampler->SetDistanceRangeMaximum(m_options.sampling_range_max);
}

PotentialAnalysisCommand::~PotentialAnalysisCommand()
{
    
}

void PotentialAnalysisCommand::SetGlobalOptions(const GlobalOptions & globals)
{
    m_globals = globals;
    m_sphere_sampler->SetThreadSize(static_cast<unsigned int>(m_globals.thread_size));
}

void PotentialAnalysisCommand::RegisterCLIOptions(CLI::App * cmd)
{
    cmd->add_option("-a,--model", m_options.model_file_path,
        "Model file path")->required();
    cmd->add_option("-m,--map", m_options.map_file_path,
        "Map file path")->required();
    cmd->add_option("--simulation", m_options.is_simulation,
        "Simulation flag")->default_val(false);
    cmd->add_option("-r,--sim-resolution", m_options.resolution_simulation,
        "Set simulated map's resolution (blurring width)")->default_val(0.0);
    cmd->add_option("-k,--save-key", m_options.saved_key_tag,
        "New key tag for saving ModelObject results into database")->default_val("");
    cmd->add_option("--asymmetry", m_options.is_asymmetry,
        "Turn On/Off asymmetry flag")->default_val(false);
    cmd->add_option("-s,--sampling", m_options.sampling_size,
        "Number of sampling points per atom")->default_val(1500);
    cmd->add_option("--sampling-min", m_options.sampling_range_min,
        "Minimum sampling range")->default_val(0.0);
    cmd->add_option("--sampling-max", m_options.sampling_range_max,
        "Maximum sampling range")->default_val(1.5);
    cmd->add_option("--fit-min", m_options.fit_range_min,
        "Minimum fitting range")->default_val(0.0);
    cmd->add_option("--fit-max", m_options.fit_range_max,
        "Maximum fitting range")->default_val(1.0);
    cmd->add_option("--alpha-r", m_options.alpha_r,
        "Alpha value for R")->default_val(0.1);
    cmd->add_option("--alpha-g", m_options.alpha_g,
        "Alpha value for G")->default_val(0.2);
}

void PotentialAnalysisCommand::Execute(void)
{
    Logger::Log(LogLevel::Info, "PotentialAnalysisCommand::Execute() called.");

    auto data_manager{ std::make_unique<DataObjectManager>(m_globals.database_path) };
    try
    {
        data_manager->ProcessFile(m_options.model_file_path, "model");
        data_manager->ProcessFile(m_options.map_file_path, "map");
        if (m_options.is_simulation == true)
        {
            auto model_object{
                dynamic_cast<ModelObject *>(data_manager->GetDataObjectRef("model")) };
            if (model_object == nullptr)
            {
                throw std::runtime_error(
                    "PotentialAnalysisCommand::Execute(): invalid model object");
            }
            UpdateModelObjectForSimulation(model_object);
        }
    }
    catch (const std::exception & e)
    {
        Logger::Log(LogLevel::Error, e.what());
        return;
    }

    auto analyzer{ std::make_unique<PotentialAnalysisVisitor>(m_sphere_sampler.get()) };
    analyzer->SetAsymmetryFlag(m_options.is_asymmetry);
    analyzer->SetThreadSize(static_cast<unsigned int>(m_globals.thread_size));
    analyzer->SetModelObjectKeyTag("model");
    analyzer->SetMapObjectKeyTag("map");
    analyzer->SetFitRange(m_options.fit_range_min, m_options.fit_range_max);
    analyzer->SetAlphaR(m_options.alpha_r);
    analyzer->SetAlphaG(m_options.alpha_g);

    data_manager->Accept(analyzer.get());
    data_manager->SaveDataObject("model", m_options.saved_key_tag);
}

void PotentialAnalysisCommand::SetThreadSize(int value)
{
    m_globals.thread_size = value;
    m_sphere_sampler->SetThreadSize(static_cast<unsigned int>(value));
}

void PotentialAnalysisCommand::SetSamplingSize(int value)
{
    m_options.sampling_size = value;
    m_sphere_sampler->SetSamplingSize(static_cast<unsigned int>(value));
}

void PotentialAnalysisCommand::SetSamplingRangeMinimum(double value)
{
    m_options.sampling_range_min = value;
    m_sphere_sampler->SetDistanceRangeMinimum(value);
}

void PotentialAnalysisCommand::SetSamplingRangeMaximum(double value)
{
    m_options.sampling_range_max = value;
    m_sphere_sampler->SetDistanceRangeMaximum(value);
}

void PotentialAnalysisCommand::UpdateModelObjectForSimulation(ModelObject * model_object)
{
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
