#include "CommandBase.hpp"
#include "FilePathHelper.hpp"
#include "Logger.hpp"

#include <CLI/CLI.hpp>
#include <system_error>

namespace rhbm_gem {

void CommandBase::RegisterCLIOptions(CLI::App * command)
{
    RegisterCLIOptionsBasic(command);
    RegisterCLIOptionsExtend(command);
}

void CommandBase::RegisterCLIOptionsBasic(CLI::App * command)
{
    auto & options{ GetOptions() };
    command->add_option_function<int>("-j,--jobs",
        [&](int value) { SetThreadSize(value); },
        "Number of threads")->default_val(options.thread_size);
    command->add_option_function<int>("-v,--verbose",
        [&](int value) { SetVerboseLevel(value); },
        "Verbose level")->default_val(options.verbose_level);
    command->add_option_function<std::string>("-d,--database",
        [&](const std::string & value) { SetDatabasePath(value); },
        "Database file path")->default_val(options.database_path.string());
    command->add_option_function<std::string>("-o,--folder",
        [&](const std::string & value) { SetFolderPath(value); },
        "folder path for output files")->default_val(options.folder_path.string());
}

void CommandBase::SetThreadSize(int value)
{
    auto & options{ GetOptions() };
    if (value < 1)
    {
        Logger::Log(LogLevel::Warning,
            "Thread size must be positive. Using 1 instead.");
        options.thread_size = 1;
    }
    else
    {
        options.thread_size = value;
    }
}

void CommandBase::SetVerboseLevel(int value)
{
    auto & options{ GetOptions() };
    if (value < static_cast<int>(LogLevel::Error) || value > static_cast<int>(LogLevel::Debug))
    {
        Logger::Log(LogLevel::Warning,
            "Invalid verbose level: " + std::to_string(value) +
            ", using default level 3 [Info]");
        options.verbose_level = static_cast<int>(LogLevel::Info);
    }
    else
    {
        options.verbose_level = value;
    }
}

void CommandBase::SetDatabasePath(const std::filesystem::path & path)
{
    auto & options{ GetOptions() };
    options.database_path = path;

    auto parent_path{ options.database_path.parent_path() };
    if (!parent_path.empty() && !std::filesystem::exists(parent_path))
    {
        std::error_code error_code;
        std::filesystem::create_directories(parent_path, error_code);
        if (error_code)
        {
            Logger::Log(LogLevel::Error,
                "Failed to create directory: " + parent_path.string() +
                " error code = (" + error_code.message() + ")");
            m_valiate_options = false;
            return;
        }
    }
}

void CommandBase::SetFolderPath(const std::filesystem::path & path)
{
    auto & options{ GetOptions() };
    options.folder_path = std::filesystem::path(FilePathHelper::EnsureTrailingSlash(path));

    if (!options.folder_path.empty() && !std::filesystem::exists(options.folder_path))
    {
        std::error_code error_code;
        std::filesystem::create_directories(options.folder_path, error_code);
        if (error_code)
        {
            Logger::Log(LogLevel::Error,
                "Failed to create directory: " + options.folder_path.string() +
                " error code = (" + error_code.message() + ")");
            m_valiate_options = false;
            return;
        }
    }
}

} // namespace rhbm_gem
