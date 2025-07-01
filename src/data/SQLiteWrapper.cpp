#include "SQLiteWrapper.hpp"

#include <stdexcept>

SQLiteWrapper::SQLiteWrapper(const std::filesystem::path & database_path) :
    m_database_ptr{ nullptr }, m_statement_ptr{ nullptr }
{
#ifdef _WIN32
    std::wstring wpath{ database_path.wstring() };
    auto return_code{ sqlite3_open16(wpath.c_str(), &m_database_ptr) };
#else
    std::string path_utf8{ database_path.u8string() };
    auto return_code{ sqlite3_open_v2(path_utf8.c_str(), &m_database_ptr,
                                      SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE,
                                      nullptr) };
#endif
    if (return_code != SQLITE_OK)
    {
        sqlite3_close(m_database_ptr);
        m_database_ptr = nullptr;
        throw std::runtime_error("Cannot open database: " + ErrorMessage());
    }
}

SQLiteWrapper::~SQLiteWrapper()
{
    if (m_statement_ptr)
    {
        Finalize();
    }
    if (m_database_ptr)
    {
        sqlite3_close(m_database_ptr);
        m_database_ptr = nullptr;
    }
}

void SQLiteWrapper::Execute(const std::string & sql)
{
    char * error_message{ nullptr };
    auto return_code{ sqlite3_exec(m_database_ptr, sql.c_str(), nullptr, nullptr, &error_message) };
    if (return_code != SQLITE_OK)
    {
        std::string msg{ (error_message ? error_message : "Unknown error") };
        sqlite3_free(error_message);
        throw std::runtime_error("Execute SQL failed: " + msg);
    }
}

void SQLiteWrapper::Prepare(const std::string & sql)
{
    Finalize();
    auto return_code{ sqlite3_prepare_v2(m_database_ptr, sql.c_str(), -1, &m_statement_ptr, nullptr) };
    if (return_code != SQLITE_OK)
    {
        throw std::runtime_error("Prepared failed: " + ErrorMessage());
    }
}

void SQLiteWrapper::StepOnce(void)
{
    auto return_code{ sqlite3_step(m_statement_ptr) };
    if (return_code != SQLITE_DONE)
    {
        Finalize();
        throw std::runtime_error("Statement step failed: " + ErrorMessage());
    }
}

void SQLiteWrapper::Finalize(void)
{
    if (m_statement_ptr)
    {
        sqlite3_finalize(m_statement_ptr);
    }
    m_statement_ptr = nullptr;
}

void SQLiteWrapper::ClearTable(const std::string & table_name)
{
    std::string sql{ "DELETE FROM " + table_name + ";" };
    this->Execute(sql);
}
