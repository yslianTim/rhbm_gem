#include "SQLiteWrapper.hpp"

#include <stdexcept>

using std::vector;
using std::tuple;
using std::string;

SQLiteWrapper::SQLiteWrapper(const std::filesystem::path & database_path) :
    m_database_ptr{ nullptr }, m_statement_ptr{ nullptr }
{
    auto utf8_path{ database_path.u8string() };
    auto return_code{ sqlite3_open(utf8_path.c_str(), &m_database_ptr) };
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

int SQLiteWrapper::StepDone(void) { return SQLITE_DONE; }
int SQLiteWrapper::StepRow(void) { return SQLITE_ROW; }

void SQLiteWrapper::Execute(const string & sql)
{
    char * error_message{ nullptr };
    auto return_code{ sqlite3_exec(m_database_ptr, sql.c_str(), nullptr, nullptr, &error_message) };
    if (return_code != SQLITE_OK)
    {
        string msg{ (error_message ? error_message : "Unknown error") };
        sqlite3_free(error_message);
        throw std::runtime_error("Execute SQL failed: " + msg);
    }
}

void SQLiteWrapper::Prepare(const string & sql)
{
    Finalize();
    auto return_code{ sqlite3_prepare_v2(m_database_ptr, sql.c_str(), -1, &m_statement_ptr, nullptr) };
    if (return_code != SQLITE_OK)
    {
        throw std::runtime_error("Prepared failed: " + ErrorMessage());
    }
}

int SQLiteWrapper::StepNext(void)
{
    return sqlite3_step(m_statement_ptr);
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

void SQLiteWrapper::Reset(void)
{
    sqlite3_reset(m_statement_ptr);
}

void SQLiteWrapper::Finalize(void)
{
    if (m_statement_ptr)
    {
        sqlite3_finalize(m_statement_ptr);
    }
    m_statement_ptr = nullptr;
}

string SQLiteWrapper::ErrorMessage(void) const
{
    return string(sqlite3_errmsg(m_database_ptr));
}

void SQLiteWrapper::BeginTransaction(void)
{
    this->Execute("BEGIN TRANSACTION;");
}

void SQLiteWrapper::CommitTransaction(void)
{
    this->Execute("COMMIT;");
}

void SQLiteWrapper::RollbackTransaction(void)
{
    this->Execute("ROLLBACK;");
}

void SQLiteWrapper::ClearTable(const std::string & table_name)
{
    std::string sql{ "DELETE FROM " + table_name + ";" };
    this->Execute(sql);
}
