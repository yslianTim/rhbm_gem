#pragma once

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>
#include <tuple>
#include <utility>
#include <memory>
#include <exception>
#include <filesystem>

#include "SQLiteBinder.hpp"
#include "SQLiteColumnReader.hpp"
#include "Logger.hpp"

namespace rhbm_gem {

class SQLiteWrapper
{
    sqlite3 * m_database_ptr;
    sqlite3_stmt * m_statement_ptr;

public:
    /**
     * @brief Helper class to ensure prepared statements are finalized.
     */
    class StatementGuard
    {
        SQLiteWrapper & m_db;

    public:
        explicit StatementGuard(SQLiteWrapper & db) : m_db{ db } {}
        ~StatementGuard() { m_db.Finalize(); }

        StatementGuard(const StatementGuard &) = delete;
        StatementGuard & operator=(const StatementGuard &) = delete;
    };

    /**
     * @brief Helper class to manage transactions using RAII.
     */
    class TransactionGuard
    {
        SQLiteWrapper & m_db;
        int m_exception_count;

    public:
        explicit TransactionGuard(SQLiteWrapper & db) :
            m_db{ db }, m_exception_count{ std::uncaught_exceptions() }
        {
            m_db.Execute("BEGIN TRANSACTION;");
        }

        ~TransactionGuard() noexcept
        {
            try
            {
                if (std::uncaught_exceptions() > m_exception_count)
                {
                    m_db.Execute("ROLLBACK;");
                }
                else
                {
                    m_db.Execute("COMMIT;");
                }
            }
            catch (const std::exception & e)
            {
                Logger::Log(LogLevel::Error,
                            "TransactionGuard: Failed to commit or rollback transaction: "
                            + std::string(e.what()));
            }
        }

        TransactionGuard(const TransactionGuard &) = delete;
        TransactionGuard & operator=(const TransactionGuard &) = delete;
    };

    template<typename... Ts>
    class QueryIterator
    {
        SQLiteWrapper & m_db;
        StatementGuard m_guard;

    public:
        QueryIterator(SQLiteWrapper & db, const std::string & sql) : m_db{ db }, m_guard{ db }
        {
            m_db.Prepare(sql);
        }

        bool Next(std::tuple<Ts...> & row)
        {
            auto rc{ m_db.StepNext() };
            if (rc == StepRow())
            {
                row = m_db.ReadRow<Ts...>(m_db, std::make_index_sequence<sizeof...(Ts)>{});
                return true;
            }
            else if (rc == StepDone())
            {
                return false;
            }
            else
            {
                throw std::runtime_error("Query step failed: " + m_db.ErrorMessage());
            }
        }

        QueryIterator(const QueryIterator &) = delete;
        QueryIterator & operator=(const QueryIterator &) = delete;
    };

public:
    explicit SQLiteWrapper(const std::filesystem::path & database_path) :
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

    ~SQLiteWrapper()
    {
        if (m_statement_ptr)
        {
            Finalize();
        }
        if (m_database_ptr)
        {
            int rc{ sqlite3_close(m_database_ptr) };
            if (rc == SQLITE_BUSY)
            {
                Finalize();
                rc = sqlite3_close(m_database_ptr);
            }
            if (rc != SQLITE_OK)
            {
                Logger::Log(LogLevel::Error,
                            "SQLiteWrapper: Failed to close database: "
                            + std::string(sqlite3_errstr(rc)));
            }
            m_database_ptr = nullptr;
        }
    }

    SQLiteWrapper(const SQLiteWrapper &) = delete;
    SQLiteWrapper & operator=(const SQLiteWrapper &) = delete;

    static int StepDone() { return SQLITE_DONE; }
    static int StepRow() { return SQLITE_ROW; }
    std::string ErrorMessage() const { return std::string(sqlite3_errmsg(m_database_ptr)); }

    void Reset()
    {
        if (m_statement_ptr == nullptr)
        {
            throw std::runtime_error("Reset called but statement pointer is null!");
        }
        sqlite3_reset(m_statement_ptr);
        sqlite3_clear_bindings(m_statement_ptr);
    }

    int StepNext()
    {
        if (m_statement_ptr == nullptr)
        {
            throw std::runtime_error("StepNext called but statement pointer is null!");
        }
        return sqlite3_step(m_statement_ptr);
    }

    void Execute(const std::string & sql)
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

    void Prepare(const std::string & sql)
    {
        Finalize();
        auto return_code{ sqlite3_prepare_v2(m_database_ptr, sql.c_str(), -1, &m_statement_ptr, nullptr) };
        if (return_code != SQLITE_OK)
        {
            throw std::runtime_error("Prepared failed: " + ErrorMessage());
        }
    }

    void StepOnce()
    {
        if (m_statement_ptr == nullptr)
        {
            throw std::runtime_error("StepOnce called but statement pointer is null!");
        }
        auto return_code{ sqlite3_step(m_statement_ptr) };
        if (return_code != SQLITE_DONE)
        {
            Finalize();
            throw std::runtime_error("Statement step failed: " + ErrorMessage());
        }
    }
    
    void Finalize()
    {
        if (m_statement_ptr)
        {
            sqlite3_finalize(m_statement_ptr);
        }
        m_statement_ptr = nullptr;
    }

    void ClearTable(const std::string & table_name)
    {
        std::string sql{ "DELETE FROM " + table_name + ";" };
        this->Execute(sql);
    }

    template<typename... Ts>
    std::vector<std::tuple<Ts...>> Query(const std::string & sql)
    {
        std::vector<std::tuple<Ts...>> results;
        Prepare(sql);
        StatementGuard guard(*this);
        while (true)
        {
            auto return_code{ StepNext() };
            if (return_code == StepDone())
            {
                break;
            }
            else if (return_code == StepRow())
            {
                std::tuple<Ts...> row{ ReadRow<Ts...>(*this, std::make_index_sequence<sizeof...(Ts)>{}) };
                results.emplace_back(std::move(row));
            }
            else
            {
                throw std::runtime_error("Query step failed: " + ErrorMessage());
            }
        }
        return results;
    }

    template<typename... Ts>
    QueryIterator<Ts...> IterateQuery(const std::string & sql)
    {
        return QueryIterator<Ts...>(*this, sql);
    }

    template<typename T>
    void Bind(int index, const T & value)
    {
        if (!m_statement_ptr)
        {
            throw std::runtime_error("Statement pointer is null!");
        }
        auto return_code{ SQLiteBinder<T>::Bind(m_statement_ptr, index, value) };
        if (return_code != SQLITE_OK)
        {
            Finalize();
            throw std::runtime_error("Bind failed: " + ErrorMessage());
        }
    }

    template<typename T>
    T GetColumn(int index)
    {
        if (!m_statement_ptr)
        {
            throw std::runtime_error("Statement pointer is null!");
        }
        return SQLiteColumnReader<T>::Get(m_statement_ptr, index);
    }

private:
    template<typename... Ts, std::size_t... Is>
    std::tuple<Ts...> ReadRow(SQLiteWrapper & db, std::index_sequence<Is...>)
    {
        return std::make_tuple(db.GetColumn<Ts>(Is)...);
    }
    
};

} // namespace rhbm_gem
