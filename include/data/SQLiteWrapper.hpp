#pragma once

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
            m_db.BeginTransaction();
        }

        ~TransactionGuard()
        {
            if (std::uncaught_exceptions() > m_exception_count)
            {
                m_db.RollbackTransaction();
            }
            else
            {
                m_db.CommitTransaction();
            }
        }

        TransactionGuard(const TransactionGuard &) = delete;
        TransactionGuard & operator=(const TransactionGuard &) = delete;
    };

public:
    explicit SQLiteWrapper(const std::filesystem::path & database_path);
    ~SQLiteWrapper();

    SQLiteWrapper(const SQLiteWrapper &) = delete;
    SQLiteWrapper & operator=(const SQLiteWrapper &) = delete;

    static int StepDone(void) { return SQLITE_DONE; }
    static int StepRow(void) { return SQLITE_ROW; }

    void Execute(const std::string & sql);
    void Prepare(const std::string & sql);
    int StepNext(void);
    void StepOnce(void);
    void Reset(void);
    void Finalize(void);
    std::string ErrorMessage(void) const;
    void BeginTransaction(void);
    void CommitTransaction(void);
    void RollbackTransaction(void);
    void ClearTable(const std::string & table_name);

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
