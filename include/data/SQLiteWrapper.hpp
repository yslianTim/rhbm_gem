#pragma once

#include <string>
#include <vector>
#include <tuple>
#include <utility>
#include "SQLiteBinder.hpp"
#include "SQLiteColumnReader.hpp"

class SQLiteWrapper
{
    sqlite3 * m_database_ptr;
    sqlite3_stmt * m_statement_ptr;

public:
    explicit SQLiteWrapper(const std::string & database_path);
    ~SQLiteWrapper();

    SQLiteWrapper(const SQLiteWrapper &) = delete;
    SQLiteWrapper & operator=(const SQLiteWrapper &) = delete;

    static int StepDone(void);
    static int StepRow(void);

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
                Finalize();
                throw std::runtime_error("Query step failed: " + ErrorMessage());
            }
        }
        Finalize();
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
