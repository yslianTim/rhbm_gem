#pragma once

#include <memory>
#include <string>

#include "DataObjectDAOBase.hpp"

class SQLiteWrapper;
class MapObject;

class MapObjectDAO : public DataObjectDAOBase
{
    SQLiteWrapper * m_database;

public:
    explicit MapObjectDAO(SQLiteWrapper * database);
    ~MapObjectDAO();
    void Save(const DataObjectBase * data_object, const std::string & key_tag) override;
    std::unique_ptr<DataObjectBase> Load(const std::string & key_tag) override;

};
