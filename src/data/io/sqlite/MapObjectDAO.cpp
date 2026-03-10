#include "internal/io/sqlite/MapObjectDAO.hpp"
#include "internal/io/sqlite/SQLiteWrapper.hpp"
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/DataObjectBase.hpp>
#include <rhbm_gem/data/dispatch/DataObjectDispatch.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include "io/sqlite/MapStoreSql.hpp"

#include <stdexcept>
#include <vector>
#include <cstring>

namespace rhbm_gem {

MapObjectDAO::MapObjectDAO(SQLiteWrapper * database) : m_database{ database }
{
}

MapObjectDAO::~MapObjectDAO()
{
}

void MapObjectDAO::EnsureSchema(SQLiteWrapper & database)
{
    database.Execute(std::string(persistence::kCreateMapTableSql));
}

void MapObjectDAO::Save(const DataObjectBase & data_object, const std::string & key_tag)
{
    const auto & map_object{ ExpectMapObject(data_object, "MapObjectDAO::Save()") };

    m_database->Prepare(std::string(persistence::kInsertMapSql));
    SQLiteWrapper::StatementGuard guard(*m_database);

    auto grid_size{ map_object.GetGridSize() };
    auto grid_spacing{ map_object.GetGridSpacing() };
    auto origin{ map_object.GetOrigin() };
    std::vector<float> values(map_object.GetMapValueArraySize());
    std::memcpy(values.data(), map_object.GetMapValueArray(), values.size() * sizeof(float));

    m_database->Bind<std::string>(1, key_tag);
    m_database->Bind<int>(2, grid_size[0]);
    m_database->Bind<int>(3, grid_size[1]);
    m_database->Bind<int>(4, grid_size[2]);
    m_database->Bind<double>(5, grid_spacing[0]);
    m_database->Bind<double>(6, grid_spacing[1]);
    m_database->Bind<double>(7, grid_spacing[2]);
    m_database->Bind<double>(8, origin[0]);
    m_database->Bind<double>(9, origin[1]);
    m_database->Bind<double>(10, origin[2]);
    m_database->Bind<std::vector<float>>(11, values);

    m_database->StepOnce();
}

std::unique_ptr<DataObjectBase> MapObjectDAO::Load(const std::string & key_tag)
{
    m_database->Prepare(std::string(persistence::kSelectMapSql));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);

    auto rc{ m_database->StepNext() };
    if (rc == SQLiteWrapper::StepDone())
    {
        throw std::runtime_error("Cannot find the row with key_tag = " + key_tag);
    }
    else if (rc != SQLiteWrapper::StepRow())
    {
        throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
    }

    std::array<int,3> grid_size{
        m_database->GetColumn<int>(1),
        m_database->GetColumn<int>(2),
        m_database->GetColumn<int>(3)
    };
    std::array<float,3> grid_spacing{
        static_cast<float>(m_database->GetColumn<double>(4)),
        static_cast<float>(m_database->GetColumn<double>(5)),
        static_cast<float>(m_database->GetColumn<double>(6))
    };
    std::array<float,3> origin{
        static_cast<float>(m_database->GetColumn<double>(7)),
        static_cast<float>(m_database->GetColumn<double>(8)),
        static_cast<float>(m_database->GetColumn<double>(9))
    };
    auto values{ m_database->GetColumn<std::vector<float>>(10) };
    std::unique_ptr<float[]> data_array{ std::make_unique<float[]>(values.size()) };
    std::memcpy(data_array.get(), values.data(), values.size() * sizeof(float));

    auto map_object{
        std::make_unique<MapObject>(grid_size, grid_spacing, origin, std::move(data_array))
    };
    map_object->SetKeyTag(m_database->GetColumn<std::string>(0));
    return map_object;
}

} // namespace rhbm_gem
