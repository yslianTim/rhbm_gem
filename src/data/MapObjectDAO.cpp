#include "MapObjectDAO.hpp"
#include "SQLiteWrapper.hpp"
#include "MapObject.hpp"
#include "DataObjectBase.hpp"
#include "DataObjectDAOFactoryRegistry.hpp"
#include "Logger.hpp"

#include <stdexcept>
#include <vector>
#include <cstring>

namespace { DataObjectDAORegistrar<MapObject, MapObjectDAO> registrar_map_dao("map"); }

namespace
{
    constexpr std::string_view CREATE_TABLE_SQL = R"sql(
        CREATE TABLE IF NOT EXISTS map_list (
            key_tag TEXT PRIMARY KEY,
            grid_size_x INTEGER,
            grid_size_y INTEGER,
            grid_size_z INTEGER,
            grid_spacing_x DOUBLE,
            grid_spacing_y DOUBLE,
            grid_spacing_z DOUBLE,
            origin_x DOUBLE,
            origin_y DOUBLE,
            origin_z DOUBLE,
            map_value_array BLOB
        ) )sql";

    constexpr std::string_view INSERT_SQL = R"sql(
        INSERT INTO map_list (
            key_tag,
            grid_size_x, grid_size_y, grid_size_z,
            grid_spacing_x, grid_spacing_y, grid_spacing_z,
            origin_x, origin_y, origin_z,
            map_value_array
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(key_tag) DO UPDATE SET
            key_tag = excluded.key_tag,
            grid_size_x = excluded.grid_size_x,
            grid_size_y = excluded.grid_size_y,
            grid_size_z = excluded.grid_size_z,
            grid_spacing_x = excluded.grid_spacing_x,
            grid_spacing_y = excluded.grid_spacing_y,
            grid_spacing_z = excluded.grid_spacing_z,
            origin_x = excluded.origin_x,
            origin_y = excluded.origin_y,
            origin_z = excluded.origin_z,
            map_value_array = excluded.map_value_array
        )sql";

    constexpr std::string_view SELECT_SQL = R"sql(
        SELECT
            key_tag,
            grid_size_x, grid_size_y, grid_size_z,
            grid_spacing_x, grid_spacing_y, grid_spacing_z,
            origin_x, origin_y, origin_z,
            map_value_array
        FROM map_list WHERE key_tag = ? LIMIT 1; )sql";
}

MapObjectDAO::MapObjectDAO(SQLiteWrapper * database) : m_database{ database }
{
    Logger::Log(LogLevel::Debug, "MapObjectDAO::MapObjectDAO() called");
}

MapObjectDAO::~MapObjectDAO()
{
    Logger::Log(LogLevel::Debug, "MapObjectDAO::~MapObjectDAO() called");
}

void MapObjectDAO::Save(const DataObjectBase * data_object, const std::string & key_tag)
{
    Logger::Log(LogLevel::Debug, "MapObjectDAO::Save() called");
    auto map_object{ dynamic_cast<const MapObject *>(data_object) };
    if (!map_object)
    {
        throw std::runtime_error("MapObjectDAO::Save() failed: object is not a MapObject instance.");
    }

    SQLiteWrapper::TransactionGuard transaction(*m_database);
    m_database->Execute(std::string(CREATE_TABLE_SQL));
    m_database->Prepare(std::string(INSERT_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);

    auto grid_size{ map_object->GetGridSize() };
    auto grid_spacing{ map_object->GetGridSpacing() };
    auto origin{ map_object->GetOrigin() };
    std::vector<float> values(map_object->GetMapValueArraySize());
    std::memcpy(values.data(), map_object->GetMapValueArray(), values.size() * sizeof(float));

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
    Logger::Log(LogLevel::Debug, "MapObjectDAO::Load() called");
    m_database->Prepare(std::string(SELECT_SQL));
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
