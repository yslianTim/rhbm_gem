#pragma once

#include <array>
#include <string_view>

namespace rhbm_gem::persistence {

inline constexpr std::string_view kMapTableName = "map_list";

inline constexpr std::string_view kCreateMapTableSql = R"sql(
    CREATE TABLE IF NOT EXISTS map_list (
        key_tag TEXT PRIMARY KEY REFERENCES object_catalog(key_tag) ON DELETE CASCADE,
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
    )
)sql";

inline constexpr std::string_view kInsertMapSql = R"sql(
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

inline constexpr std::string_view kSelectMapSql = R"sql(
    SELECT
        key_tag,
        grid_size_x, grid_size_y, grid_size_z,
        grid_spacing_x, grid_spacing_y, grid_spacing_z,
        origin_x, origin_y, origin_z,
        map_value_array
    FROM map_list WHERE key_tag = ? LIMIT 1;
)sql";

inline constexpr std::array<std::string_view, 1> kMapCanonicalTableNames{
    kMapTableName
};

} // namespace rhbm_gem::persistence
