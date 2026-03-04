#pragma once

#include <array>
#include <string_view>

namespace rhbm_gem::model_io {

inline constexpr std::string_view kCreateModelObjectTableSql = R"sql(
    CREATE TABLE IF NOT EXISTS model_object (
        key_tag TEXT PRIMARY KEY,
        atom_size INTEGER,
        pdb_id TEXT,
        emd_id TEXT,
        map_resolution DOUBLE,
        resolution_method TEXT
    )
)sql";

inline constexpr std::string_view kCreateModelChainMapTableSql = R"sql(
    CREATE TABLE IF NOT EXISTS model_chain_map (
        key_tag TEXT,
        entity_id TEXT,
        chain_ordinal INTEGER,
        chain_id TEXT,
        PRIMARY KEY (key_tag, entity_id, chain_ordinal)
    )
)sql";

inline constexpr std::string_view kCreateModelComponentTableSql = R"sql(
    CREATE TABLE IF NOT EXISTS model_component (
        key_tag TEXT,
        component_key INTEGER,
        id TEXT,
        name TEXT,
        type TEXT,
        formula TEXT,
        molecular_weight DOUBLE,
        is_standard_monomer INTEGER,
        PRIMARY KEY (key_tag, component_key)
    )
)sql";

inline constexpr std::string_view kCreateModelComponentAtomTableSql = R"sql(
    CREATE TABLE IF NOT EXISTS model_component_atom (
        key_tag TEXT,
        component_key INTEGER,
        atom_key INTEGER,
        atom_id TEXT,
        element_type INTEGER,
        aromatic_atom_flag INTEGER,
        stereo_config INTEGER,
        PRIMARY KEY (key_tag, component_key, atom_key)
    )
)sql";

inline constexpr std::string_view kCreateModelComponentBondTableSql = R"sql(
    CREATE TABLE IF NOT EXISTS model_component_bond (
        key_tag TEXT,
        component_key INTEGER,
        bond_key INTEGER,
        bond_id TEXT,
        bond_type INTEGER,
        bond_order INTEGER,
        aromatic_atom_flag INTEGER,
        stereo_config INTEGER,
        PRIMARY KEY (key_tag, component_key, bond_key)
    )
)sql";

inline constexpr std::string_view kCreateModelAtomTableSql = R"sql(
    CREATE TABLE IF NOT EXISTS model_atom (
        key_tag TEXT,
        serial_id INTEGER,
        sequence_id INTEGER,
        component_id TEXT,
        atom_id TEXT,
        chain_id TEXT,
        indicator TEXT,
        occupancy DOUBLE,
        temperature DOUBLE,
        element INTEGER,
        structure INTEGER,
        is_special_atom INTEGER,
        position_x DOUBLE,
        position_y DOUBLE,
        position_z DOUBLE,
        component_key INTEGER,
        atom_key INTEGER,
        PRIMARY KEY (key_tag, serial_id)
    )
)sql";

inline constexpr std::string_view kCreateModelBondTableSql = R"sql(
    CREATE TABLE IF NOT EXISTS model_bond (
        key_tag TEXT,
        atom_serial_id_1 INTEGER,
        atom_serial_id_2 INTEGER,
        bond_key INTEGER,
        bond_type INTEGER,
        bond_order INTEGER,
        is_special_bond INTEGER,
        PRIMARY KEY (key_tag, atom_serial_id_1, atom_serial_id_2)
    )
)sql";

inline constexpr std::string_view kCreateModelAtomLocalTableSql = R"sql(
    CREATE TABLE IF NOT EXISTS model_atom_local_potential (
        key_tag TEXT,
        serial_id INTEGER,
        sampling_size INTEGER,
        distance_and_map_value_list BLOB,
        amplitude_estimate_ols DOUBLE,
        width_estimate_ols DOUBLE,
        amplitude_estimate_mdpde DOUBLE,
        width_estimate_mdpde DOUBLE,
        alpha_r DOUBLE,
        PRIMARY KEY (key_tag, serial_id)
    )
)sql";

inline constexpr std::string_view kCreateModelBondLocalTableSql = R"sql(
    CREATE TABLE IF NOT EXISTS model_bond_local_potential (
        key_tag TEXT,
        atom_serial_id_1 INTEGER,
        atom_serial_id_2 INTEGER,
        sampling_size INTEGER,
        distance_and_map_value_list BLOB,
        amplitude_estimate_ols DOUBLE,
        width_estimate_ols DOUBLE,
        amplitude_estimate_mdpde DOUBLE,
        width_estimate_mdpde DOUBLE,
        alpha_r DOUBLE,
        PRIMARY KEY (key_tag, atom_serial_id_1, atom_serial_id_2)
    )
)sql";

inline constexpr std::string_view kCreateModelAtomPosteriorTableSql = R"sql(
    CREATE TABLE IF NOT EXISTS model_atom_posterior (
        key_tag TEXT,
        class_key TEXT,
        serial_id INTEGER,
        amplitude_estimate_posterior DOUBLE,
        width_estimate_posterior DOUBLE,
        amplitude_variance_posterior DOUBLE,
        width_variance_posterior DOUBLE,
        outlier_tag INTEGER,
        statistical_distance DOUBLE,
        PRIMARY KEY (key_tag, class_key, serial_id)
    )
)sql";

inline constexpr std::string_view kCreateModelBondPosteriorTableSql = R"sql(
    CREATE TABLE IF NOT EXISTS model_bond_posterior (
        key_tag TEXT,
        class_key TEXT,
        atom_serial_id_1 INTEGER,
        atom_serial_id_2 INTEGER,
        amplitude_estimate_posterior DOUBLE,
        width_estimate_posterior DOUBLE,
        amplitude_variance_posterior DOUBLE,
        width_variance_posterior DOUBLE,
        outlier_tag INTEGER,
        statistical_distance DOUBLE,
        PRIMARY KEY (key_tag, class_key, atom_serial_id_1, atom_serial_id_2)
    )
)sql";

inline constexpr std::string_view kCreateModelAtomGroupTableSql = R"sql(
    CREATE TABLE IF NOT EXISTS model_atom_group_potential (
        key_tag TEXT,
        class_key TEXT,
        group_key INTEGER,
        member_size INTEGER,
        amplitude_estimate_mean DOUBLE,
        width_estimate_mean DOUBLE,
        amplitude_estimate_mdpde DOUBLE,
        width_estimate_mdpde DOUBLE,
        amplitude_estimate_prior DOUBLE,
        width_estimate_prior DOUBLE,
        amplitude_variance_prior DOUBLE,
        width_variance_prior DOUBLE,
        alpha_g DOUBLE,
        PRIMARY KEY (key_tag, class_key, group_key)
    )
)sql";

inline constexpr std::string_view kCreateModelBondGroupTableSql = R"sql(
    CREATE TABLE IF NOT EXISTS model_bond_group_potential (
        key_tag TEXT,
        class_key TEXT,
        group_key INTEGER,
        member_size INTEGER,
        amplitude_estimate_mean DOUBLE,
        width_estimate_mean DOUBLE,
        amplitude_estimate_mdpde DOUBLE,
        width_estimate_mdpde DOUBLE,
        amplitude_estimate_prior DOUBLE,
        width_estimate_prior DOUBLE,
        amplitude_variance_prior DOUBLE,
        width_variance_prior DOUBLE,
        alpha_g DOUBLE,
        PRIMARY KEY (key_tag, class_key, group_key)
    )
)sql";

inline constexpr std::string_view kUpsertModelObjectSql = R"sql(
    INSERT INTO model_object (
        key_tag, atom_size, pdb_id, emd_id, map_resolution, resolution_method
    ) VALUES (?, ?, ?, ?, ?, ?)
    ON CONFLICT(key_tag) DO UPDATE SET
        atom_size = excluded.atom_size,
        pdb_id = excluded.pdb_id,
        emd_id = excluded.emd_id,
        map_resolution = excluded.map_resolution,
        resolution_method = excluded.resolution_method
)sql";

inline constexpr std::string_view kInsertModelChainMapSql = R"sql(
    INSERT OR REPLACE INTO model_chain_map (
        key_tag, entity_id, chain_ordinal, chain_id
    ) VALUES (?, ?, ?, ?)
)sql";

inline constexpr std::string_view kInsertModelComponentSql = R"sql(
    INSERT OR REPLACE INTO model_component (
        key_tag, component_key, id, name, type, formula,
        molecular_weight, is_standard_monomer
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
)sql";

inline constexpr std::string_view kInsertModelComponentAtomSql = R"sql(
    INSERT OR REPLACE INTO model_component_atom (
        key_tag, component_key, atom_key, atom_id,
        element_type, aromatic_atom_flag, stereo_config
    ) VALUES (?, ?, ?, ?, ?, ?, ?)
)sql";

inline constexpr std::string_view kInsertModelComponentBondSql = R"sql(
    INSERT OR REPLACE INTO model_component_bond (
        key_tag, component_key, bond_key, bond_id,
        bond_type, bond_order, aromatic_atom_flag, stereo_config
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
)sql";

inline constexpr std::string_view kInsertModelAtomSql = R"sql(
    INSERT OR REPLACE INTO model_atom (
        key_tag, serial_id, sequence_id, component_id, atom_id, chain_id, indicator,
        occupancy, temperature, element, structure, is_special_atom,
        position_x, position_y, position_z, component_key, atom_key
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
)sql";

inline constexpr std::string_view kInsertModelBondSql = R"sql(
    INSERT OR REPLACE INTO model_bond (
        key_tag, atom_serial_id_1, atom_serial_id_2,
        bond_key, bond_type, bond_order, is_special_bond
    ) VALUES (?, ?, ?, ?, ?, ?, ?)
)sql";

inline constexpr std::string_view kInsertModelAtomLocalSql = R"sql(
    INSERT OR REPLACE INTO model_atom_local_potential (
        key_tag, serial_id, sampling_size, distance_and_map_value_list,
        amplitude_estimate_ols, width_estimate_ols,
        amplitude_estimate_mdpde, width_estimate_mdpde, alpha_r
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
)sql";

inline constexpr std::string_view kInsertModelBondLocalSql = R"sql(
    INSERT OR REPLACE INTO model_bond_local_potential (
        key_tag, atom_serial_id_1, atom_serial_id_2, sampling_size,
        distance_and_map_value_list, amplitude_estimate_ols, width_estimate_ols,
        amplitude_estimate_mdpde, width_estimate_mdpde, alpha_r
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
)sql";

inline constexpr std::string_view kInsertModelAtomPosteriorSql = R"sql(
    INSERT OR REPLACE INTO model_atom_posterior (
        key_tag, class_key, serial_id, amplitude_estimate_posterior, width_estimate_posterior,
        amplitude_variance_posterior, width_variance_posterior, outlier_tag, statistical_distance
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
)sql";

inline constexpr std::string_view kInsertModelBondPosteriorSql = R"sql(
    INSERT OR REPLACE INTO model_bond_posterior (
        key_tag, class_key, atom_serial_id_1, atom_serial_id_2,
        amplitude_estimate_posterior, width_estimate_posterior,
        amplitude_variance_posterior, width_variance_posterior, outlier_tag, statistical_distance
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
)sql";

inline constexpr std::string_view kInsertModelAtomGroupSql = R"sql(
    INSERT OR REPLACE INTO model_atom_group_potential (
        key_tag, class_key, group_key, member_size,
        amplitude_estimate_mean, width_estimate_mean,
        amplitude_estimate_mdpde, width_estimate_mdpde,
        amplitude_estimate_prior, width_estimate_prior,
        amplitude_variance_prior, width_variance_prior, alpha_g
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
)sql";

inline constexpr std::string_view kInsertModelBondGroupSql = R"sql(
    INSERT OR REPLACE INTO model_bond_group_potential (
        key_tag, class_key, group_key, member_size,
        amplitude_estimate_mean, width_estimate_mean,
        amplitude_estimate_mdpde, width_estimate_mdpde,
        amplitude_estimate_prior, width_estimate_prior,
        amplitude_variance_prior, width_variance_prior, alpha_g
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
)sql";

inline constexpr std::string_view kDeleteRowsForKeySqlPrefix = "DELETE FROM ";
inline constexpr std::string_view kDeleteRowsForKeySqlSuffix = " WHERE key_tag = ?;";

inline constexpr std::string_view kSelectModelObjectSql = R"sql(
    SELECT key_tag, atom_size, pdb_id, emd_id, map_resolution, resolution_method
    FROM model_object WHERE key_tag = ? LIMIT 1;
)sql";

inline constexpr std::string_view kSelectModelChainMapSql = R"sql(
    SELECT entity_id, chain_ordinal, chain_id
    FROM model_chain_map WHERE key_tag = ?
    ORDER BY entity_id, chain_ordinal;
)sql";

inline constexpr std::string_view kSelectModelComponentSql = R"sql(
    SELECT component_key, id, name, type, formula, molecular_weight, is_standard_monomer
    FROM model_component WHERE key_tag = ?
    ORDER BY component_key;
)sql";

inline constexpr std::string_view kSelectModelComponentAtomSql = R"sql(
    SELECT component_key, atom_key, atom_id, element_type, aromatic_atom_flag, stereo_config
    FROM model_component_atom WHERE key_tag = ?
    ORDER BY component_key, atom_key;
)sql";

inline constexpr std::string_view kSelectModelComponentBondSql = R"sql(
    SELECT component_key, bond_key, bond_id, bond_type, bond_order, aromatic_atom_flag, stereo_config
    FROM model_component_bond WHERE key_tag = ?
    ORDER BY component_key, bond_key;
)sql";

inline constexpr std::string_view kSelectModelAtomSql = R"sql(
    SELECT
        serial_id, sequence_id, component_id, atom_id, chain_id, indicator,
        occupancy, temperature, element, structure, is_special_atom,
        position_x, position_y, position_z, component_key, atom_key
    FROM model_atom WHERE key_tag = ?
    ORDER BY serial_id;
)sql";

inline constexpr std::string_view kSelectModelBondSql = R"sql(
    SELECT atom_serial_id_1, atom_serial_id_2, bond_key, bond_type, bond_order, is_special_bond
    FROM model_bond WHERE key_tag = ?
    ORDER BY atom_serial_id_1, atom_serial_id_2;
)sql";

inline constexpr std::string_view kSelectModelAtomLocalSql = R"sql(
    SELECT
        serial_id, sampling_size, distance_and_map_value_list,
        amplitude_estimate_ols, width_estimate_ols,
        amplitude_estimate_mdpde, width_estimate_mdpde, alpha_r
    FROM model_atom_local_potential WHERE key_tag = ?;
)sql";

inline constexpr std::string_view kSelectModelBondLocalSql = R"sql(
    SELECT
        atom_serial_id_1, atom_serial_id_2, sampling_size, distance_and_map_value_list,
        amplitude_estimate_ols, width_estimate_ols,
        amplitude_estimate_mdpde, width_estimate_mdpde, alpha_r
    FROM model_bond_local_potential WHERE key_tag = ?;
)sql";

inline constexpr std::string_view kSelectModelAtomPosteriorSql = R"sql(
    SELECT
        serial_id, amplitude_estimate_posterior, width_estimate_posterior,
        amplitude_variance_posterior, width_variance_posterior, outlier_tag, statistical_distance
    FROM model_atom_posterior WHERE key_tag = ? AND class_key = ?;
)sql";

inline constexpr std::string_view kSelectModelBondPosteriorSql = R"sql(
    SELECT
        atom_serial_id_1, atom_serial_id_2, amplitude_estimate_posterior, width_estimate_posterior,
        amplitude_variance_posterior, width_variance_posterior, outlier_tag, statistical_distance
    FROM model_bond_posterior WHERE key_tag = ? AND class_key = ?;
)sql";

inline constexpr std::string_view kSelectModelAtomGroupSql = R"sql(
    SELECT
        group_key, member_size, amplitude_estimate_mean, width_estimate_mean,
        amplitude_estimate_mdpde, width_estimate_mdpde,
        amplitude_estimate_prior, width_estimate_prior,
        amplitude_variance_prior, width_variance_prior, alpha_g
    FROM model_atom_group_potential WHERE key_tag = ? AND class_key = ?;
)sql";

inline constexpr std::string_view kSelectModelBondGroupSql = R"sql(
    SELECT
        group_key, member_size, amplitude_estimate_mean, width_estimate_mean,
        amplitude_estimate_mdpde, width_estimate_mdpde,
        amplitude_estimate_prior, width_estimate_prior,
        amplitude_variance_prior, width_variance_prior, alpha_g
    FROM model_bond_group_potential WHERE key_tag = ? AND class_key = ?;
)sql";

inline constexpr std::array<std::string_view, 13> kCreateModelTableSqlList{
    kCreateModelObjectTableSql,
    kCreateModelChainMapTableSql,
    kCreateModelComponentTableSql,
    kCreateModelComponentAtomTableSql,
    kCreateModelComponentBondTableSql,
    kCreateModelAtomTableSql,
    kCreateModelBondTableSql,
    kCreateModelAtomLocalTableSql,
    kCreateModelBondLocalTableSql,
    kCreateModelAtomPosteriorTableSql,
    kCreateModelBondPosteriorTableSql,
    kCreateModelAtomGroupTableSql,
    kCreateModelBondGroupTableSql
};

inline constexpr std::array<std::string_view, 12> kModelTablesScopedByKey{
    "model_chain_map",
    "model_component",
    "model_component_atom",
    "model_component_bond",
    "model_atom",
    "model_bond",
    "model_atom_local_potential",
    "model_bond_local_potential",
    "model_atom_posterior",
    "model_bond_posterior",
    "model_atom_group_potential",
    "model_bond_group_potential"
};

} // namespace rhbm_gem::model_io
