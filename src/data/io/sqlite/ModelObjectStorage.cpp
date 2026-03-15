#include "internal/io/sqlite/ModelObjectStorage.hpp"

#include "internal/io/sqlite/SQLiteWrapper.hpp"
#include <rhbm_gem/data/object/AtomClassifier.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondClassifier.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ChemicalComponentEntry.hpp>
#include <rhbm_gem/data/object/GroupPotentialEntry.hpp>
#include <rhbm_gem/data/object/LocalPotentialEntry.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/AtomKeySystem.hpp>
#include <rhbm_gem/utils/domain/BondKeySystem.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/ComponentKeySystem.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>

#include <array>
#include <map>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace {

using namespace std::literals;

inline constexpr std::string_view kCreateModelObjectTableSql = R"sql(
    CREATE TABLE IF NOT EXISTS model_object (
        key_tag TEXT PRIMARY KEY REFERENCES object_catalog(key_tag) ON DELETE CASCADE,
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
        PRIMARY KEY (key_tag, entity_id, chain_ordinal),
        FOREIGN KEY(key_tag) REFERENCES model_object(key_tag) ON DELETE CASCADE
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
        PRIMARY KEY (key_tag, component_key),
        FOREIGN KEY(key_tag) REFERENCES model_object(key_tag) ON DELETE CASCADE
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
        PRIMARY KEY (key_tag, component_key, atom_key),
        FOREIGN KEY(key_tag) REFERENCES model_object(key_tag) ON DELETE CASCADE
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
        PRIMARY KEY (key_tag, component_key, bond_key),
        FOREIGN KEY(key_tag) REFERENCES model_object(key_tag) ON DELETE CASCADE
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
        PRIMARY KEY (key_tag, serial_id),
        FOREIGN KEY(key_tag) REFERENCES model_object(key_tag) ON DELETE CASCADE
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
        PRIMARY KEY (key_tag, atom_serial_id_1, atom_serial_id_2),
        FOREIGN KEY(key_tag) REFERENCES model_object(key_tag) ON DELETE CASCADE
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
        PRIMARY KEY (key_tag, serial_id),
        FOREIGN KEY(key_tag) REFERENCES model_object(key_tag) ON DELETE CASCADE
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
        PRIMARY KEY (key_tag, atom_serial_id_1, atom_serial_id_2),
        FOREIGN KEY(key_tag) REFERENCES model_object(key_tag) ON DELETE CASCADE
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
        PRIMARY KEY (key_tag, class_key, serial_id),
        FOREIGN KEY(key_tag) REFERENCES model_object(key_tag) ON DELETE CASCADE
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
        PRIMARY KEY (key_tag, class_key, atom_serial_id_1, atom_serial_id_2),
        FOREIGN KEY(key_tag) REFERENCES model_object(key_tag) ON DELETE CASCADE
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
        PRIMARY KEY (key_tag, class_key, group_key),
        FOREIGN KEY(key_tag) REFERENCES model_object(key_tag) ON DELETE CASCADE
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
        PRIMARY KEY (key_tag, class_key, group_key),
        FOREIGN KEY(key_tag) REFERENCES model_object(key_tag) ON DELETE CASCADE
    )
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

inline constexpr auto kUpsertModelObjectSql = R"sql(
    INSERT INTO model_object (
        key_tag, atom_size, pdb_id, emd_id, map_resolution, resolution_method
    ) VALUES (?, ?, ?, ?, ?, ?)
    ON CONFLICT(key_tag) DO UPDATE SET
        atom_size = excluded.atom_size,
        pdb_id = excluded.pdb_id,
        emd_id = excluded.emd_id,
        map_resolution = excluded.map_resolution,
        resolution_method = excluded.resolution_method
)sql"sv;

inline constexpr auto kInsertModelChainMapSql = R"sql(
    INSERT OR REPLACE INTO model_chain_map (
        key_tag, entity_id, chain_ordinal, chain_id
    ) VALUES (?, ?, ?, ?)
)sql"sv;

inline constexpr auto kInsertModelComponentSql = R"sql(
    INSERT OR REPLACE INTO model_component (
        key_tag, component_key, id, name, type, formula,
        molecular_weight, is_standard_monomer
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
)sql"sv;

inline constexpr auto kInsertModelComponentAtomSql = R"sql(
    INSERT OR REPLACE INTO model_component_atom (
        key_tag, component_key, atom_key, atom_id,
        element_type, aromatic_atom_flag, stereo_config
    ) VALUES (?, ?, ?, ?, ?, ?, ?)
)sql"sv;

inline constexpr auto kInsertModelComponentBondSql = R"sql(
    INSERT OR REPLACE INTO model_component_bond (
        key_tag, component_key, bond_key, bond_id,
        bond_type, bond_order, aromatic_atom_flag, stereo_config
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
)sql"sv;

inline constexpr auto kInsertModelAtomSql = R"sql(
    INSERT OR REPLACE INTO model_atom (
        key_tag, serial_id, sequence_id, component_id, atom_id, chain_id, indicator,
        occupancy, temperature, element, structure, is_special_atom,
        position_x, position_y, position_z, component_key, atom_key
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
)sql"sv;

inline constexpr auto kInsertModelBondSql = R"sql(
    INSERT OR REPLACE INTO model_bond (
        key_tag, atom_serial_id_1, atom_serial_id_2,
        bond_key, bond_type, bond_order, is_special_bond
    ) VALUES (?, ?, ?, ?, ?, ?, ?)
)sql"sv;

inline constexpr auto kInsertModelAtomLocalSql = R"sql(
    INSERT OR REPLACE INTO model_atom_local_potential (
        key_tag, serial_id, sampling_size, distance_and_map_value_list,
        amplitude_estimate_ols, width_estimate_ols,
        amplitude_estimate_mdpde, width_estimate_mdpde, alpha_r
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
)sql"sv;

inline constexpr auto kInsertModelBondLocalSql = R"sql(
    INSERT OR REPLACE INTO model_bond_local_potential (
        key_tag, atom_serial_id_1, atom_serial_id_2, sampling_size,
        distance_and_map_value_list, amplitude_estimate_ols, width_estimate_ols,
        amplitude_estimate_mdpde, width_estimate_mdpde, alpha_r
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
)sql"sv;

inline constexpr auto kInsertModelAtomPosteriorSql = R"sql(
    INSERT OR REPLACE INTO model_atom_posterior (
        key_tag, class_key, serial_id, amplitude_estimate_posterior, width_estimate_posterior,
        amplitude_variance_posterior, width_variance_posterior, outlier_tag, statistical_distance
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
)sql"sv;

inline constexpr auto kInsertModelBondPosteriorSql = R"sql(
    INSERT OR REPLACE INTO model_bond_posterior (
        key_tag, class_key, atom_serial_id_1, atom_serial_id_2,
        amplitude_estimate_posterior, width_estimate_posterior,
        amplitude_variance_posterior, width_variance_posterior, outlier_tag, statistical_distance
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
)sql"sv;

inline constexpr auto kInsertModelAtomGroupSql = R"sql(
    INSERT OR REPLACE INTO model_atom_group_potential (
        key_tag, class_key, group_key, member_size,
        amplitude_estimate_mean, width_estimate_mean,
        amplitude_estimate_mdpde, width_estimate_mdpde,
        amplitude_estimate_prior, width_estimate_prior,
        amplitude_variance_prior, width_variance_prior, alpha_g
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
)sql"sv;

inline constexpr auto kInsertModelBondGroupSql = R"sql(
    INSERT OR REPLACE INTO model_bond_group_potential (
        key_tag, class_key, group_key, member_size,
        amplitude_estimate_mean, width_estimate_mean,
        amplitude_estimate_mdpde, width_estimate_mdpde,
        amplitude_estimate_prior, width_estimate_prior,
        amplitude_variance_prior, width_variance_prior, alpha_g
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
)sql"sv;

inline constexpr auto kDeleteRowsForKeySqlPrefix = "DELETE FROM "sv;
inline constexpr auto kDeleteRowsForKeySqlSuffix = " WHERE key_tag = ?;"sv;

inline constexpr auto kSelectModelObjectSql = R"sql(
    SELECT key_tag, atom_size, pdb_id, emd_id, map_resolution, resolution_method
    FROM model_object WHERE key_tag = ? LIMIT 1;
)sql"sv;

inline constexpr auto kSelectModelChainMapSql = R"sql(
    SELECT entity_id, chain_ordinal, chain_id
    FROM model_chain_map WHERE key_tag = ?
    ORDER BY entity_id, chain_ordinal;
)sql"sv;

inline constexpr auto kSelectModelComponentSql = R"sql(
    SELECT component_key, id, name, type, formula, molecular_weight, is_standard_monomer
    FROM model_component WHERE key_tag = ?
    ORDER BY component_key;
)sql"sv;

inline constexpr auto kSelectModelComponentAtomSql = R"sql(
    SELECT component_key, atom_key, atom_id, element_type, aromatic_atom_flag, stereo_config
    FROM model_component_atom WHERE key_tag = ?
    ORDER BY component_key, atom_key;
)sql"sv;

inline constexpr auto kSelectModelComponentBondSql = R"sql(
    SELECT component_key, bond_key, bond_id, bond_type, bond_order, aromatic_atom_flag, stereo_config
    FROM model_component_bond WHERE key_tag = ?
    ORDER BY component_key, bond_key;
)sql"sv;

inline constexpr auto kSelectModelAtomSql = R"sql(
    SELECT
        serial_id, sequence_id, component_id, atom_id, chain_id, indicator,
        occupancy, temperature, element, structure, is_special_atom,
        position_x, position_y, position_z, component_key, atom_key
    FROM model_atom WHERE key_tag = ?
    ORDER BY serial_id;
)sql"sv;

inline constexpr auto kSelectModelBondSql = R"sql(
    SELECT atom_serial_id_1, atom_serial_id_2, bond_key, bond_type, bond_order, is_special_bond
    FROM model_bond WHERE key_tag = ?
    ORDER BY atom_serial_id_1, atom_serial_id_2;
)sql"sv;

inline constexpr auto kSelectModelAtomLocalSql = R"sql(
    SELECT
        serial_id, sampling_size, distance_and_map_value_list,
        amplitude_estimate_ols, width_estimate_ols,
        amplitude_estimate_mdpde, width_estimate_mdpde, alpha_r
    FROM model_atom_local_potential WHERE key_tag = ?;
)sql"sv;

inline constexpr auto kSelectModelBondLocalSql = R"sql(
    SELECT
        atom_serial_id_1, atom_serial_id_2, sampling_size, distance_and_map_value_list,
        amplitude_estimate_ols, width_estimate_ols,
        amplitude_estimate_mdpde, width_estimate_mdpde, alpha_r
    FROM model_bond_local_potential WHERE key_tag = ?;
)sql"sv;

inline constexpr auto kSelectModelAtomPosteriorSql = R"sql(
    SELECT
        serial_id, amplitude_estimate_posterior, width_estimate_posterior,
        amplitude_variance_posterior, width_variance_posterior, outlier_tag, statistical_distance
    FROM model_atom_posterior WHERE key_tag = ? AND class_key = ?;
)sql"sv;

inline constexpr auto kSelectModelBondPosteriorSql = R"sql(
    SELECT
        atom_serial_id_1, atom_serial_id_2, amplitude_estimate_posterior, width_estimate_posterior,
        amplitude_variance_posterior, width_variance_posterior, outlier_tag, statistical_distance
    FROM model_bond_posterior WHERE key_tag = ? AND class_key = ?;
)sql"sv;

inline constexpr auto kSelectModelAtomGroupSql = R"sql(
    SELECT
        group_key, member_size, amplitude_estimate_mean, width_estimate_mean,
        amplitude_estimate_mdpde, width_estimate_mdpde,
        amplitude_estimate_prior, width_estimate_prior,
        amplitude_variance_prior, width_variance_prior, alpha_g
    FROM model_atom_group_potential WHERE key_tag = ? AND class_key = ?;
)sql"sv;

inline constexpr auto kSelectModelBondGroupSql = R"sql(
    SELECT
        group_key, member_size, amplitude_estimate_mean, width_estimate_mean,
        amplitude_estimate_mdpde, width_estimate_mdpde,
        amplitude_estimate_prior, width_estimate_prior,
        amplitude_variance_prior, width_variance_prior, alpha_g
    FROM model_bond_group_potential WHERE key_tag = ? AND class_key = ?;
)sql"sv;

class SQLiteStatementBatch
{
    rhbm_gem::SQLiteWrapper & m_database;
    rhbm_gem::SQLiteWrapper::StatementGuard m_guard;

public:
    SQLiteStatementBatch(rhbm_gem::SQLiteWrapper & database, const std::string & sql) :
        m_database{ database },
        m_guard{ database }
    {
        m_database.Prepare(sql);
    }

    template <typename Binder>
    void Execute(Binder && binder)
    {
        std::forward<Binder>(binder)(m_database);
        m_database.StepOnce();
        m_database.Reset();
    }
};

void DeleteRowsForKey(
    rhbm_gem::SQLiteWrapper & database,
    const std::string & table_name,
    const std::string & key_tag)
{
    database.Prepare(
        std::string(kDeleteRowsForKeySqlPrefix) + table_name + std::string(kDeleteRowsForKeySqlSuffix));
    rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    database.StepOnce();
}

void SaveModelObjectRow(
    rhbm_gem::SQLiteWrapper & database,
    const rhbm_gem::ModelObject & model_obj,
    const std::string & key_tag)
{
    SQLiteStatementBatch batch{ database, std::string(kUpsertModelObjectSql) };
    batch.Execute([&](rhbm_gem::SQLiteWrapper & statement_db)
    {
        statement_db.Bind<std::string>(1, key_tag);
        statement_db.Bind<int>(2, static_cast<int>(model_obj.GetNumberOfAtom()));
        statement_db.Bind<std::string>(3, model_obj.GetPdbID());
        statement_db.Bind<std::string>(4, model_obj.GetEmdID());
        statement_db.Bind<double>(5, model_obj.GetResolution());
        statement_db.Bind<std::string>(6, model_obj.GetResolutionMethod());
    });
}

void SaveChainMap(
    rhbm_gem::SQLiteWrapper & database,
    const rhbm_gem::ModelObject & model_obj,
    const std::string & key_tag)
{
    SQLiteStatementBatch batch{ database, std::string(kInsertModelChainMapSql) };
    for (const auto & chain_entry : model_obj.GetChainIDListMap())
    {
        const auto & entity_id{ chain_entry.first };
        const auto & chain_id_list{ chain_entry.second };
        for (size_t ordinal = 0; ordinal < chain_id_list.size(); ordinal++)
        {
            batch.Execute([&](rhbm_gem::SQLiteWrapper & statement_db)
            {
                statement_db.Bind<std::string>(1, key_tag);
                statement_db.Bind<std::string>(2, entity_id);
                statement_db.Bind<int>(3, static_cast<int>(ordinal));
                statement_db.Bind<std::string>(4, chain_id_list.at(ordinal));
            });
        }
    }
}

void SaveChemicalComponentEntryList(
    rhbm_gem::SQLiteWrapper & database,
    const rhbm_gem::ModelObject & model_obj,
    const std::string & key_tag)
{
    SQLiteStatementBatch batch{ database, std::string(kInsertModelComponentSql) };
    for (const auto & component_map_entry : model_obj.GetChemicalComponentEntryMap())
    {
        const auto & component_key{ component_map_entry.first };
        const auto & component_entry{ component_map_entry.second };
        batch.Execute([&](rhbm_gem::SQLiteWrapper & statement_db)
        {
            statement_db.Bind<std::string>(1, key_tag);
            statement_db.Bind<ComponentKey>(2, component_key);
            statement_db.Bind<std::string>(3, component_entry->GetComponentId());
            statement_db.Bind<std::string>(4, component_entry->GetComponentName());
            statement_db.Bind<std::string>(5, component_entry->GetComponentType());
            statement_db.Bind<std::string>(6, component_entry->GetComponentFormula());
            statement_db.Bind<double>(
                7, static_cast<double>(component_entry->GetComponentMolecularWeight()));
            statement_db.Bind<int>(8, static_cast<int>(component_entry->IsStandardMonomer()));
        });
    }
}

void SaveComponentAtomEntryList(
    rhbm_gem::SQLiteWrapper & database,
    const rhbm_gem::ModelObject & model_obj,
    const std::string & key_tag)
{
    SQLiteStatementBatch batch{ database, std::string(kInsertModelComponentAtomSql) };
    for (const auto & component_map_entry : model_obj.GetChemicalComponentEntryMap())
    {
        const auto & component_key{ component_map_entry.first };
        const auto & component_entry{ component_map_entry.second };
        for (const auto & atom_map_entry : component_entry->GetComponentAtomEntryMap())
        {
            const auto & atom_key{ atom_map_entry.first };
            const auto & atom_entry{ atom_map_entry.second };
            batch.Execute([&](rhbm_gem::SQLiteWrapper & statement_db)
            {
                statement_db.Bind<std::string>(1, key_tag);
                statement_db.Bind<ComponentKey>(2, component_key);
                statement_db.Bind<AtomKey>(3, atom_key);
                statement_db.Bind<std::string>(4, atom_entry.atom_id);
                statement_db.Bind<int>(5, static_cast<int>(atom_entry.element_type));
                statement_db.Bind<int>(6, static_cast<int>(atom_entry.aromatic_atom_flag));
                statement_db.Bind<int>(7, static_cast<int>(atom_entry.stereo_config));
            });
        }
    }
}

void SaveComponentBondEntryList(
    rhbm_gem::SQLiteWrapper & database,
    const rhbm_gem::ModelObject & model_obj,
    const std::string & key_tag)
{
    SQLiteStatementBatch batch{ database, std::string(kInsertModelComponentBondSql) };
    for (const auto & component_map_entry : model_obj.GetChemicalComponentEntryMap())
    {
        const auto & component_key{ component_map_entry.first };
        const auto & component_entry{ component_map_entry.second };
        for (const auto & bond_map_entry : component_entry->GetComponentBondEntryMap())
        {
            const auto & bond_key{ bond_map_entry.first };
            const auto & bond_entry{ bond_map_entry.second };
            batch.Execute([&](rhbm_gem::SQLiteWrapper & statement_db)
            {
                statement_db.Bind<std::string>(1, key_tag);
                statement_db.Bind<ComponentKey>(2, component_key);
                statement_db.Bind<BondKey>(3, bond_key);
                statement_db.Bind<std::string>(4, bond_entry.bond_id);
                statement_db.Bind<int>(5, static_cast<int>(bond_entry.bond_type));
                statement_db.Bind<int>(6, static_cast<int>(bond_entry.bond_order));
                statement_db.Bind<int>(7, static_cast<int>(bond_entry.aromatic_atom_flag));
                statement_db.Bind<int>(8, static_cast<int>(bond_entry.stereo_config));
            });
        }
    }
}

void SaveAtomObjectList(
    rhbm_gem::SQLiteWrapper & database,
    const rhbm_gem::ModelObject & model_obj,
    const std::string & key_tag)
{
    SQLiteStatementBatch batch{ database, std::string(kInsertModelAtomSql) };
    for (const auto & atom_object : model_obj.GetAtomList())
    {
        batch.Execute([&](rhbm_gem::SQLiteWrapper & statement_db)
        {
            statement_db.Bind<std::string>(1, key_tag);
            statement_db.Bind<int>(2, atom_object->GetSerialID());
            statement_db.Bind<int>(3, atom_object->GetSequenceID());
            statement_db.Bind<std::string>(4, atom_object->GetComponentID());
            statement_db.Bind<std::string>(5, atom_object->GetAtomID());
            statement_db.Bind<std::string>(6, atom_object->GetChainID());
            statement_db.Bind<std::string>(7, atom_object->GetIndicator());
            statement_db.Bind<double>(8, static_cast<double>(atom_object->GetOccupancy()));
            statement_db.Bind<double>(9, static_cast<double>(atom_object->GetTemperature()));
            statement_db.Bind<int>(10, static_cast<int>(atom_object->GetElement()));
            statement_db.Bind<int>(11, static_cast<int>(atom_object->GetStructure()));
            statement_db.Bind<int>(12, static_cast<int>(atom_object->GetSpecialAtomFlag()));
            statement_db.Bind<double>(13, static_cast<double>(atom_object->GetPosition().at(0)));
            statement_db.Bind<double>(14, static_cast<double>(atom_object->GetPosition().at(1)));
            statement_db.Bind<double>(15, static_cast<double>(atom_object->GetPosition().at(2)));
            statement_db.Bind<int>(16, static_cast<int>(atom_object->GetComponentKey()));
            statement_db.Bind<int>(17, static_cast<int>(atom_object->GetAtomKey()));
        });
    }
}

void SaveBondObjectList(
    rhbm_gem::SQLiteWrapper & database,
    const rhbm_gem::ModelObject & model_obj,
    const std::string & key_tag)
{
    SQLiteStatementBatch batch{ database, std::string(kInsertModelBondSql) };
    for (const auto & bond_object : model_obj.GetBondList())
    {
        batch.Execute([&](rhbm_gem::SQLiteWrapper & statement_db)
        {
            statement_db.Bind<std::string>(1, key_tag);
            statement_db.Bind<int>(2, bond_object->GetAtomSerialID1());
            statement_db.Bind<int>(3, bond_object->GetAtomSerialID2());
            statement_db.Bind<int>(4, bond_object->GetBondKey());
            statement_db.Bind<int>(5, static_cast<int>(bond_object->GetBondType()));
            statement_db.Bind<int>(6, static_cast<int>(bond_object->GetBondOrder()));
            statement_db.Bind<int>(7, static_cast<int>(bond_object->GetSpecialBondFlag()));
        });
    }
}

void LoadModelObjectRow(
    rhbm_gem::SQLiteWrapper & database,
    rhbm_gem::ModelObject & model_obj,
    const std::string & key_tag)
{
    database.Prepare(std::string(kSelectModelObjectSql));
    rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    const auto rc{ database.StepNext() };
    if (rc == rhbm_gem::SQLiteWrapper::StepDone())
    {
        throw std::runtime_error("Cannot find the row with key_tag = " + key_tag);
    }
    if (rc != rhbm_gem::SQLiteWrapper::StepRow())
    {
        throw std::runtime_error("Step failed: " + database.ErrorMessage());
    }

    const auto atom_size{ database.GetColumn<int>(1) };
    model_obj.SetKeyTag(database.GetColumn<std::string>(0));
    model_obj.SetPdbID(database.GetColumn<std::string>(2));
    model_obj.SetEmdID(database.GetColumn<std::string>(3));
    model_obj.SetResolution(database.GetColumn<double>(4));
    model_obj.SetResolutionMethod(database.GetColumn<std::string>(5));
    if (atom_size != static_cast<int>(model_obj.GetNumberOfAtom()))
    {
        throw std::runtime_error(
            "The number of atoms in the model object does not match the database record.");
    }
}

void LoadChainMap(
    rhbm_gem::SQLiteWrapper & database,
    rhbm_gem::ModelObject & model_obj,
    const std::string & key_tag)
{
    std::unordered_map<std::string, std::vector<std::string>> chain_map;
    database.Prepare(std::string(kSelectModelChainMapSql));
    rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ database.StepNext() };
        if (rc == rhbm_gem::SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != rhbm_gem::SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + database.ErrorMessage());
        }
        auto & chain_id_list{ chain_map[database.GetColumn<std::string>(0)] };
        const auto ordinal{ database.GetColumn<int>(1) };
        if (static_cast<int>(chain_id_list.size()) <= ordinal)
        {
            chain_id_list.resize(static_cast<size_t>(ordinal + 1));
        }
        chain_id_list.at(static_cast<size_t>(ordinal)) = database.GetColumn<std::string>(2);
    }
    model_obj.SetChainIDListMap(chain_map);
}

void LoadChemicalComponentEntryList(
    rhbm_gem::SQLiteWrapper & database,
    rhbm_gem::ModelObject & model_obj,
    const std::string & key_tag)
{
    database.Prepare(std::string(kSelectModelComponentSql));
    rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ database.StepNext() };
        if (rc == rhbm_gem::SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != rhbm_gem::SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + database.ErrorMessage());
        }

        auto component_entry{ std::make_unique<rhbm_gem::ChemicalComponentEntry>() };
        const auto component_key{ database.GetColumn<ComponentKey>(0) };
        const auto component_id{ database.GetColumn<std::string>(1) };
        component_entry->SetComponentId(component_id);
        component_entry->SetComponentName(database.GetColumn<std::string>(2));
        component_entry->SetComponentType(database.GetColumn<std::string>(3));
        component_entry->SetComponentFormula(database.GetColumn<std::string>(4));
        component_entry->SetComponentMolecularWeight(
            static_cast<float>(database.GetColumn<double>(5)));
        component_entry->SetStandardMonomerFlag(static_cast<bool>(database.GetColumn<int>(6)));
        model_obj.AddChemicalComponentEntry(component_key, std::move(component_entry));
        model_obj.GetComponentKeySystemPtr()->RegisterComponent(component_id, component_key);
    }
}

void LoadComponentAtomEntryList(
    rhbm_gem::SQLiteWrapper & database,
    rhbm_gem::ModelObject & model_obj,
    const std::string & key_tag)
{
    database.Prepare(std::string(kSelectModelComponentAtomSql));
    rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ database.StepNext() };
        if (rc == rhbm_gem::SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != rhbm_gem::SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + database.ErrorMessage());
        }

        const auto component_key{ database.GetColumn<ComponentKey>(0) };
        if (model_obj.GetChemicalComponentEntryMap().find(component_key)
            == model_obj.GetChemicalComponentEntryMap().end())
        {
            continue;
        }
        const rhbm_gem::ComponentAtomEntry atom_entry{
            database.GetColumn<std::string>(2),
            static_cast<Element>(database.GetColumn<int>(3)),
            static_cast<bool>(database.GetColumn<int>(4)),
            static_cast<StereoChemistry>(database.GetColumn<int>(5))
        };
        const auto atom_key{ database.GetColumn<AtomKey>(1) };
        model_obj.GetChemicalComponentEntryMap().at(component_key)->AddComponentAtomEntry(
            atom_key, atom_entry);
        model_obj.GetAtomKeySystemPtr()->RegisterAtom(atom_entry.atom_id, atom_key);
    }
}

void LoadComponentBondEntryList(
    rhbm_gem::SQLiteWrapper & database,
    rhbm_gem::ModelObject & model_obj,
    const std::string & key_tag)
{
    database.Prepare(std::string(kSelectModelComponentBondSql));
    rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ database.StepNext() };
        if (rc == rhbm_gem::SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != rhbm_gem::SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + database.ErrorMessage());
        }

        const auto component_key{ database.GetColumn<ComponentKey>(0) };
        if (model_obj.GetChemicalComponentEntryMap().find(component_key)
            == model_obj.GetChemicalComponentEntryMap().end())
        {
            continue;
        }
        rhbm_gem::ComponentBondEntry bond_entry;
        const auto bond_key{ database.GetColumn<BondKey>(1) };
        bond_entry.bond_id = database.GetColumn<std::string>(2);
        bond_entry.bond_type = static_cast<BondType>(database.GetColumn<int>(3));
        bond_entry.bond_order = static_cast<BondOrder>(database.GetColumn<int>(4));
        bond_entry.aromatic_atom_flag = static_cast<bool>(database.GetColumn<int>(5));
        bond_entry.stereo_config =
            static_cast<StereoChemistry>(database.GetColumn<int>(6));
        model_obj.GetChemicalComponentEntryMap().at(component_key)->AddComponentBondEntry(
            bond_key, bond_entry);
        model_obj.GetBondKeySystemPtr()->RegisterBond(bond_entry.bond_id, bond_key);
    }
}

std::vector<std::unique_ptr<rhbm_gem::AtomObject>> LoadAtomObjectList(
    rhbm_gem::SQLiteWrapper & database,
    const std::string & key_tag)
{
    std::vector<std::unique_ptr<rhbm_gem::AtomObject>> atom_object_list;

    database.Prepare(std::string(kSelectModelAtomSql));
    rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ database.StepNext() };
        if (rc == rhbm_gem::SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != rhbm_gem::SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + database.ErrorMessage());
        }

        auto atom_object{ std::make_unique<rhbm_gem::AtomObject>() };
        atom_object->SetSerialID(database.GetColumn<int>(0));
        atom_object->SetSequenceID(database.GetColumn<int>(1));
        atom_object->SetComponentID(database.GetColumn<std::string>(2));
        atom_object->SetAtomID(database.GetColumn<std::string>(3));
        atom_object->SetChainID(database.GetColumn<std::string>(4));
        atom_object->SetIndicator(database.GetColumn<std::string>(5));
        atom_object->SetOccupancy(static_cast<float>(database.GetColumn<double>(6)));
        atom_object->SetTemperature(static_cast<float>(database.GetColumn<double>(7)));
        atom_object->SetElement(static_cast<Element>(database.GetColumn<int>(8)));
        atom_object->SetStructure(static_cast<Structure>(database.GetColumn<int>(9)));
        atom_object->SetSpecialAtomFlag(static_cast<bool>(database.GetColumn<int>(10)));
        atom_object->SetPosition(
            static_cast<float>(database.GetColumn<double>(11)),
            static_cast<float>(database.GetColumn<double>(12)),
            static_cast<float>(database.GetColumn<double>(13)));
        atom_object->SetComponentKey(database.GetColumn<ComponentKey>(14));
        atom_object->SetAtomKey(database.GetColumn<AtomKey>(15));
        atom_object->SetSelectedFlag(false);
        atom_object_list.emplace_back(std::move(atom_object));
    }
    return atom_object_list;
}

std::vector<std::unique_ptr<rhbm_gem::BondObject>> LoadBondObjectList(
    rhbm_gem::SQLiteWrapper & database,
    const std::string & key_tag,
    const rhbm_gem::ModelObject & model_obj)
{
    std::vector<std::unique_ptr<rhbm_gem::BondObject>> bond_object_list;

    database.Prepare(std::string(kSelectModelBondSql));
    rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ database.StepNext() };
        if (rc == rhbm_gem::SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != rhbm_gem::SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + database.ErrorMessage());
        }

        auto atom_object_1{ model_obj.GetAtomPtr(database.GetColumn<int>(0)) };
        auto atom_object_2{ model_obj.GetAtomPtr(database.GetColumn<int>(1)) };
        auto bond_object{ std::make_unique<rhbm_gem::BondObject>(atom_object_1, atom_object_2) };
        bond_object->SetBondKey(database.GetColumn<BondKey>(2));
        bond_object->SetBondType(static_cast<BondType>(database.GetColumn<int>(3)));
        bond_object->SetBondOrder(static_cast<BondOrder>(database.GetColumn<int>(4)));
        bond_object->SetSpecialBondFlag(static_cast<bool>(database.GetColumn<int>(5)));
        bond_object->SetSelectedFlag(false);
        bond_object_list.emplace_back(std::move(bond_object));
    }
    return bond_object_list;
}

void SaveStructure(
    rhbm_gem::SQLiteWrapper & database,
    const rhbm_gem::ModelObject & model_obj,
    const std::string & key_tag)
{
    SaveModelObjectRow(database, model_obj, key_tag);
    SaveChainMap(database, model_obj, key_tag);
    SaveChemicalComponentEntryList(database, model_obj, key_tag);
    SaveComponentAtomEntryList(database, model_obj, key_tag);
    SaveComponentBondEntryList(database, model_obj, key_tag);
    SaveAtomObjectList(database, model_obj, key_tag);
    SaveBondObjectList(database, model_obj, key_tag);
}

void LoadStructure(
    rhbm_gem::SQLiteWrapper & database,
    rhbm_gem::ModelObject & model_obj,
    const std::string & key_tag)
{
    LoadChemicalComponentEntryList(database, model_obj, key_tag);
    LoadComponentAtomEntryList(database, model_obj, key_tag);
    LoadComponentBondEntryList(database, model_obj, key_tag);
    model_obj.SetAtomList(LoadAtomObjectList(database, key_tag));
    model_obj.SetBondList(LoadBondObjectList(database, key_tag, model_obj));
    LoadModelObjectRow(database, model_obj, key_tag);
    LoadChainMap(database, model_obj, key_tag);
}

void SaveAtomLocalPotentialEntryList(
    rhbm_gem::SQLiteWrapper & database,
    const rhbm_gem::ModelObject & model_obj,
    const std::string & key_tag)
{
    SQLiteStatementBatch batch{ database, std::string(kInsertModelAtomLocalSql) };
    for (const auto & atom_object : model_obj.GetAtomList())
    {
        auto entry{ atom_object->GetLocalPotentialEntry() };
        if (entry == nullptr) continue;

        batch.Execute([&](rhbm_gem::SQLiteWrapper & statement_db)
        {
            statement_db.Bind<std::string>(1, key_tag);
            statement_db.Bind<int>(2, atom_object->GetSerialID());
            statement_db.Bind<int>(3, entry->GetDistanceAndMapValueListSize());
            statement_db.Bind<std::vector<std::tuple<float, float>>>(
                4, entry->GetDistanceAndMapValueList());
            statement_db.Bind<double>(5, entry->GetAmplitudeEstimateOLS());
            statement_db.Bind<double>(6, entry->GetWidthEstimateOLS());
            statement_db.Bind<double>(7, entry->GetAmplitudeEstimateMDPDE());
            statement_db.Bind<double>(8, entry->GetWidthEstimateMDPDE());
            statement_db.Bind<double>(9, entry->GetAlphaR());
        });
    }
}

void SaveBondLocalPotentialEntryList(
    rhbm_gem::SQLiteWrapper & database,
    const rhbm_gem::ModelObject & model_obj,
    const std::string & key_tag)
{
    SQLiteStatementBatch batch{ database, std::string(kInsertModelBondLocalSql) };
    for (const auto & bond_object : model_obj.GetBondList())
    {
        auto entry{ bond_object->GetLocalPotentialEntry() };
        if (entry == nullptr) continue;

        batch.Execute([&](rhbm_gem::SQLiteWrapper & statement_db)
        {
            statement_db.Bind<std::string>(1, key_tag);
            statement_db.Bind<int>(2, bond_object->GetAtomSerialID1());
            statement_db.Bind<int>(3, bond_object->GetAtomSerialID2());
            statement_db.Bind<int>(4, entry->GetDistanceAndMapValueListSize());
            statement_db.Bind<std::vector<std::tuple<float, float>>>(
                5, entry->GetDistanceAndMapValueList());
            statement_db.Bind<double>(6, entry->GetAmplitudeEstimateOLS());
            statement_db.Bind<double>(7, entry->GetWidthEstimateOLS());
            statement_db.Bind<double>(8, entry->GetAmplitudeEstimateMDPDE());
            statement_db.Bind<double>(9, entry->GetWidthEstimateMDPDE());
            statement_db.Bind<double>(10, entry->GetAlphaR());
        });
    }
}

void SaveAtomLocalPotentialEntrySubList(
    rhbm_gem::SQLiteWrapper & database,
    const rhbm_gem::ModelObject & model_obj,
    const std::string & key_tag,
    const std::string & class_key)
{
    SQLiteStatementBatch batch{ database, std::string(kInsertModelAtomPosteriorSql) };
    for (const auto & atom_object : model_obj.GetAtomList())
    {
        auto entry{ atom_object->GetLocalPotentialEntry() };
        if (entry == nullptr) continue;

        batch.Execute([&](rhbm_gem::SQLiteWrapper & statement_db)
        {
            statement_db.Bind<std::string>(1, key_tag);
            statement_db.Bind<std::string>(2, class_key);
            statement_db.Bind<int>(3, atom_object->GetSerialID());
            statement_db.Bind<double>(4, entry->GetAmplitudeEstimatePosterior(class_key));
            statement_db.Bind<double>(5, entry->GetWidthEstimatePosterior(class_key));
            statement_db.Bind<double>(6, entry->GetAmplitudeVariancePosterior(class_key));
            statement_db.Bind<double>(7, entry->GetWidthVariancePosterior(class_key));
            statement_db.Bind<int>(8, static_cast<int>(entry->GetOutlierTag(class_key)));
            statement_db.Bind<double>(9, entry->GetStatisticalDistance(class_key));
        });
    }
}

void SaveBondLocalPotentialEntrySubList(
    rhbm_gem::SQLiteWrapper & database,
    const rhbm_gem::ModelObject & model_obj,
    const std::string & key_tag,
    const std::string & class_key)
{
    SQLiteStatementBatch batch{ database, std::string(kInsertModelBondPosteriorSql) };
    for (const auto & bond_object : model_obj.GetBondList())
    {
        auto entry{ bond_object->GetLocalPotentialEntry() };
        if (entry == nullptr) continue;

        batch.Execute([&](rhbm_gem::SQLiteWrapper & statement_db)
        {
            statement_db.Bind<std::string>(1, key_tag);
            statement_db.Bind<std::string>(2, class_key);
            statement_db.Bind<int>(3, bond_object->GetAtomSerialID1());
            statement_db.Bind<int>(4, bond_object->GetAtomSerialID2());
            statement_db.Bind<double>(5, entry->GetAmplitudeEstimatePosterior(class_key));
            statement_db.Bind<double>(6, entry->GetWidthEstimatePosterior(class_key));
            statement_db.Bind<double>(7, entry->GetAmplitudeVariancePosterior(class_key));
            statement_db.Bind<double>(8, entry->GetWidthVariancePosterior(class_key));
            statement_db.Bind<int>(9, static_cast<int>(entry->GetOutlierTag(class_key)));
            statement_db.Bind<double>(10, entry->GetStatisticalDistance(class_key));
        });
    }
}

void SaveAtomGroupPotentialEntryList(
    rhbm_gem::SQLiteWrapper & database,
    const rhbm_gem::GroupPotentialEntry & group_entry,
    const std::string & key_tag,
    const std::string & class_key)
{
    SQLiteStatementBatch batch{ database, std::string(kInsertModelAtomGroupSql) };
    for (const auto & group_key : group_entry.GetGroupKeySet())
    {
        batch.Execute([&](rhbm_gem::SQLiteWrapper & statement_db)
        {
            statement_db.Bind<std::string>(1, key_tag);
            statement_db.Bind<std::string>(2, class_key);
            statement_db.Bind<GroupKey>(3, group_key);
            statement_db.Bind<int>(4, group_entry.GetAtomObjectPtrListSize(group_key));
            statement_db.Bind<double>(5, std::get<0>(group_entry.GetGausEstimateMean(group_key)));
            statement_db.Bind<double>(6, std::get<1>(group_entry.GetGausEstimateMean(group_key)));
            statement_db.Bind<double>(7, std::get<0>(group_entry.GetGausEstimateMDPDE(group_key)));
            statement_db.Bind<double>(8, std::get<1>(group_entry.GetGausEstimateMDPDE(group_key)));
            statement_db.Bind<double>(9, std::get<0>(group_entry.GetGausEstimatePrior(group_key)));
            statement_db.Bind<double>(10, std::get<1>(group_entry.GetGausEstimatePrior(group_key)));
            statement_db.Bind<double>(11, std::get<0>(group_entry.GetGausVariancePrior(group_key)));
            statement_db.Bind<double>(12, std::get<1>(group_entry.GetGausVariancePrior(group_key)));
            statement_db.Bind<double>(13, group_entry.GetAlphaG(group_key));
        });
    }
}

void SaveBondGroupPotentialEntryList(
    rhbm_gem::SQLiteWrapper & database,
    const rhbm_gem::GroupPotentialEntry & group_entry,
    const std::string & key_tag,
    const std::string & class_key)
{
    SQLiteStatementBatch batch{ database, std::string(kInsertModelBondGroupSql) };
    for (const auto & group_key : group_entry.GetGroupKeySet())
    {
        batch.Execute([&](rhbm_gem::SQLiteWrapper & statement_db)
        {
            statement_db.Bind<std::string>(1, key_tag);
            statement_db.Bind<std::string>(2, class_key);
            statement_db.Bind<GroupKey>(3, group_key);
            statement_db.Bind<int>(4, group_entry.GetBondObjectPtrListSize(group_key));
            statement_db.Bind<double>(5, std::get<0>(group_entry.GetGausEstimateMean(group_key)));
            statement_db.Bind<double>(6, std::get<1>(group_entry.GetGausEstimateMean(group_key)));
            statement_db.Bind<double>(7, std::get<0>(group_entry.GetGausEstimateMDPDE(group_key)));
            statement_db.Bind<double>(8, std::get<1>(group_entry.GetGausEstimateMDPDE(group_key)));
            statement_db.Bind<double>(9, std::get<0>(group_entry.GetGausEstimatePrior(group_key)));
            statement_db.Bind<double>(10, std::get<1>(group_entry.GetGausEstimatePrior(group_key)));
            statement_db.Bind<double>(11, std::get<0>(group_entry.GetGausVariancePrior(group_key)));
            statement_db.Bind<double>(12, std::get<1>(group_entry.GetGausVariancePrior(group_key)));
            statement_db.Bind<double>(13, group_entry.GetAlphaG(group_key));
        });
    }
}

void LoadAtomLocalPotentialEntrySubList(
    rhbm_gem::SQLiteWrapper & database,
    const std::string & key_tag,
    const std::string & class_key,
    std::unordered_map<int, std::unique_ptr<rhbm_gem::LocalPotentialEntry>> & entry_map)
{
    database.Prepare(std::string(kSelectModelAtomPosteriorSql));
    rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    database.Bind<std::string>(2, class_key);
    while (true)
    {
        const auto rc{ database.StepNext() };
        if (rc == rhbm_gem::SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != rhbm_gem::SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + database.ErrorMessage());
        }

        const auto serial_id{ database.GetColumn<int>(0) };
        auto iter{ entry_map.find(serial_id) };
        if (iter == entry_map.end()) continue;
        auto & entry{ iter->second };
        entry->AddGausEstimatePosterior(
            class_key, database.GetColumn<double>(1), database.GetColumn<double>(2));
        entry->AddGausVariancePosterior(
            class_key, database.GetColumn<double>(3), database.GetColumn<double>(4));
        entry->AddOutlierTag(class_key, static_cast<bool>(database.GetColumn<int>(5)));
        entry->AddStatisticalDistance(class_key, database.GetColumn<double>(6));
    }
}

void LoadBondLocalPotentialEntrySubList(
    rhbm_gem::SQLiteWrapper & database,
    const std::string & key_tag,
    const std::string & class_key,
    std::map<std::pair<int, int>, std::unique_ptr<rhbm_gem::LocalPotentialEntry>> & entry_map)
{
    database.Prepare(std::string(kSelectModelBondPosteriorSql));
    rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    database.Bind<std::string>(2, class_key);
    while (true)
    {
        const auto rc{ database.StepNext() };
        if (rc == rhbm_gem::SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != rhbm_gem::SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + database.ErrorMessage());
        }

        const auto atom_pair{
            std::make_pair(database.GetColumn<int>(0), database.GetColumn<int>(1)) };
        auto iter{ entry_map.find(atom_pair) };
        if (iter == entry_map.end()) continue;
        auto & entry{ iter->second };
        entry->AddGausEstimatePosterior(
            class_key, database.GetColumn<double>(2), database.GetColumn<double>(3));
        entry->AddGausVariancePosterior(
            class_key, database.GetColumn<double>(4), database.GetColumn<double>(5));
        entry->AddOutlierTag(class_key, static_cast<bool>(database.GetColumn<int>(6)));
        entry->AddStatisticalDistance(class_key, database.GetColumn<double>(7));
    }
}

std::unordered_map<int, std::unique_ptr<rhbm_gem::LocalPotentialEntry>> LoadAtomLocalPotentialEntryMap(
    rhbm_gem::SQLiteWrapper & database,
    const std::string & key_tag)
{
    std::unordered_map<int, std::unique_ptr<rhbm_gem::LocalPotentialEntry>> entry_map;
    database.Prepare(std::string(kSelectModelAtomLocalSql));
    rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ database.StepNext() };
        if (rc == rhbm_gem::SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != rhbm_gem::SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + database.ErrorMessage());
        }

        auto entry{ std::make_unique<rhbm_gem::LocalPotentialEntry>() };
        const auto serial_id{ database.GetColumn<int>(0) };
        entry->AddDistanceAndMapValueList(
            database.GetColumn<std::vector<std::tuple<float, float>>>(2));
        entry->AddGausEstimateOLS(database.GetColumn<double>(3), database.GetColumn<double>(4));
        entry->AddGausEstimateMDPDE(database.GetColumn<double>(5), database.GetColumn<double>(6));
        entry->SetAlphaR(database.GetColumn<double>(7));
        entry_map[serial_id] = std::move(entry);
    }

    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        LoadAtomLocalPotentialEntrySubList(
            database, key_tag, ChemicalDataHelper::GetGroupAtomClassKey(i), entry_map);
    }
    return entry_map;
}

std::map<std::pair<int, int>, std::unique_ptr<rhbm_gem::LocalPotentialEntry>> LoadBondLocalPotentialEntryMap(
    rhbm_gem::SQLiteWrapper & database,
    const std::string & key_tag)
{
    std::map<std::pair<int, int>, std::unique_ptr<rhbm_gem::LocalPotentialEntry>> entry_map;
    database.Prepare(std::string(kSelectModelBondLocalSql));
    rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ database.StepNext() };
        if (rc == rhbm_gem::SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != rhbm_gem::SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + database.ErrorMessage());
        }

        auto entry{ std::make_unique<rhbm_gem::LocalPotentialEntry>() };
        const auto key{ std::make_pair(database.GetColumn<int>(0), database.GetColumn<int>(1)) };
        entry->AddDistanceAndMapValueList(
            database.GetColumn<std::vector<std::tuple<float, float>>>(3));
        entry->AddGausEstimateOLS(database.GetColumn<double>(4), database.GetColumn<double>(5));
        entry->AddGausEstimateMDPDE(database.GetColumn<double>(6), database.GetColumn<double>(7));
        entry->SetAlphaR(database.GetColumn<double>(8));
        entry_map[key] = std::move(entry);
    }

    for (size_t i = 0; i < ChemicalDataHelper::GetGroupBondClassCount(); i++)
    {
        LoadBondLocalPotentialEntrySubList(
            database, key_tag, ChemicalDataHelper::GetGroupBondClassKey(i), entry_map);
    }
    return entry_map;
}

void LoadAtomGroupPotentialEntryList(
    rhbm_gem::SQLiteWrapper & database,
    rhbm_gem::ModelObject & model_obj,
    const std::string & key_tag,
    const std::string & class_key)
{
    auto group_entry{ model_obj.GetAtomGroupPotentialEntry(class_key) };
    database.Prepare(std::string(kSelectModelAtomGroupSql));
    rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    database.Bind<std::string>(2, class_key);
    while (true)
    {
        const auto rc{ database.StepNext() };
        if (rc == rhbm_gem::SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != rhbm_gem::SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + database.ErrorMessage());
        }

        const auto group_key{ database.GetColumn<GroupKey>(0) };
        group_entry->InsertGroupKey(group_key);
        group_entry->ReserveAtomObjectPtrList(group_key, database.GetColumn<int>(1));
        group_entry->AddGausEstimateMean(group_key, database.GetColumn<double>(2), database.GetColumn<double>(3));
        group_entry->AddGausEstimateMDPDE(group_key, database.GetColumn<double>(4), database.GetColumn<double>(5));
        group_entry->AddGausEstimatePrior(group_key, database.GetColumn<double>(6), database.GetColumn<double>(7));
        group_entry->AddGausVariancePrior(group_key, database.GetColumn<double>(8), database.GetColumn<double>(9));
        group_entry->AddAlphaG(group_key, database.GetColumn<double>(10));
    }

    for (auto & atom : model_obj.GetSelectedAtomList())
    {
        group_entry->AddAtomObjectPtr(
            rhbm_gem::AtomClassifier::GetGroupKeyInClass(atom, class_key), atom);
    }
}

void LoadBondGroupPotentialEntryList(
    rhbm_gem::SQLiteWrapper & database,
    rhbm_gem::ModelObject & model_obj,
    const std::string & key_tag,
    const std::string & class_key)
{
    auto group_entry{ model_obj.GetBondGroupPotentialEntry(class_key) };
    database.Prepare(std::string(kSelectModelBondGroupSql));
    rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    database.Bind<std::string>(2, class_key);
    while (true)
    {
        const auto rc{ database.StepNext() };
        if (rc == rhbm_gem::SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != rhbm_gem::SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + database.ErrorMessage());
        }

        const auto group_key{ database.GetColumn<GroupKey>(0) };
        group_entry->InsertGroupKey(group_key);
        group_entry->ReserveBondObjectPtrList(group_key, database.GetColumn<int>(1));
        group_entry->AddGausEstimateMean(group_key, database.GetColumn<double>(2), database.GetColumn<double>(3));
        group_entry->AddGausEstimateMDPDE(group_key, database.GetColumn<double>(4), database.GetColumn<double>(5));
        group_entry->AddGausEstimatePrior(group_key, database.GetColumn<double>(6), database.GetColumn<double>(7));
        group_entry->AddGausVariancePrior(group_key, database.GetColumn<double>(8), database.GetColumn<double>(9));
        group_entry->AddAlphaG(group_key, database.GetColumn<double>(10));
    }

    for (auto & bond : model_obj.GetSelectedBondList())
    {
        group_entry->AddBondObjectPtr(
            rhbm_gem::BondClassifier::GetGroupKeyInClass(bond, class_key), bond);
    }
}

void ApplyAtomLocalPotentialEntries(
    rhbm_gem::ModelObject & model_obj,
    std::unordered_map<int, std::unique_ptr<rhbm_gem::LocalPotentialEntry>> entry_map)
{
    for (const auto & atom_object : model_obj.GetAtomList())
    {
        const auto serial_id{ atom_object->GetSerialID() };
        auto iter{ entry_map.find(serial_id) };
        if (iter != entry_map.end())
        {
            atom_object->SetSelectedFlag(true);
            atom_object->AddLocalPotentialEntry(std::move(iter->second));
        }
        else
        {
            atom_object->SetSelectedFlag(false);
        }
    }
}

void ApplyBondLocalPotentialEntries(
    rhbm_gem::ModelObject & model_obj,
    std::map<std::pair<int, int>, std::unique_ptr<rhbm_gem::LocalPotentialEntry>> entry_map)
{
    for (const auto & bond_object : model_obj.GetBondList())
    {
        const auto serial_id_pair{
            std::make_pair(bond_object->GetAtomSerialID1(), bond_object->GetAtomSerialID2()) };
        auto iter{ entry_map.find(serial_id_pair) };
        if (iter != entry_map.end())
        {
            bond_object->SetSelectedFlag(true);
            bond_object->AddLocalPotentialEntry(std::move(iter->second));
        }
        else
        {
            bond_object->SetSelectedFlag(false);
        }
    }
}

void SaveAnalysis(
    rhbm_gem::SQLiteWrapper & database,
    const rhbm_gem::ModelObject & model_obj,
    const std::string & key_tag)
{
    SaveAtomLocalPotentialEntryList(database, model_obj, key_tag);
    SaveBondLocalPotentialEntryList(database, model_obj, key_tag);

    for (const auto & [class_key, group_entry] : model_obj.GetAtomGroupPotentialEntryMap())
    {
        if (group_entry == nullptr || group_entry->GetGroupKeySet().empty())
        {
            continue;
        }
        SaveAtomLocalPotentialEntrySubList(database, model_obj, key_tag, class_key);
        SaveAtomGroupPotentialEntryList(database, *group_entry, key_tag, class_key);
    }
    for (const auto & [class_key, group_entry] : model_obj.GetBondGroupPotentialEntryMap())
    {
        if (group_entry == nullptr || group_entry->GetGroupKeySet().empty())
        {
            continue;
        }
        SaveBondLocalPotentialEntrySubList(database, model_obj, key_tag, class_key);
        SaveBondGroupPotentialEntryList(database, *group_entry, key_tag, class_key);
    }
}

void LoadAnalysis(
    rhbm_gem::SQLiteWrapper & database,
    rhbm_gem::ModelObject & model_obj,
    const std::string & key_tag)
{
    ApplyAtomLocalPotentialEntries(model_obj, LoadAtomLocalPotentialEntryMap(database, key_tag));
    ApplyBondLocalPotentialEntries(model_obj, LoadBondLocalPotentialEntryMap(database, key_tag));
    model_obj.Update();

    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        auto class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        auto group_entry{ std::make_unique<rhbm_gem::GroupPotentialEntry>() };
        model_obj.AddAtomGroupPotentialEntry(class_key, group_entry);
        LoadAtomGroupPotentialEntryList(database, model_obj, key_tag, class_key);
    }

    for (size_t i = 0; i < ChemicalDataHelper::GetGroupBondClassCount(); i++)
    {
        auto class_key{ ChemicalDataHelper::GetGroupBondClassKey(i) };
        auto group_entry{ std::make_unique<rhbm_gem::GroupPotentialEntry>() };
        model_obj.AddBondGroupPotentialEntry(class_key, group_entry);
        LoadBondGroupPotentialEntryList(database, model_obj, key_tag, class_key);
    }
}

} // namespace

namespace rhbm_gem {

ModelObjectStorage::ModelObjectStorage(SQLiteWrapper * db_manager) :
    m_database{ db_manager }
{
}

ModelObjectStorage::~ModelObjectStorage() = default;

void ModelObjectStorage::CreateTables(SQLiteWrapper & database)
{
    for (const auto create_sql : kCreateModelTableSqlList)
    {
        database.Execute(std::string(create_sql));
    }
}

void ModelObjectStorage::Save(const ModelObject & model_obj, const std::string & key_tag)
{
    for (const auto table_name : kModelTablesScopedByKey)
    {
        DeleteRowsForKey(*m_database, std::string(table_name), key_tag);
    }

    SaveStructure(*m_database, model_obj, key_tag);
    SaveAnalysis(*m_database, model_obj, key_tag);
}

std::unique_ptr<ModelObject> ModelObjectStorage::Load(const std::string & key_tag)
{
    auto model_object{ std::make_unique<ModelObject>() };
    LoadStructure(*m_database, *model_object, key_tag);
    LoadAnalysis(*m_database, *model_object, key_tag);
    return model_object;
}

} // namespace rhbm_gem
