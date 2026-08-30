#pragma once

#include <array>
#include <cstddef>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>

#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include "data/detail/ModelObjectParts.hpp"
#include "io/sqlite/SQLiteWrapper.hpp"

namespace data_test {

namespace rg = rhbm_gem;

inline std::shared_ptr<rg::ModelObject> LoadFixtureModel(
    const std::filesystem::path & model_path,
    const std::string & key_tag = "model")
{
    auto model{ rg::ReadModel(model_path) };
    model->SetKeyTag(key_tag);
    return std::shared_ptr<rg::ModelObject>{ std::move(model) };
}

inline std::unique_ptr<rg::ModelObject> MakeModelWithBond()
{
    rg::ModelObjectParts parts;
    parts.atom_list.reserve(2);

    auto atom_1{ std::make_unique<rg::AtomObject>() };
    atom_1->SetSerialID(1);
    atom_1->SetPosition(0.0f, 0.0f, 0.0f);

    auto atom_2{ std::make_unique<rg::AtomObject>() };
    atom_2->SetSerialID(2);
    atom_2->SetPosition(1.0f, 0.0f, 0.0f);

    auto * atom_1_ptr{ atom_1.get() };
    auto * atom_2_ptr{ atom_2.get() };
    parts.atom_list.emplace_back(std::move(atom_1));
    parts.atom_list.emplace_back(std::move(atom_2));
    parts.bond_list.emplace_back(std::make_unique<rg::BondObject>(
        atom_1_ptr, atom_2_ptr));
    return std::make_unique<rg::ModelObject>(
        rg::AssembleModelObject(std::move(parts)));
}

inline rg::MapObject MakeMapObject()
{
    std::array<int, 3> grid_size{ 2, 2, 2 };
    std::array<float, 3> grid_spacing{ 1.0f, 1.0f, 1.0f };
    std::array<float, 3> origin{ 0.0f, 0.0f, 0.0f };
    auto values{ std::make_unique<float[]>(8) };
    for (size_t i = 0; i < 8; ++i)
    {
        values[i] = static_cast<float>(i + 1);
    }
    return rg::MapObject{
        grid_size, grid_spacing, origin, std::move(values) };
}

inline std::filesystem::path CopyFixtureWithNewName(
    const std::filesystem::path & source_path,
    const std::filesystem::path & output_path)
{
    std::filesystem::copy_file(
        source_path,
        output_path,
        std::filesystem::copy_options::overwrite_existing);
    return output_path;
}

inline int GetUserVersion(const std::filesystem::path & database_path)
{
    rg::SQLiteWrapper database{ database_path };
    database.Prepare("PRAGMA user_version;");
    rg::SQLiteWrapper::StatementGuard guard(database);
    if (database.StepNext() != rg::SQLiteWrapper::StepRow())
    {
        throw std::runtime_error("Failed to read SQLite user_version.");
    }
    return database.GetColumn<int>(0);
}

inline void ExecuteSql(
    const std::filesystem::path & database_path,
    const std::string & sql)
{
    rg::SQLiteWrapper database{ database_path };
    database.Execute(sql);
}

inline void ExecuteSqlWithForeignKeysOff(
    const std::filesystem::path & database_path,
    const std::string & sql)
{
    rg::SQLiteWrapper database{ database_path };
    database.Execute("PRAGMA foreign_keys = OFF;");
    database.Execute(sql);
}

inline bool HasTable(
    const std::filesystem::path & database_path,
    const std::string & table_name)
{
    rg::SQLiteWrapper database{ database_path };
    database.Prepare(
        "SELECT name FROM sqlite_master WHERE type = 'table' AND name = ? LIMIT 1;");
    rg::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, table_name);
    return database.StepNext() == rg::SQLiteWrapper::StepRow();
}

inline int CountRows(
    const std::filesystem::path & database_path,
    const std::string & table_name,
    const std::string & key_tag = "")
{
    rg::SQLiteWrapper database{ database_path };
    const auto sql{
        key_tag.empty()
            ? "SELECT COUNT(*) FROM " + table_name + ";"
            : "SELECT COUNT(*) FROM " + table_name + " WHERE key_tag = ?;" };
    database.Prepare(sql);
    rg::SQLiteWrapper::StatementGuard guard(database);
    if (!key_tag.empty())
    {
        database.Bind<std::string>(1, key_tag);
    }
    if (database.StepNext() != rg::SQLiteWrapper::StepRow())
    {
        throw std::runtime_error("Failed to count rows in " + table_name);
    }
    return database.GetColumn<int>(0);
}

inline bool HasColumn(
    const std::filesystem::path & database_path,
    const std::string & table_name,
    const std::string & column_name)
{
    rg::SQLiteWrapper database{ database_path };
    database.Prepare("PRAGMA table_info(" + table_name + ");");
    rg::SQLiteWrapper::StatementGuard guard(database);
    while (true)
    {
        const auto rc{ database.StepNext() };
        if (rc == rg::SQLiteWrapper::StepDone())
        {
            return false;
        }
        if (rc != rg::SQLiteWrapper::StepRow())
        {
            throw std::runtime_error(
                "Failed to inspect table_info for " + table_name);
        }
        if (database.GetColumn<std::string>(1) == column_name)
        {
            return true;
        }
    }
}

} // namespace data_test
