#pragma once

#include <functional>
#include <string>
#include <string_view>
#include <vector>

namespace rhbm_gem {

class SQLiteWrapper;

namespace persistence {

struct ManagedStoreDescriptor
{
    std::string_view object_type;
    std::vector<std::string_view> managed_table_names;
    std::function<void(SQLiteWrapper &)> ensure_schema_v2;
    std::function<void(SQLiteWrapper &)> validate_schema_v2;
    std::function<std::vector<std::string>(SQLiteWrapper &)> list_keys;
    std::function<std::vector<std::string>(SQLiteWrapper &, const std::string &)> list_keys_with_suffix;
    std::function<void(SQLiteWrapper &, const std::string &)> create_shadow_tables_v2;
    std::function<void(SQLiteWrapper &, const std::string &)> copy_into_shadow_tables_v2;
    std::function<void(SQLiteWrapper &, const std::string &)> drop_old_and_rename_shadow_tables_v2;
};

const std::vector<ManagedStoreDescriptor> & GetAllManagedStoreDescriptors();
const ManagedStoreDescriptor & LookupManagedStoreDescriptor(std::string_view object_type);

} // namespace persistence

} // namespace rhbm_gem
