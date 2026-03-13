#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/core/command/ResultDumpCommand.hpp>

namespace rhbm_gem {

class DataObjectManager;
class AtomObject;
class ModelObject;
class MapObject;

namespace detail {

struct ResultDumpWorkflowContext
{
    DataObjectManager & data_manager;
    const ResultDumpCommandOptions & options;
    std::string & map_key_tag;
    std::unordered_map<std::string, std::vector<AtomObject *>> & selected_atom_list_map;
    std::vector<std::shared_ptr<ModelObject>> & model_object_list;
    std::shared_ptr<MapObject> & map_object;
};

bool ExecuteResultDumpWorkflow(ResultDumpWorkflowContext & context);

} // namespace detail

} // namespace rhbm_gem
