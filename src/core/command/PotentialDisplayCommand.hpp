#pragma once

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include <filesystem>

#include "detail/CommandBase.hpp"

class AtomSelector;

namespace rhbm_gem {

class DataObjectManager;
class ModelObject;

class PotentialDisplayCommand : public CommandWithRequest<PotentialDisplayRequest>
{
private:
    std::unique_ptr<::AtomSelector> m_atom_selector;
    std::vector<std::shared_ptr<ModelObject>> m_model_object_list;
    std::unordered_map<std::string, std::vector<std::shared_ptr<ModelObject>>> m_ref_model_object_list_map;

public:
    PotentialDisplayCommand();
    ~PotentialDisplayCommand() override;

private:
    void NormalizeRequest() override;
    void ResetRuntimeState() override;
    bool ExecuteImpl() override;
    bool BuildDataObject();
    void RunDataObjectSelection();
    void RunDisplay();

};

} // namespace rhbm_gem
