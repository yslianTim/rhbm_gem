#pragma once

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include <filesystem>

#include "CommandBase.hpp"
#include <rhbm_gem/core/command/CommandEnumClass.hpp>

class AtomSelector;
namespace CLI
{
    class App;
}

namespace rhbm_gem {

class DataObjectManager;
class ModelObject;
struct PotentialDisplayRequest;
void BindPotentialDisplayRequestOptions(CLI::App * command, PotentialDisplayRequest & request);

struct PotentialDisplayCommandOptions
{
    PainterType painter_choice{ PainterType::MODEL };
    std::vector<std::string> model_key_tag_list{};
    std::string ref_model_key_tag_list{ "" };
    std::string pick_chain_id{ "" };
    std::string veto_chain_id{ "" };
    std::string pick_residue{ "" };
    std::string veto_residue{ "" };
    std::string pick_element{ "" };
    std::string veto_element{ "" };
};

class PotentialDisplayCommand : public CommandBase
{
private:
    PotentialDisplayCommandOptions m_options{};
    std::unordered_map<std::string, std::vector<std::string>> m_ref_model_key_tag_list_map;
    std::unique_ptr<::AtomSelector> m_atom_selector;
    std::vector<std::shared_ptr<ModelObject>> m_model_object_list;
    std::unordered_map<std::string, std::vector<std::shared_ptr<ModelObject>>> m_ref_model_object_list_map;

public:
    explicit PotentialDisplayCommand(CommonOptionProfile profile);
    ~PotentialDisplayCommand() override;
    void ApplyRequest(const PotentialDisplayRequest & request);

private:
    void SetPainterChoice(PainterType value);
    void SetModelKeyTagList(const std::string & value);
    void SetRefModelKeyTagListMap(const std::string & value);
    void ResetRuntimeState() override;
    bool ExecuteImpl() override;
    bool BuildDataObject();
    void RunDataObjectSelection();
    void RunDisplay();

};

} // namespace rhbm_gem
