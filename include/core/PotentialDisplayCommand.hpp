#pragma once

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include <filesystem>
#include <CLI/CLI.hpp>

#include "CommandBase.hpp"
#include "OptionEnumClass.hpp"

class AtomSelector;

namespace rhbm_gem {

class DataObjectManager;
class ModelObject;

class PotentialDisplayCommand : public CommandBase
{
public:
    struct Options : public CommandOptions
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

private:
    Options m_options;
    std::unordered_map<std::string, std::vector<std::string>> m_ref_model_key_tag_list_map;
    std::unique_ptr<::AtomSelector> m_atom_selector;
    std::vector<std::shared_ptr<ModelObject>> m_model_object_list;
    std::unordered_map<std::string, std::vector<std::shared_ptr<ModelObject>>> m_ref_model_object_list_map;

public:
    PotentialDisplayCommand();
    ~PotentialDisplayCommand();
    bool Execute() override;
    void RegisterCLIOptionsExtend(CLI::App * cmd) override;
    void ResetRuntimeState() override;
    const CommandOptions & GetOptions() const override { return m_options; }
    CommandOptions & GetOptions() override { return m_options; }

    void SetPainterChoice(PainterType value);
    void SetModelKeyTagList(const std::string & value);
    void SetRefModelKeyTagListMap(const std::string & value);
    void SetPickChainID(const std::string & value);
    void SetVetoChainID(const std::string & value);
    void SetPickResidueType(const std::string & value);
    void SetVetoResidueType(const std::string & value);
    void SetPickElementType(const std::string & value);
    void SetVetoElementType(const std::string & value);

private:
    bool BuildDataObject();
    void RunDataObjectSelection();
    void RunDisplay();

};

} // namespace rhbm_gem
