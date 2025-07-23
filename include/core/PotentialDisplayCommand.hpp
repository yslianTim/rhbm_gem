#pragma once

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include <filesystem>
#include <CLI/CLI.hpp>

#include "CommandBase.hpp"
#include "OptionEnumClass.hpp"

class DataObjectManager;
class ModelObject;
class AtomSelector;

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
        std::string pick_remoteness{ "" };
        std::string veto_remoteness{ "" };
    };

private:
    Options m_options{};
    std::unordered_map<std::string, std::vector<std::string>> m_ref_model_key_tag_list_map;
    std::unique_ptr<AtomSelector> m_atom_selector;
    std::vector<ModelObject *> m_model_object_list;
    std::vector<ModelObject *> m_ordered_model_object_list;
    std::unordered_map<std::string, std::vector<ModelObject *>> m_ref_model_object_list_map;
    std::unordered_map<std::string, std::vector<ModelObject *>> m_ordered_ref_model_object_list_map;

public:
    PotentialDisplayCommand(void);
    ~PotentialDisplayCommand();
    bool Execute(void) override;
    bool ValidateOptions(void) const override;
    void RegisterCLIOptionsExtend(CLI::App * cmd) override;
    const CommandOptions & GetOptions(void) const override { return m_options; }

    void SetPainterChoice(PainterType value) { m_options.painter_choice = value; }
    void SetDatabasePath(const std::filesystem::path & path) { m_options.database_path = path; }
    void SetFolderPath(const std::filesystem::path & path) { m_options.folder_path = path; }
    void SetModelKeyTagList(const std::string & value);
    void SetRefModelKeyTagListMap(const std::string & value);
    void SetPickChainID(const std::string & value);
    void SetPickResidueType(const std::string & value);
    void SetPickElementType(const std::string & value);
    void SetPickRemotenessType(const std::string & value);
    void SetVetoChainID(const std::string & value);
    void SetVetoResidueType(const std::string & value);
    void SetVetoElementType(const std::string & value);
    void SetVetoRemotenessType(const std::string & value);

private:
    void LoadModelObjects(DataObjectManager * data_manager);
    void LoadRefModelObjects(DataObjectManager * data_manager);
    void BuildOrderedModelObjectList(void);
    void BuildOrderedRefModelObjectListMap(void);
    void RunDisplay(void);

};
