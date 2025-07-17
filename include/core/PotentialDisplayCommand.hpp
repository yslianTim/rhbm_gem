#pragma once

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include <CLI/CLI.hpp>

#include "CommandBase.hpp"
#include "GlobalOptions.hpp"

class DataObjectManager;
class AtomSelector;

class PotentialDisplayCommand : public CommandBase
{
public:
    struct Options
    {
        int painter_choice{ 1 };
        std::string model_key_tag_list;
        std::string ref_model_key_tag_list;
        std::string pick_chain_id;
        std::string veto_chain_id;
        std::string pick_residue;
        std::string veto_residue;
        std::string pick_element;
        std::string veto_element;
        std::string pick_remoteness;
        std::string veto_remoteness;
    };

private:
    Options m_options{};
    GlobalOptions m_globals{};
    std::vector<std::string> m_model_key_tag_list;
    std::unordered_map<std::string, std::vector<std::string>> m_ref_model_key_tag_list_map;
    std::unique_ptr<AtomSelector> m_atom_selector;

public:
    PotentialDisplayCommand(void);
    ~PotentialDisplayCommand();
    void Execute(void) override;
    void RegisterCLIOptions(CLI::App * cmd) override;
    void SetGlobalOptions(const GlobalOptions & globals) override { m_globals = globals; }

    void SetPainterChoice(int value) { m_options.painter_choice = value; }
    void SetDatabasePath(const std::string & path) { m_globals.database_path = path; }
    void SetFolderPath(const std::string & path) { m_globals.folder_path = path; }
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

};
