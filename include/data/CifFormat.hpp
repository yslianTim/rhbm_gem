#pragma once

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>
#include <memory>
#include <functional>
#include <unordered_map>
#include <fstream>
#include <array>

#include "GlobalEnumClass.hpp"
#include "ModelFileFormatBase.hpp"

class AtomicModelDataBlock;
class ModelObject;
class AtomObject;

class CifFormat : public ModelFileFormatBase
{
    static constexpr float m_bond_searching_radius{ 2.0f };
    std::unique_ptr<AtomicModelDataBlock> m_data_block;
    bool m_find_chemical_component_entry{ false };
    bool m_find_component_atom_entry{ false };
    bool m_find_component_bond_entry{ false };

public:
    CifFormat(void);
    ~CifFormat();
    void LoadHeader(const std::string & filename) override;
    void PrintHeader(void) const override;
    void LoadDataArray(const std::string & filename) override;
    void SaveHeader(const ModelObject * model_object, std::ostream & stream) override;
    void SaveDataArray(const ModelObject * model_object, std::ostream & stream, int par) override;
    AtomicModelDataBlock * GetDataBlockPtr(void) override;

private:
    void LoadChemicalComponentBlock(std::ifstream & infile);
    void LoadChemicalComponentAtomBlock(std::ifstream & infile);
    void LoadChemicalComponentBondBlock(std::ifstream & infile);
    void LoadDatabaseBlock(std::ifstream & infile);
    void LoadEntityBlock(std::ifstream & infile);
    void LoadPdbxData(std::ifstream & infile);
    void LoadXRayResolutionInfo(std::ifstream & infile);
    void LoadAtomTypeBlock(std::ifstream & infile);
    void LoadStructureConformationBlock(std::ifstream & infile);
    void LoadStructureConnectionBlock(std::ifstream & infile);
    void LoadStructureSheetBlock(std::ifstream & infile);
    void LoadAtomSiteBlock(std::ifstream & infile);
    void ConstructBondList(void);
    void ParseLoopBlock(std::ifstream & infile,
        std::string_view data_block_prefix,
        const std::function<void(const std::unordered_map<std::string, size_t> &,
                                 const std::vector<std::string> &)> & row_handler
    );
    void WriteAtomSiteBlock(const ModelObject * model_object, std::ostream & stream, int model_par);
    void WriteAtomSiteBlockEntry(const AtomObject * atom,
        const std::array<float, 3> & position,
        const std::string & alt_id,
        float occupancy,
        float temperature,
        int model_number,
        std::ostream & stream
    );
    void BuildDefaultChemicalComponentEntry(const std::string & comp_id);
    void BuildDefaultComponentAtomEntry(
        const std::string & comp_id,
        const std::string & atom_id,
        const std::string & element_symbol
    );
    void BuildDefaultComponentBondEntry(void);
    void BuildPepetideBondEntry(void);
    void BuildPhosphodiesterBondEntry(void);

};
