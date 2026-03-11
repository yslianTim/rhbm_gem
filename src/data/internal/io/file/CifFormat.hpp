#pragma once

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>
#include <memory>
#include <functional>
#include <unordered_map>
#include <array>
#include <map>
#include <optional>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include "ICifCategoryParser.hpp"
#include "ModelFileFormatBase.hpp"

namespace rhbm_gem {

class AtomicModelDataBlock;
class ModelObject;
class AtomObject;
class ICifCategoryParser;

class CifFormat : public ModelFileFormatBase
{
    using ColumnIndexMap = CifColumnIndexMap;

    static constexpr float m_bond_searching_radius{ 2.0f };
    std::unique_ptr<AtomicModelDataBlock> m_data_block;
    std::unique_ptr<ICifCategoryParser> m_category_parser;
    bool m_find_chemical_component_entry{ false };
    bool m_find_component_atom_entry{ false };
    bool m_find_component_bond_entry{ false };
    std::unordered_map<std::string, std::vector<MmCifParsedLoopCategory>> m_loop_category_map;
    std::unordered_map<std::string, std::vector<std::string>> m_data_item_map;

public:
    CifFormat();
    ~CifFormat();
    void Read(std::istream & stream, const std::string & source_name) override;
    void Write(const ModelObject & model_object, std::ostream & stream, int par) override;
    AtomicModelDataBlock * GetDataBlockPtr() override;

private:
    void ResetReadState();
    void LogHeaderSummary() const;
    void LoadChemicalComponentBlock();
    void LoadChemicalComponentAtomBlock();
    void LoadChemicalComponentBondBlock();
    void LoadDatabaseBlock();
    void LoadEntityBlock();
    void LoadPdbxData();
    void LoadXRayResolutionInfo();
    void LoadAtomTypeBlock();
    void LoadStructureConformationBlock();
    void LoadStructureConnectionBlock();
    void LoadStructureSheetBlock();
    void LoadAtomSiteBlock();
    void ConstructBondList();
    void ParseLoopBlock(
        std::string_view data_block_prefix,
        const CifLoopRowHandler & row_handler
    );
    void ParseMmCifDocument(std::istream & stream, const std::string & source_name);
    std::optional<std::string> GetFirstDataItemValue(std::string_view key) const;
    void WriteAtomSiteBlock(const ModelObject & model_object, std::ostream & stream, int model_par);
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
    void BuildDefaultComponentBondEntry();
    void BuildPepetideBondEntry();
    void BuildPhosphodiesterBondEntry();

};

} // namespace rhbm_gem
