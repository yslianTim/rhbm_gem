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

#include "GlobalEnumClass.hpp"
#include "ModelFileFormatBase.hpp"

namespace rhbm_gem {

class AtomicModelDataBlock;
class ModelObject;
class AtomObject;

class CifFormat : public ModelFileFormatBase
{
    struct ParsedLoopRow
    {
        std::vector<std::string> token_list;
        size_t line_number;
    };

    struct ParsedLoopCategory
    {
        std::vector<std::string> column_name_list;
        std::vector<ParsedLoopRow> row_list;
    };

    using ColumnIndexMap = std::map<std::string, size_t, std::less<>>;

    static constexpr float m_bond_searching_radius{ 2.0f };
    std::unique_ptr<AtomicModelDataBlock> m_data_block;
    bool m_find_chemical_component_entry{ false };
    bool m_find_component_atom_entry{ false };
    bool m_find_component_bond_entry{ false };
    std::string m_cached_filename;
    bool m_has_parsed_document{ false };
    std::unordered_map<std::string, std::vector<ParsedLoopCategory>> m_loop_category_map;
    std::unordered_map<std::string, std::vector<std::string>> m_data_item_map;

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
    void LoadChemicalComponentBlock(void);
    void LoadChemicalComponentAtomBlock(void);
    void LoadChemicalComponentBondBlock(void);
    void LoadDatabaseBlock(void);
    void LoadEntityBlock(void);
    void LoadPdbxData(void);
    void LoadXRayResolutionInfo(void);
    void LoadAtomTypeBlock(void);
    void LoadStructureConformationBlock(void);
    void LoadStructureConnectionBlock(void);
    void LoadStructureSheetBlock(void);
    void LoadAtomSiteBlock(void);
    void ConstructBondList(void);
    void ParseLoopBlock(
        std::string_view data_block_prefix,
        const std::function<void(const ColumnIndexMap &,
                                 const std::vector<std::string> &)> & row_handler
    );
    void EnsureParsedDocument(const std::string & filename);
    void ParseMmCifDocument(const std::string & filename);
    void ResetParsedDocument(void);
    std::optional<std::string> GetFirstDataItemValue(std::string_view key) const;
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

} // namespace rhbm_gem
