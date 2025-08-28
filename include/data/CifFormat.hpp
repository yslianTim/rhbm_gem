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
#include "ModelFileFormatBase.hpp"

enum class Element : uint16_t;
class AtomicModelDataBlock;
class ModelObject;
class AtomObject;

class CifFormat : public ModelFileFormatBase
{
    std::unique_ptr<AtomicModelDataBlock> m_data_block;

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
    void LoadChemicalComponentInfo(std::ifstream & infile);
    void LoadDatabaseInfo(std::ifstream & infile);
    void LoadEntityInfo(std::ifstream & infile);
    void LoadPdbxData(std::ifstream & infile);
    void LoadXRayResolutionInfo(std::ifstream & infile);
    void LoadElementTypeList(std::ifstream & infile);
    void LoadStructHelixInfo(std::ifstream & infile);
    void LoadStructSheetInfo(std::ifstream & infile);
    void LoadAtomSiteData(std::ifstream & infile);
    void ParseLoopBlock(std::ifstream & infile,
        std::string_view data_block_prefix,
        const std::function<void(const std::unordered_map<std::string, size_t> &,
                                 const std::vector<std::string> &)> & row_handler);
    void SaveAtomSiteData(const ModelObject * model_object, std::ostream & stream, int model_par);
    void WriteAtomSiteBlock(const AtomObject * atom,
        const std::array<float, 3> & position,
        const std::string & alt_id,
        float occupancy,
        float temperature,
        int model_number,
        std::ostream & stream);

};
