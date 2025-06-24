#pragma once

#include <string>
#include <string_view>
#include <vector>
#include <memory>
#include <functional>
#include <unordered_map>
#include "ModelFileFormatBase.hpp"

enum class Element : uint16_t;
class AtomicModelDataBlock;

class CifFormat : public ModelFileFormatBase
{
    std::unique_ptr<AtomicModelDataBlock> m_data_block;

public:
    CifFormat(void);
    ~CifFormat();
    void LoadHeader(const std::string & filename) override;
    void PrintHeader(void) const override;
    void LoadDataArray(const std::string & filename) override;
    AtomicModelDataBlock * GetDataBlockPtr(void) override;

private:
    void LoadDatabaseInfo(const std::string & filename);
    void LoadEntityInfo(const std::string & filename);
    void LoadPdbxData(const std::string & filename);
    void LoadElementTypeList(const std::string & filename);
    void LoadStructHelixInfo(const std::string & filename);
    void LoadStructSheetInfo(const std::string & filename);
    void LoadAtomSiteData(const std::string & filename);
    void ParseLoopBlock(std::ifstream & infile,
        std::string_view data_block_prefix,
        const std::function<void(const std::unordered_map<std::string, size_t> &,
                                 const std::vector<std::string> &)> & row_handler);

};
