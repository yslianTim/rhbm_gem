#pragma once

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include "ModelFileFormatBase.hpp"

class AtomicModelDataBlock;

class PdbFormat : public ModelFileFormatBase
{
    enum class PDB_HEADER
    {
        UNK,
        // Title Section
        HEADER, OBSLTE, TITLE, SPLT, CAVEAT, COMPND, SOURCE, KEYWDS, EXPDTA,
        NUMMDL, MDLTYP, AUTHOR, REVDAT, SPRSDE, JRNL, REMARK,
        // Primary Structure Section
        DBREF, SEQADV, SEQRES, MODRES,
        // Heterogen Section
        HET, HETNAM, HETSYN, FORMUL,
        // Secondary Structure Section
        HELIX, SHEET,
        // Connectivity Annotation Section
        SSBOND, LINK, CISPEP,
        // Miscellaneous Features Section
        SITE,
        // Crystallographic and Coordinate Transformation Section
        CRYST1, ORIGX1, ORIGX2, ORIGX3, SCALE1, SCALE2, SCALE3, MTRIX1, MTRIX2, MTRIX3,
        // Coordinate Section
        MODEL, ATOM, HETATM, ANISOU, TER, ENDMDL,
        // Connectivity Section
        CONECT,
        // Bookkeeping Section
        MASTER, END
    };

    enum class ATOM : unsigned int
    {
        SERIAL_ID  =  6, ATOM_NAME   = 12, INDICATOR  = 16,  RESIDUE_NAME = 17,
        CHAIN_ID   = 21, RESIDUE_ID  = 22, CODE       = 26,
        POSITION_X = 30, POSITION_Y  = 38, POSITION_Z = 46,
        OCCUPANCY  = 54, TEMPERATURE = 60, SEGMENT_ID = 72,  ELEMENT = 76,  CHARGE = 78
    };

    inline static std::unordered_map<std::string_view, PDB_HEADER> header_map
    {
        {"HEADER", PDB_HEADER::HEADER}, {"OBSLTE", PDB_HEADER::OBSLTE}, {"TITLE ", PDB_HEADER::TITLE },
        {"SPLT  ", PDB_HEADER::SPLT  }, {"CAVEAT", PDB_HEADER::CAVEAT}, {"COMPND", PDB_HEADER::COMPND},
        {"SOURCE", PDB_HEADER::SOURCE}, {"KEYWDS", PDB_HEADER::KEYWDS}, {"EXPDTA", PDB_HEADER::EXPDTA},
        {"NUMMDL", PDB_HEADER::NUMMDL}, {"MDLTYP", PDB_HEADER::MDLTYP}, {"AUTHOR", PDB_HEADER::AUTHOR},
        {"REVDAT", PDB_HEADER::REVDAT}, {"SPRSDE", PDB_HEADER::SPRSDE}, {"JRNL  ", PDB_HEADER::JRNL  },
        {"REMARK", PDB_HEADER::REMARK}, {"DBREF ", PDB_HEADER::DBREF }, {"SEQADV", PDB_HEADER::SEQADV},
        {"SEQRES", PDB_HEADER::SEQRES}, {"MODRES", PDB_HEADER::MODRES}, {"HET   ", PDB_HEADER::HET   },
        {"HETNAM", PDB_HEADER::HETNAM}, {"HETSYN", PDB_HEADER::HETSYN}, {"FORMUL", PDB_HEADER::FORMUL},
        {"HELIX ", PDB_HEADER::HELIX }, {"SHEET ", PDB_HEADER::SHEET }, {"SSBOND", PDB_HEADER::SSBOND},
        {"LINK  ", PDB_HEADER::LINK  }, {"CISPEP", PDB_HEADER::CISPEP}, {"SITE  ", PDB_HEADER::SITE  },
        {"CRYST1", PDB_HEADER::CRYST1}, {"ORIGX1", PDB_HEADER::ORIGX1}, {"ORIGX2", PDB_HEADER::ORIGX2},
        {"ORIGX3", PDB_HEADER::ORIGX3}, {"SCALE1", PDB_HEADER::SCALE1}, {"SCALE2", PDB_HEADER::SCALE2},
        {"SCALE3", PDB_HEADER::SCALE3}, {"MTRIX1", PDB_HEADER::MTRIX1}, {"MTRIX2", PDB_HEADER::MTRIX2},
        {"MTRIX3", PDB_HEADER::MTRIX3}, {"MODEL ", PDB_HEADER::MODEL }, {"ATOM  ", PDB_HEADER::ATOM  },
        {"HETATM", PDB_HEADER::HETATM}, {"ANISOU", PDB_HEADER::ANISOU}, {"TER   ", PDB_HEADER::TER   },
        {"ENDMDL", PDB_HEADER::ENDMDL}, {"CONECT", PDB_HEADER::CONECT}, {"MASTER", PDB_HEADER::MASTER},
        {"END   ", PDB_HEADER::END   }
    };

    struct AtomSite
    {
        int   serial_id;
        char  atom_name[5];
        char  indicator;
        char  residue_name[4];
        char  chain_id;
        int   sequence_id;
        char  code;
        float position_x, position_y, position_z;
        float occupancy;
        float temperature;
        char  segment_id[5];
        char  element[3];
        char  charge[3];
    };

    const int header_string_length{ 6 };
    std::unique_ptr<AtomicModelDataBlock> m_data_block;

public:
    PdbFormat(void);
    ~PdbFormat();
    void LoadHeader(const std::string & filename) override;
    void PrintHeader(void) const override;
    void LoadDataArray(const std::string & filename) override;
    void SaveHeader(const ModelObject * model_object, std::ostream & stream) override;
    void SaveDataArray(const ModelObject * model_object, std::ostream & stream, int par) override;
    AtomicModelDataBlock * GetDataBlockPtr(void) override;

private:
    void LoadAtomSiteData(const std::string & filename);
    void ScanAtomEntry(char * line, bool is_special, int model_number);
    PDB_HEADER MapToHeaderType(const std::string & name) const;

};
