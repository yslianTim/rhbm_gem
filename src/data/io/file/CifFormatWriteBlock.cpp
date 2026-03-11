#include "internal/io/file/CifFormat.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/LocalPotentialEntry.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>

#include <iomanip>
#include <ostream>

namespace rhbm_gem {

AtomicModelDataBlock * CifFormat::GetDataBlockPtr()
{
    return m_data_block.get();
}

void CifFormat::Write(const ModelObject & model_object, std::ostream & stream, int model_par)
{
    stream << "data_" << model_object.GetKeyTag() << '\n';
    stream << "#\n";
    WriteAtomSiteBlock(model_object, stream, model_par);
}

void CifFormat::WriteAtomSiteBlock(
    const ModelObject & model_object, std::ostream & stream, int model_par)
{
    stream << "loop_\n";
    stream << "_atom_site.group_PDB\n";
    stream << "_atom_site.id\n";
    stream << "_atom_site.type_symbol\n";
    stream << "_atom_site.label_atom_id\n";
    stream << "_atom_site.label_alt_id\n";
    stream << "_atom_site.label_comp_id\n";
    stream << "_atom_site.label_asym_id\n";
    stream << "_atom_site.label_seq_id\n";
    stream << "_atom_site.Cartn_x\n";
    stream << "_atom_site.Cartn_y\n";
    stream << "_atom_site.Cartn_z\n";
    stream << "_atom_site.occupancy\n";
    stream << "_atom_site.B_iso_or_equiv\n";
    stream << "_atom_site.pdbx_PDB_model_num\n";

    const int model_number{ 1 };
    for (const auto & atom_ptr : model_object.GetAtomList())
    {
        const AtomObject * atom{ atom_ptr.get() };
        if (atom->GetLocalPotentialEntry() == nullptr) continue;
        auto model_entry{ atom->GetLocalPotentialEntry() };
        auto gaus_estimate{ model_entry->GetGausEstimateMDPDE(model_par) };
        auto position{ atom->GetPosition() };
        WriteAtomSiteBlockEntry(
            atom, position, atom->GetIndicator(), atom->GetOccupancy(),
            static_cast<float>(gaus_estimate), model_number, stream
        );
    }
    stream << "#\n";
}

void CifFormat::WriteAtomSiteBlockEntry(
    const AtomObject * atom,
    const std::array<float, 3> & position,
    const std::string & alt_id,
    float occupancy,
    float temperature,
    int model_number,
    std::ostream & stream)
{
    std::string group_PDB{ atom->GetSpecialAtomFlag() ? "HETATM" : "ATOM" };
    std::string type_symbol{ ChemicalDataHelper::GetLabel(atom->GetElement()) };
    std::string label_atom_id{ atom->GetAtomID() };
    std::string label_alt_id{ alt_id.empty() ? "." : alt_id };
    std::string label_comp_id{ atom->GetComponentID() };
    std::string label_asym_id{ atom->GetChainID() };
    std::string label_seq_id{
        (atom->GetSequenceID() == -1) ? "." : std::to_string(atom->GetSequenceID())
    };

    stream << group_PDB << ' '
           << atom->GetSerialID() << ' '
           << type_symbol << ' '
           << label_atom_id << ' '
           << label_alt_id << ' '
           << label_comp_id << ' '
           << label_asym_id << ' '
           << label_seq_id << ' '
           << std::fixed << std::setprecision(3)
           << position[0] << ' ' << position[1] << ' ' << position[2] << ' '
           << occupancy << ' ' << temperature << ' '
           << model_number << '\n';
}

} // namespace rhbm_gem
