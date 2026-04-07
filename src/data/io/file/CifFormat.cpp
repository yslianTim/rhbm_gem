#include "CifFormat.hpp"

#include "ModelImportState.hpp"
#include "data/detail/LocalPotentialEntry.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ChemicalComponentEntry.hpp>
#include "data/detail/ModelAnalysisData.hpp"
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/ComponentHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>
#include <rhbm_gem/utils/math/KDTreeAlgorithm.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <initializer_list>
#include <istream>
#include <map>
#include <memory>
#include <optional>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace rhbm_gem {

struct CifFormatState {
    std::unique_ptr<ModelImportState> data_block{std::make_unique<ModelImportState>()};
    bool find_chemical_component_entry{false};
    bool find_component_atom_entry{false};
    bool find_component_bond_entry{false};
    std::unordered_map<std::string, std::vector<std::string>> data_item_map;
    std::unordered_map<std::string, std::unordered_set<std::string>>
        observed_component_atom_id_set_map;

    struct ParsedLoopRow {
        std::vector<std::string> token_list;
        size_t line_number;
    };

    struct ParsedLoopCategory {
        std::vector<std::string> column_name_list;
        std::vector<ParsedLoopRow> row_list;
    };

    std::unordered_map<std::string, std::vector<ParsedLoopCategory>> loop_category_map;
};

namespace {

using CifColumnIndexMap = std::map<std::string, size_t, std::less<>>;
using CifLoopRowHandler = std::function<void(
    const CifColumnIndexMap&,
    const std::vector<std::string>&)>;

constexpr float kBondSearchingRadius{2.0f};

bool IsMmCifMissingValue(const std::string& value) {
    return value.empty() || value == "." || value == "?";
}

template <typename IndexMap>
std::optional<std::string> GetTokenOptional(
    const IndexMap& index_map,
    const std::vector<std::string>& token_list,
    std::initializer_list<std::string_view> name_candidates) {
    for (const auto& name : name_candidates) {
        auto iter{index_map.find(name)};
        if (iter == index_map.end())
            continue;
        if (iter->second >= token_list.size())
            continue;
        return token_list[iter->second];
    }
    return std::nullopt;
}

int ParseIntOrDefault(
    const std::string& value,
    int default_value,
    const std::string& field_name,
    const std::string& log_context) {
    if (IsMmCifMissingValue(value))
        return default_value;
    try {
        return std::stoi(value);
    } catch (const std::exception&) {
        Logger::Log(
            LogLevel::Warning,
            log_context + " Invalid integer in " + field_name + ": " + value + ", fallback = " + std::to_string(default_value));
        return default_value;
    }
}

std::optional<int> TryParseInt(const std::string& value) {
    if (IsMmCifMissingValue(value))
        return std::nullopt;
    try {
        return std::stoi(value);
    } catch (const std::exception&) {
        return std::nullopt;
    }
}

std::optional<float> TryParseFloat(const std::string& value) {
    if (IsMmCifMissingValue(value))
        return std::nullopt;
    try {
        return std::stof(value);
    } catch (const std::exception&) {
        return std::nullopt;
    }
}

float ParseFloatOrDefault(
    const std::string& value,
    float default_value,
    const std::string& field_name,
    const std::string& log_context) {
    auto parsed_value{TryParseFloat(value)};
    if (parsed_value.has_value())
        return *parsed_value;
    if (!IsMmCifMissingValue(value)) {
        Logger::Log(
            LogLevel::Warning,
            log_context + " Invalid float in " + field_name + ": " + value + ", fallback = " + std::to_string(default_value));
    }
    return default_value;
}

std::string_view TrimLeft(std::string_view value) {
    size_t pos{0};
    while (pos < value.size() && std::isspace(static_cast<unsigned char>(value[pos]))) {
        ++pos;
    }
    return value.substr(pos);
}

bool StartsWithToken(std::string_view line, std::string_view token) {
    return line.rfind(token, 0) == 0;
}

std::string BuildLoopCategoryPrefix(const std::string& column_name) {
    auto dot_pos{column_name.find('.')};
    if (dot_pos == std::string::npos)
        return column_name;
    return column_name.substr(0, dot_pos + 1);
}

std::vector<std::string> SplitMmCifTokens(const std::string& line) {
    std::vector<std::string> token_list;
    std::string current_token;
    enum class State {
        InSpace,
        InUnquoted,
        InSingleQuote,
        InDoubleQuote
    };
    State state{State::InSpace};

    auto flush_token =
        [&]() {
            if (current_token.empty())
                return;
            token_list.emplace_back(std::move(current_token));
            current_token.clear();
        };

    for (size_t pos = 0; pos < line.size(); ++pos) {
        const char current_char{line[pos]};
        switch (state) {
        case State::InSpace:
            if (std::isspace(static_cast<unsigned char>(current_char))) {
                continue;
            }
            if (current_char == '\'') {
                state = State::InSingleQuote;
                continue;
            }
            if (current_char == '"') {
                state = State::InDoubleQuote;
                continue;
            }
            current_token.push_back(current_char);
            state = State::InUnquoted;
            continue;
        case State::InUnquoted:
            if (std::isspace(static_cast<unsigned char>(current_char))) {
                flush_token();
                state = State::InSpace;
                continue;
            }
            current_token.push_back(current_char);
            continue;
        case State::InSingleQuote:
            if (current_char == '\'') {
                flush_token();
                state = State::InSpace;
                continue;
            }
            current_token.push_back(current_char);
            continue;
        case State::InDoubleQuote:
            if (current_char == '"') {
                flush_token();
                state = State::InSpace;
                continue;
            }
            current_token.push_back(current_char);
            continue;
        }
    }
    flush_token();
    return token_list;
}

std::string BuildMmCifTokenPreview(
    const std::vector<std::string>& token_list,
    size_t max_items = 8) {
    std::ostringstream oss;
    oss << "[";
    for (size_t i = 0; i < token_list.size() && i < max_items; ++i) {
        if (i > 0)
            oss << ", ";
        oss << token_list[i];
    }
    if (token_list.size() > max_items)
        oss << ", ...";
    oss << "]";
    return oss.str();
}

void ResetReadState(CifFormatState& state) {
    state.data_block = std::make_unique<ModelImportState>();
    state.find_chemical_component_entry = false;
    state.find_component_atom_entry = false;
    state.find_component_bond_entry = false;
    state.loop_category_map.clear();
    state.data_item_map.clear();
    state.observed_component_atom_id_set_map.clear();
}

std::string ToLowerCopy(std::string value) {
    std::transform(
        value.begin(),
        value.end(),
        value.begin(),
        [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    return value;
}

bool ContainsCaseInsensitive(std::string haystack, std::string_view needle) {
    const auto lowered_haystack{ToLowerCopy(std::move(haystack))};
    const auto lowered_needle{ToLowerCopy(std::string{needle})};
    return lowered_haystack.find(lowered_needle) != std::string::npos;
}

std::pair<std::string, std::string> SplitBondIdAtomPair(const std::string& bond_id) {
    const auto pos{bond_id.find('_')};
    if (pos == std::string::npos) {
        return {"", ""};
    }
    return {bond_id.substr(0, pos), bond_id.substr(pos + 1)};
}

std::optional<Residue> ResolveObservedPolymerTemplateResidue(
    const ChemicalComponentEntry& component_entry) {
    const auto residue{ChemicalDataHelper::GetResidueFromString(component_entry.GetComponentId(), false)};
    if (residue != Residue::UNK) {
        return residue;
    }

    if (!ContainsCaseInsensitive(component_entry.GetComponentType(), "peptide linking")) {
        return std::nullopt;
    }

    static const std::unordered_map<std::string_view, Residue> kModifiedPeptideAliasMap{
        {"CSX", Residue::CYS},
    };
    const auto iter{kModifiedPeptideAliasMap.find(component_entry.GetComponentId())};
    if (iter == kModifiedPeptideAliasMap.end()) {
        return std::nullopt;
    }
    return iter->second;
}

void BuildObservedPolymerComponentBondEntry(CifFormatState& state) {
    auto* bond_key_system{state.data_block->GetBondKeySystemPtr()};
    for (const auto& [comp_id, observed_atom_id_set] : state.observed_component_atom_id_set_map) {
        if (observed_atom_id_set.empty()) {
            continue;
        }

        auto component_key{state.data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id)};
        auto* component_entry{state.data_block->GetChemicalComponentEntryPtr(component_key)};
        if (component_entry == nullptr) {
            continue;
        }

        const auto template_residue{ResolveObservedPolymerTemplateResidue(*component_entry)};
        if (!template_residue.has_value()) {
            continue;
        }

        for (const auto link : ComponentHelper::GetLinkList(*template_residue)) {
            const auto bond_id{ChemicalDataHelper::GetLabel(link)};
            const auto [atom_id_1, atom_id_2]{SplitBondIdAtomPair(bond_id)};
            if (atom_id_1.empty() || atom_id_2.empty()) {
                continue;
            }
            if (observed_atom_id_set.find(atom_id_1) == observed_atom_id_set.end() ||
                observed_atom_id_set.find(atom_id_2) == observed_atom_id_set.end()) {
                continue;
            }

            bond_key_system->RegisterBond(bond_id);
            const auto bond_key{bond_key_system->GetBondKey(bond_id)};
            if (component_entry->HasComponentBondEntry(bond_key)) {
                continue;
            }

            ComponentBondEntry bond_entry;
            bond_entry.bond_id = bond_id;
            bond_entry.bond_type = BondType::COVALENT;
            bond_entry.bond_order = BondOrder::UNK;
            bond_entry.aromatic_atom_flag = false;
            bond_entry.stereo_config = StereoChemistry::NONE;
            component_entry->AddComponentBondEntry(bond_key, bond_entry);
        }
    }
}

std::optional<std::string> GetFirstDataItemValue(
    const CifFormatState& state,
    std::string_view key) {
    auto iter{state.data_item_map.find(std::string{key})};
    if (iter == state.data_item_map.end())
        return std::nullopt;
    if (iter->second.empty())
        return std::nullopt;
    return iter->second.front();
}

void ParseMmCifDocument(std::istream& stream, const std::string& source_name, CifFormatState& state) {
    if (!stream) {
        Logger::Log(LogLevel::Error, "Cannot open the file: " + source_name);
        throw std::runtime_error("ParseMmCifDocument failed!");
    }

    std::vector<std::string> line_list;
    std::string line;
    while (std::getline(stream, line)) {
        StringHelper::StripCarriageReturn(line);
        line_list.emplace_back(std::move(line));
    }

    ResetReadState(state);

    size_t line_idx{0};
    while (line_idx < line_list.size()) {
        const auto& raw_line{line_list[line_idx]};
        const auto trimmed_line{TrimLeft(raw_line)};
        if (trimmed_line.empty() || trimmed_line.front() == '#') {
            ++line_idx;
            continue;
        }

        if (trimmed_line == "loop_") {
            ++line_idx;
            std::vector<std::string> column_name_list;
            while (line_idx < line_list.size()) {
                const auto& header_raw_line{line_list[line_idx]};
                const auto header_trimmed_line{TrimLeft(header_raw_line)};
                if (header_trimmed_line.empty()) {
                    ++line_idx;
                    continue;
                }
                if (header_trimmed_line.front() != '_')
                    break;
                auto header_token_list{SplitMmCifTokens(std::string{header_trimmed_line})};
                if (!header_token_list.empty()) {
                    column_name_list.emplace_back(header_token_list.front());
                }
                ++line_idx;
            }
            if (column_name_list.empty())
                continue;

            CifFormatState::ParsedLoopCategory loop_category;
            loop_category.column_name_list = column_name_list;
            const auto expected_column_size{loop_category.column_name_list.size()};
            std::vector<std::string> row_token_list;
            size_t row_start_line_number{line_idx + 1};

            while (line_idx < line_list.size()) {
                const auto& row_raw_line{line_list[line_idx]};
                const auto row_trimmed_line{TrimLeft(row_raw_line)};
                if (row_trimmed_line.empty()) {
                    ++line_idx;
                    continue;
                }
                if (row_trimmed_line.front() == '#') {
                    ++line_idx;
                    break;
                }
                if (row_token_list.empty() &&
                    (row_trimmed_line == "loop_" ||
                     StartsWithToken(row_trimmed_line, "data_") ||
                     row_trimmed_line.front() == '_')) {
                    break;
                }
                if (row_token_list.empty())
                    row_start_line_number = line_idx + 1;

                if (!row_raw_line.empty() && row_raw_line.front() == ';') {
                    std::string multiline_value{row_raw_line.substr(1)};
                    ++line_idx;
                    auto terminated{false};
                    while (line_idx < line_list.size()) {
                        const auto& multiline_raw_line{line_list[line_idx]};
                        if (!multiline_raw_line.empty() && multiline_raw_line.front() == ';') {
                            terminated = true;
                            break;
                        }
                        if (!multiline_value.empty())
                            multiline_value += "\n";
                        multiline_value += multiline_raw_line;
                        ++line_idx;
                    }
                    if (!terminated) {
                        Logger::Log(
                            LogLevel::Warning,
                            "ParseMmCifDocument() unterminated multiline token in loop category " + BuildLoopCategoryPrefix(loop_category.column_name_list.front()) + " near file line " + std::to_string(row_start_line_number) + ".");
                        break;
                    }
                    row_token_list.emplace_back(std::move(multiline_value));
                    ++line_idx;
                } else {
                    auto token_list{SplitMmCifTokens(std::string{row_trimmed_line})};
                    row_token_list.insert(row_token_list.end(), token_list.begin(), token_list.end());
                    ++line_idx;
                }

                while (row_token_list.size() >= expected_column_size) {
                    std::vector<std::string> row_value_list(
                        row_token_list.begin(),
                        row_token_list.begin() + static_cast<std::ptrdiff_t>(expected_column_size));
                    loop_category.row_list.emplace_back(CifFormatState::ParsedLoopRow{
                        std::move(row_value_list),
                        row_start_line_number});
                    row_token_list.erase(
                        row_token_list.begin(),
                        row_token_list.begin() + static_cast<std::ptrdiff_t>(expected_column_size));
                    row_start_line_number = line_idx + 1;
                }
            }

            if (!row_token_list.empty()) {
                Logger::Log(
                    LogLevel::Warning,
                    "ParseMmCifDocument() loop category " + BuildLoopCategoryPrefix(loop_category.column_name_list.front()) + " ends with incomplete row (token count = " + std::to_string(row_token_list.size()) + ", expected = " + std::to_string(expected_column_size) + ").");
            }

            auto category_prefix{BuildLoopCategoryPrefix(loop_category.column_name_list.front())};
            state.loop_category_map[category_prefix].emplace_back(std::move(loop_category));
            continue;
        }

        if (trimmed_line.front() == '_') {
            auto token_list{SplitMmCifTokens(std::string{trimmed_line})};
            if (token_list.empty()) {
                ++line_idx;
                continue;
            }

            std::string key{token_list.front()};
            std::string value;
            if (token_list.size() >= 2) {
                value = token_list[1];
                state.data_item_map[key].emplace_back(std::move(value));
                ++line_idx;
                continue;
            }

            ++line_idx;
            while (line_idx < line_list.size()) {
                const auto& value_raw_line{line_list[line_idx]};
                const auto value_trimmed_line{TrimLeft(value_raw_line)};
                if (value_trimmed_line.empty()) {
                    ++line_idx;
                    continue;
                }
                if (value_trimmed_line.front() == '#')
                    break;
                if (!value_raw_line.empty() && value_raw_line.front() == ';') {
                    value = value_raw_line.substr(1);
                    ++line_idx;
                    while (line_idx < line_list.size()) {
                        const auto& multiline_raw_line{line_list[line_idx]};
                        if (!multiline_raw_line.empty() && multiline_raw_line.front() == ';')
                            break;
                        if (!value.empty())
                            value += "\n";
                        value += multiline_raw_line;
                        ++line_idx;
                    }
                    break;
                }
                auto value_token_list{SplitMmCifTokens(std::string{value_trimmed_line})};
                value = value_token_list.empty() ? std::string{} : value_token_list.front();
                break;
            }
            state.data_item_map[key].emplace_back(std::move(value));
            ++line_idx;
            continue;
        }

        ++line_idx;
    }
}

void ParseLoopBlock(
    const CifFormatState& state,
    std::string_view data_block_prefix,
    const CifLoopRowHandler& row_handler) {
    auto category_iter{state.loop_category_map.find(std::string{data_block_prefix})};
    if (category_iter == state.loop_category_map.end())
        return;

    for (const auto& category : category_iter->second) {
        CifColumnIndexMap column_index_map;
        for (size_t i = 0; i < category.column_name_list.size(); ++i) {
            std::string short_name{category.column_name_list[i]};
            if (short_name.rfind(data_block_prefix, 0) == 0) {
                short_name = short_name.substr(data_block_prefix.size());
            }
            column_index_map[short_name] = i;
        }

        size_t loop_row_number{0};
        for (const auto& row : category.row_list) {
            ++loop_row_number;
            if (row.token_list.size() != column_index_map.size()) {
                Logger::Log(
                    LogLevel::Warning,
                    "ParseLoopBlock(" + std::string(data_block_prefix) + ") row " + std::to_string(loop_row_number) + " has token count " + std::to_string(row.token_list.size()) + " but expects " + std::to_string(column_index_map.size()) + " at file line " + std::to_string(row.line_number) + ". tokens = " + BuildMmCifTokenPreview(row.token_list));
                continue;
            }
            try {
                row_handler(column_index_map, row.token_list);
            } catch (const std::exception& ex) {
                Logger::Log(
                    LogLevel::Warning,
                    "ParseLoopBlock(" + std::string(data_block_prefix) + ") row " + std::to_string(loop_row_number) + " failed at file line " + std::to_string(row.line_number) + ": " + std::string(ex.what()) + ". tokens = " + BuildMmCifTokenPreview(row.token_list));
            }
        }
    }
}

void LogHeaderSummary(const ModelImportState& data_block) {
    std::ostringstream oss;
    oss << "CIF Header Information:\n";
    oss << "#Entities = " << data_block.GetEntityTypeMap().size() << "\n";
    for (const auto& [entity_id, chain_id] : data_block.GetChainIDListMap()) {
        oss << "[" << entity_id << "] : ";
        for (size_t i = 0; i < chain_id.size(); ++i) {
            oss << chain_id.at(i);
            if (i < chain_id.size() - 1)
                oss << ",";
        }
        oss << "\n";
    }

    const auto element_size{data_block.GetElementTypeList().size()};
    oss << "#Elementry types = " << element_size << "\n";
    oss << "Element type list : ";
    size_t count{0};
    for (const auto element : data_block.GetElementTypeList()) {
        oss << ChemicalDataHelper::GetLabel(element);
        if (count < element_size - 1)
            oss << ",";
        ++count;
    }
    oss << "\n";
    Logger::Log(LogLevel::Info, oss.str());
}

void BuildDefaultChemicalComponentEntry(CifFormatState& state, const std::string& comp_id) {
    auto entry{std::make_unique<ChemicalComponentEntry>()};
    entry->SetComponentId(comp_id);
    entry->SetComponentName("");
    entry->SetComponentType("");
    entry->SetComponentFormula("");
    entry->SetComponentMolecularWeight(0.0f);
    auto standard_flag{false};
    if (ChemicalDataHelper::GetResidueFromString(comp_id) != Residue::UNK)
        standard_flag = true;
    entry->SetStandardMonomerFlag(standard_flag);

    state.data_block->GetComponentKeySystemPtr()->RegisterComponent(comp_id);
    auto component_key{state.data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id)};
    state.data_block->AddChemicalComponentEntry(component_key, std::move(entry));
}

void BuildDefaultComponentAtomEntry(
    CifFormatState& state,
    const std::string& comp_id,
    const std::string& atom_id,
    const std::string& element_symbol) {
    auto component_key{state.data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id)};

    state.data_block->GetAtomKeySystemPtr()->RegisterAtom(atom_id);
    auto atom_key{state.data_block->GetAtomKeySystemPtr()->GetAtomKey(atom_id)};

    ComponentAtomEntry atom_entry;
    atom_entry.atom_id = atom_id;
    atom_entry.element_type = ChemicalDataHelper::GetElementFromString(element_symbol);
    atom_entry.aromatic_atom_flag = false;
    atom_entry.stereo_config = StereoChemistry::NONE;

    state.data_block->AddComponentAtomEntry(component_key, atom_key, atom_entry);
}

void BuildDefaultComponentBondEntry(CifFormatState& state) {
    for (auto& residue : ChemicalDataHelper::GetStandardAminoAcidList()) {
        auto comp_id{ChemicalDataHelper::GetLabel(residue)};
        auto component_key{state.data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id)};
        for (auto& link : ComponentHelper::GetLinkList(residue)) {
            auto bond_id{ChemicalDataHelper::GetLabel(link)};
            state.data_block->GetBondKeySystemPtr()->RegisterBond(bond_id);
            auto bond_key{state.data_block->GetBondKeySystemPtr()->GetBondKey(bond_id)};

            ComponentBondEntry bond_entry;
            bond_entry.bond_id = bond_id;
            bond_entry.bond_type = BondType::COVALENT;
            bond_entry.bond_order = BondOrder::UNK;
            bond_entry.aromatic_atom_flag = false;
            bond_entry.stereo_config = StereoChemistry::NONE;

            state.data_block->AddComponentBondEntry(component_key, bond_key, bond_entry);
        }
    }
}

void BuildPepetideBondEntry(CifFormatState& state) {
    for (auto& residue : ChemicalDataHelper::GetStandardAminoAcidList()) {
        auto comp_id{ChemicalDataHelper::GetLabel(residue)};
        auto component_key{state.data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id)};
        if (!state.data_block->HasChemicalComponentEntry(component_key))
            continue;
        auto bond_id{ChemicalDataHelper::GetLabel(Link::C_N)};
        state.data_block->GetBondKeySystemPtr()->RegisterBond(bond_id);
        auto bond_key{state.data_block->GetBondKeySystemPtr()->GetBondKey(bond_id)};

        ComponentBondEntry bond_entry;
        bond_entry.bond_id = bond_id;
        bond_entry.bond_type = BondType::COVALENT;
        bond_entry.bond_order = BondOrder::SINGLE;
        bond_entry.aromatic_atom_flag = false;
        bond_entry.stereo_config = StereoChemistry::NONE;

        state.data_block->AddComponentBondEntry(component_key, bond_key, bond_entry);
    }
}

void BuildPhosphodiesterBondEntry(CifFormatState& state) {
    for (auto& residue : ChemicalDataHelper::GetStandardNucleotideList()) {
        auto comp_id{ChemicalDataHelper::GetLabel(residue)};
        auto component_key{state.data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id)};
        if (!state.data_block->HasChemicalComponentEntry(component_key))
            continue;
        auto bond_id{ChemicalDataHelper::GetLabel(Link::P_O3p)};
        state.data_block->GetBondKeySystemPtr()->RegisterBond(bond_id);
        auto bond_key{state.data_block->GetBondKeySystemPtr()->GetBondKey(bond_id)};

        ComponentBondEntry bond_entry;
        bond_entry.bond_id = bond_id;
        bond_entry.bond_type = BondType::COVALENT;
        bond_entry.bond_order = BondOrder::SINGLE;
        bond_entry.aromatic_atom_flag = false;
        bond_entry.stereo_config = StereoChemistry::NONE;

        state.data_block->AddComponentBondEntry(component_key, bond_key, bond_entry);
    }
}

void LoadChemicalComponentAtomBlock(CifFormatState& state);
void LoadChemicalComponentBondBlock(CifFormatState& state);

void LoadChemicalComponentBlock(CifFormatState& state) {
    ParseLoopBlock(
        state,
        "_chem_comp.",
        [&state](const CifColumnIndexMap& index_map, const std::vector<std::string>& token_list) {
            auto comp_id{token_list[index_map.at("id")]};
            auto name{token_list[index_map.at("name")]};
            auto type{token_list[index_map.at("type")]};
            auto formula{token_list[index_map.at("formula")]};
            auto formula_weight_str{token_list[index_map.at("formula_weight")]};
            auto standard_flag_str{token_list[index_map.at("mon_nstd_flag")]};
            auto formula_weight{0.0f};
            try {
                formula_weight = std::stof(formula_weight_str);
            } catch (const std::exception&) {
                formula_weight = 0.0f;
            }
            auto entry{std::make_unique<ChemicalComponentEntry>()};
            entry->SetComponentId(comp_id);
            entry->SetComponentName(name);
            entry->SetComponentType(type);
            entry->SetComponentFormula(formula);
            entry->SetComponentMolecularWeight(formula_weight);
            auto standard_flag{false};
            if (standard_flag_str == "n" || standard_flag_str == "no")
                standard_flag = false;
            else if (standard_flag_str == "y" || standard_flag_str == "yes")
                standard_flag = true;
            entry->SetStandardMonomerFlag(standard_flag);

            state.data_block->GetComponentKeySystemPtr()->RegisterComponent(comp_id);
            auto component_key{state.data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id)};
            state.data_block->AddChemicalComponentEntry(component_key, std::move(entry));

            state.find_chemical_component_entry = true;
        });

    LoadChemicalComponentAtomBlock(state);
    LoadChemicalComponentBondBlock(state);
}

void LoadChemicalComponentAtomBlock(CifFormatState& state) {
    ParseLoopBlock(
        state,
        "_chem_comp_atom.",
        [&state](const CifColumnIndexMap& index_map, const std::vector<std::string>& token_list) {
            auto comp_id{token_list[index_map.at("comp_id")]};
            auto atom_id{token_list[index_map.at("atom_id")]};
            auto element_symbol{token_list[index_map.at("type_symbol")]};
            auto pdbx_aromatic_flag{token_list[index_map.at("pdbx_aromatic_flag")]};
            auto pdbx_stereo_config{token_list[index_map.at("pdbx_stereo_config")]};

            StringHelper::EraseCharFromString(atom_id, '\"');

            ComponentAtomEntry atom_entry;
            atom_entry.atom_id = atom_id;
            atom_entry.element_type = ChemicalDataHelper::GetElementFromString(element_symbol);
            atom_entry.aromatic_atom_flag = (pdbx_aromatic_flag == "Y");
            atom_entry.stereo_config =
                ChemicalDataHelper::GetStereoChemistryFromString(pdbx_stereo_config);

            state.data_block->GetAtomKeySystemPtr()->RegisterAtom(atom_id);
            auto component_key{state.data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id)};
            auto atom_key{state.data_block->GetAtomKeySystemPtr()->GetAtomKey(atom_id)};
            state.data_block->AddComponentAtomEntry(component_key, atom_key, atom_entry);

            state.find_component_atom_entry = true;
        });
}

void LoadChemicalComponentBondBlock(CifFormatState& state) {
    ParseLoopBlock(
        state,
        "_chem_comp_bond.",
        [&state](const CifColumnIndexMap& index_map, const std::vector<std::string>& token_list) {
            auto comp_id{token_list[index_map.at("comp_id")]};
            auto atom_id_1{token_list[index_map.at("atom_id_1")]};
            auto atom_id_2{token_list[index_map.at("atom_id_2")]};
            auto bond_order{token_list[index_map.at("value_order")]};
            auto pdbx_aromatic_flag{token_list[index_map.at("pdbx_aromatic_flag")]};
            auto pdbx_stereo_config{token_list[index_map.at("pdbx_stereo_config")]};

            StringHelper::EraseCharFromString(atom_id_1, '\"');
            StringHelper::EraseCharFromString(atom_id_2, '\"');

            state.data_block->GetBondKeySystemPtr()->RegisterBond(atom_id_1, atom_id_2);
            auto component_key{state.data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id)};
            auto bond_key{state.data_block->GetBondKeySystemPtr()->GetBondKey(atom_id_1, atom_id_2)};
            auto bond_id{state.data_block->GetBondKeySystemPtr()->GetBondId(bond_key)};

            ComponentBondEntry bond_entry;
            bond_entry.bond_id = bond_id;
            bond_entry.bond_type = BondType::COVALENT;
            bond_entry.bond_order = ChemicalDataHelper::GetBondOrderFromString(bond_order);
            bond_entry.aromatic_atom_flag = (pdbx_aromatic_flag == "Y");
            bond_entry.stereo_config =
                ChemicalDataHelper::GetStereoChemistryFromString(pdbx_stereo_config);

            state.data_block->AddComponentBondEntry(component_key, bond_key, bond_entry);
            state.find_component_bond_entry = true;
        });
}

void LoadDatabaseBlock(CifFormatState& state) {
    std::unordered_map<std::string, std::string> data_map;
    ParseLoopBlock(
        state,
        "_database_2.",
        [&data_map](const CifColumnIndexMap& index_map, const std::vector<std::string>& token_list) {
            auto database_id{GetTokenOptional(index_map, token_list, {"database_id"})};
            auto database_code{GetTokenOptional(index_map, token_list, {"database_code"})};
            if (!database_id.has_value() || !database_code.has_value())
                return;
            if (IsMmCifMissingValue(*database_id) || IsMmCifMissingValue(*database_code))
                return;
            data_map[*database_id] = *database_code;
        });
    if (data_map.find("PDB") != data_map.end()) {
        state.data_block->SetPdbID(data_map.at("PDB"));
    }
    if (data_map.find("EMDB") != data_map.end()) {
        state.data_block->SetEmdID(data_map.at("EMDB"));
    } else {
        state.data_block->SetEmdID("X-RAY DIFF");
    }
}

void LoadEntityBlock(CifFormatState& state) {
    ParseLoopBlock(
        state,
        "_entity.",
        [&state](const CifColumnIndexMap& index_map, const std::vector<std::string>& token_list) {
            if (index_map.find("id") == index_map.end() || index_map.find("type") == index_map.end()) {
                return;
            }
            if (index_map.at("id") >= token_list.size() || index_map.at("type") >= token_list.size()) {
                return;
            }

            auto entity_id{token_list[index_map.at("id")]};
            auto entity_type{token_list[index_map.at("type")]};
            std::string molecules_size_string{"-1"};
            if (index_map.find("pdbx_number_of_molecules") != index_map.end() &&
                index_map.at("pdbx_number_of_molecules") < token_list.size()) {
                molecules_size_string = token_list[index_map.at("pdbx_number_of_molecules")];
            }
            int molecules_size{-1};
            try {
                molecules_size = std::stoi(molecules_size_string);
            } catch (const std::exception&) {
                Logger::Log(LogLevel::Error, "Invalid molecules size: " + molecules_size_string);
            }

            state.data_block->AddEntityTypeInEntityMap(
                entity_id, ChemicalDataHelper::GetEntityFromString(entity_type));
            state.data_block->AddMoleculesSizeInEntityMap(entity_id, molecules_size);
        });

    ParseLoopBlock(
        state,
        "_struct_asym.",
        [&state](const CifColumnIndexMap& index_map, const std::vector<std::string>& token_list) {
            if (index_map.find("entity_id") == index_map.end() || index_map.find("id") == index_map.end()) {
                return;
            }
            if (index_map.at("entity_id") >= token_list.size() || index_map.at("id") >= token_list.size()) {
                return;
            }
            auto entity_id{token_list[index_map.at("entity_id")]};
            auto chain_id{token_list[index_map.at("id")]};
            state.data_block->AddChainIDInEntityMap(entity_id, chain_id);
        });

    const auto need_entity_fallback{state.data_block->GetEntityTypeMap().empty()};
    const auto need_chain_fallback{state.data_block->GetChainIDListMap().empty()};
    if (!need_entity_fallback && !need_chain_fallback)
        return;

    Logger::Log(
        LogLevel::Warning,
        "Cannot parse complete entity/chain metadata from loop blocks. "
        "Trying key-value fallback for _entity/_struct_asym.");

    std::string entity_id{GetFirstDataItemValue(state, "_entity.id").value_or("")};
    std::string entity_type{GetFirstDataItemValue(state, "_entity.type").value_or("")};
    std::string molecules_size_raw{
        GetFirstDataItemValue(state, "_entity.pdbx_number_of_molecules").value_or("")};
    std::string struct_asym_id{GetFirstDataItemValue(state, "_struct_asym.id").value_or("")};
    std::string struct_asym_entity_id{
        GetFirstDataItemValue(state, "_struct_asym.entity_id").value_or("")};

    if (need_entity_fallback) {
        if (IsMmCifMissingValue(entity_id) || IsMmCifMissingValue(entity_type)) {
            Logger::Log(
                LogLevel::Warning,
                "Key-value fallback cannot recover required _entity.id/_entity.type fields.");
        } else {
            int molecules_size{-1};
            if (!IsMmCifMissingValue(molecules_size_raw)) {
                try {
                    molecules_size = std::stoi(molecules_size_raw);
                } catch (const std::exception&) {
                    Logger::Log(
                        LogLevel::Warning,
                        "Invalid _entity.pdbx_number_of_molecules in key-value fallback: " + molecules_size_raw);
                }
            }
            state.data_block->AddEntityTypeInEntityMap(
                entity_id, ChemicalDataHelper::GetEntityFromString(entity_type));
            state.data_block->AddMoleculesSizeInEntityMap(entity_id, molecules_size);
        }
    }

    if (need_chain_fallback) {
        if (IsMmCifMissingValue(struct_asym_id) || IsMmCifMissingValue(struct_asym_entity_id)) {
            Logger::Log(
                LogLevel::Warning,
                "Key-value fallback cannot recover required _struct_asym.id/_struct_asym.entity_id fields.");
        } else {
            state.data_block->AddChainIDInEntityMap(struct_asym_entity_id, struct_asym_id);
        }
    }
}

void LoadXRayResolutionInfo(CifFormatState& state) {
    auto found_resolution{false};
    auto found_resolution_method{false};

    if (auto resolution_value{GetFirstDataItemValue(state, "_refine.ls_d_res_high")}) {
        if (!IsMmCifMissingValue(*resolution_value)) {
            state.data_block->SetResolution(*resolution_value);
            found_resolution = true;
        }
    }
    if (auto method_value{GetFirstDataItemValue(state, "_refine.pdbx_refine_id")}) {
        if (!IsMmCifMissingValue(*method_value)) {
            state.data_block->SetResolutionMethod(*method_value);
            found_resolution_method = true;
        }
    }

    ParseLoopBlock(
        state,
        "_refine.",
        [&state, &found_resolution, &found_resolution_method](
            const CifColumnIndexMap& index_map,
            const std::vector<std::string>& token_list) {
            if (!found_resolution) {
                auto resolution_value{GetTokenOptional(index_map, token_list, {"ls_d_res_high"})};
                if (resolution_value.has_value() && !IsMmCifMissingValue(*resolution_value)) {
                    state.data_block->SetResolution(*resolution_value);
                    found_resolution = true;
                }
            }
            if (!found_resolution_method) {
                auto method_value{GetTokenOptional(index_map, token_list, {"pdbx_refine_id"})};
                if (method_value.has_value() && !IsMmCifMissingValue(*method_value)) {
                    state.data_block->SetResolutionMethod(*method_value);
                    found_resolution_method = true;
                }
            }
        });
}

void LoadPdbxData(CifFormatState& state) {
    auto found_resolution{false};
    auto found_resolution_method{false};

    if (auto resolution_value{GetFirstDataItemValue(state, "_em_3d_reconstruction.resolution")}) {
        if (!IsMmCifMissingValue(*resolution_value)) {
            state.data_block->SetResolution(*resolution_value);
            found_resolution = true;
        }
    }
    if (auto resolution_method_value{
            GetFirstDataItemValue(state, "_em_3d_reconstruction.resolution_method")}) {
        if (!IsMmCifMissingValue(*resolution_method_value)) {
            state.data_block->SetResolutionMethod(*resolution_method_value);
            found_resolution_method = true;
        }
    }

    ParseLoopBlock(
        state,
        "_em_3d_reconstruction.",
        [&state, &found_resolution, &found_resolution_method](
            const CifColumnIndexMap& index_map,
            const std::vector<std::string>& token_list) {
            if (!found_resolution) {
                auto resolution_value{GetTokenOptional(index_map, token_list, {"resolution"})};
                if (resolution_value.has_value() && !IsMmCifMissingValue(*resolution_value)) {
                    state.data_block->SetResolution(*resolution_value);
                    found_resolution = true;
                }
            }
            if (!found_resolution_method) {
                auto resolution_method_value{
                    GetTokenOptional(index_map, token_list, {"resolution_method"})};
                if (resolution_method_value.has_value() &&
                    !IsMmCifMissingValue(*resolution_method_value)) {
                    state.data_block->SetResolutionMethod(*resolution_method_value);
                    found_resolution_method = true;
                }
            }
        });

    if (!found_resolution || !found_resolution_method) {
        LoadXRayResolutionInfo(state);
    }
}

void LoadAtomTypeBlock(CifFormatState& state) {
    ParseLoopBlock(
        state,
        "_atom_type.",
        [&state](const CifColumnIndexMap& index_map, const std::vector<std::string>& token_list) {
            auto element_type_string{token_list[index_map.at("symbol")]};
            auto element{ChemicalDataHelper::GetElementFromString(element_type_string)};
            state.data_block->AddElementType(element);
        });
}

void LoadStructureConformationBlock(CifFormatState& state) {
    ParseLoopBlock(
        state,
        "_struct_conf.",
        [&state](const CifColumnIndexMap& index_map, const std::vector<std::string>& token_list) {
            auto helix_id_opt{GetTokenOptional(index_map, token_list, {"id"})};
            auto conf_type_opt{GetTokenOptional(index_map, token_list, {"conf_type_id"})};
            auto chain_id_beg_opt{GetTokenOptional(index_map, token_list, {"beg_label_asym_id"})};
            auto residue_id_beg_opt{GetTokenOptional(index_map, token_list, {"beg_label_seq_id"})};
            auto chain_id_end_opt{GetTokenOptional(index_map, token_list, {"end_label_asym_id"})};
            auto residue_id_end_opt{GetTokenOptional(index_map, token_list, {"end_label_seq_id"})};
            if (!helix_id_opt.has_value() || !conf_type_opt.has_value() ||
                !chain_id_beg_opt.has_value() || !residue_id_beg_opt.has_value() ||
                !chain_id_end_opt.has_value() || !residue_id_end_opt.has_value()) {
                return;
            }

            auto residue_id_beg{TryParseInt(*residue_id_beg_opt)};
            auto residue_id_end{TryParseInt(*residue_id_end_opt)};
            if (!residue_id_beg.has_value() || !residue_id_end.has_value()) {
                Logger::Log(
                    LogLevel::Warning,
                    "LoadStructureConformationBlock(): invalid residue range in _struct_conf. "
                    "row is skipped.");
                return;
            }

            HelixRange range;
            range.chain_id_beg = *chain_id_beg_opt;
            range.seq_id_beg = *residue_id_beg;
            range.chain_id_end = *chain_id_end_opt;
            range.seq_id_end = *residue_id_end;
            range.conf_type = *conf_type_opt;
            state.data_block->AddHelixRange(*helix_id_opt, range);
        });
}

void LoadStructureConnectionBlock(CifFormatState& state) {
    ParseLoopBlock(
        state,
        "_struct_conn.",
        [&state](const CifColumnIndexMap& index_map, const std::vector<std::string>& token_list) {
            auto conn_id{GetTokenOptional(index_map, token_list, {"id"})};
            auto conn_type_id{GetTokenOptional(index_map, token_list, {"conn_type_id"})};
            auto ptnr1_asym_id{GetTokenOptional(index_map, token_list, {"ptnr1_label_asym_id"})};
            auto ptnr1_comp_id{GetTokenOptional(index_map, token_list, {"ptnr1_label_comp_id"})};
            auto ptnr1_seq_id{GetTokenOptional(index_map, token_list, {"ptnr1_label_seq_id"})};
            auto ptnr1_atom_id{GetTokenOptional(index_map, token_list, {"ptnr1_label_atom_id"})};
            auto ptnr2_asym_id{GetTokenOptional(index_map, token_list, {"ptnr2_label_asym_id"})};
            auto ptnr2_comp_id{GetTokenOptional(index_map, token_list, {"ptnr2_label_comp_id"})};
            auto ptnr2_seq_id{GetTokenOptional(index_map, token_list, {"ptnr2_label_seq_id"})};
            auto ptnr2_atom_id{GetTokenOptional(index_map, token_list, {"ptnr2_label_atom_id"})};
            auto value_order_str{
                GetTokenOptional(index_map, token_list, {"pdbx_value_order"})};

            if (!conn_type_id.has_value() ||
                !ptnr1_asym_id.has_value() || !ptnr1_comp_id.has_value() ||
                !ptnr1_seq_id.has_value() || !ptnr1_atom_id.has_value() ||
                !ptnr2_asym_id.has_value() || !ptnr2_comp_id.has_value() ||
                !ptnr2_seq_id.has_value() || !ptnr2_atom_id.has_value()) {
                return;
            }
            auto conn_id_label{conn_id.value_or("UNKNOWN_CONN")};
            auto ptnr1_atom_id_value{*ptnr1_atom_id};
            auto ptnr2_atom_id_value{*ptnr2_atom_id};
            auto ptnr1_seq_id_value{IsMmCifMissingValue(*ptnr1_seq_id) ? "." : *ptnr1_seq_id};
            auto ptnr2_seq_id_value{IsMmCifMissingValue(*ptnr2_seq_id) ? "." : *ptnr2_seq_id};

            StringHelper::EraseCharFromString(ptnr1_atom_id_value, '\"');
            StringHelper::EraseCharFromString(ptnr2_atom_id_value, '\"');

            state.data_block->GetBondKeySystemPtr()->RegisterBond(
                ptnr1_atom_id_value, ptnr2_atom_id_value);
            auto component_key{
                state.data_block->GetComponentKeySystemPtr()->GetComponentKey(*ptnr1_comp_id)};
            auto bond_key{
                state.data_block->GetBondKeySystemPtr()->GetBondKey(
                    ptnr1_atom_id_value, ptnr2_atom_id_value)};
            auto bond_id{state.data_block->GetBondKeySystemPtr()->GetBondId(bond_key)};
            auto bond_type{ChemicalDataHelper::GetBondTypeFromString(*conn_type_id)};
            if (bond_type == BondType::HYDROGEN)
                return;

            ComponentBondEntry bond_entry;
            bond_entry.bond_id = bond_id;
            bond_entry.bond_type = bond_type;
            if (!value_order_str.has_value() || IsMmCifMissingValue(*value_order_str)) {
                bond_entry.bond_order = BondOrder::UNK;
            } else {
                bond_entry.bond_order = ChemicalDataHelper::GetBondOrderFromString(*value_order_str);
            }
            bond_entry.aromatic_atom_flag = false;
            bond_entry.stereo_config = StereoChemistry::NONE;
            state.data_block->AddComponentBondEntry(component_key, bond_key, bond_entry);

            auto ptnr1_model_num_token{GetTokenOptional(
                index_map,
                token_list,
                {"pdbx_ptnr1_PDB_model_num", "ptnr1_PDB_model_num"})};
            auto ptnr2_model_num_token{GetTokenOptional(
                index_map,
                token_list,
                {"pdbx_ptnr2_PDB_model_num", "ptnr2_PDB_model_num"})};
            auto model_number_1{
                ParseIntOrDefault(
                    ptnr1_model_num_token.value_or("1"),
                    1,
                    "_struct_conn.pdbx_ptnr1_PDB_model_num",
                    "LoadStructureConnectionBlock()[" + conn_id_label + "]")};
            auto model_number_2{
                ParseIntOrDefault(
                    ptnr2_model_num_token.value_or(std::to_string(model_number_1)),
                    model_number_1,
                    "_struct_conn.pdbx_ptnr2_PDB_model_num",
                    "LoadStructureConnectionBlock()[" + conn_id_label + "]")};

            auto atom_1{state.data_block->GetAtomObjectPtrInTuple(
                model_number_1,
                *ptnr1_asym_id,
                *ptnr1_comp_id,
                ptnr1_seq_id_value,
                ptnr1_atom_id_value)};
            auto atom_2{state.data_block->GetAtomObjectPtrInTuple(
                model_number_2,
                *ptnr2_asym_id,
                *ptnr2_comp_id,
                ptnr2_seq_id_value,
                ptnr2_atom_id_value)};
            if (atom_1 == nullptr) {
                atom_1 = state.data_block->GetAtomObjectPtrInAnyModel(
                    *ptnr1_asym_id,
                    *ptnr1_comp_id,
                    ptnr1_seq_id_value,
                    ptnr1_atom_id_value,
                    &model_number_1);
            }
            if (atom_2 == nullptr) {
                atom_2 = state.data_block->GetAtomObjectPtrInAnyModel(
                    *ptnr2_asym_id,
                    *ptnr2_comp_id,
                    ptnr2_seq_id_value,
                    ptnr2_atom_id_value,
                    &model_number_2);
            }
            if (atom_1 == nullptr || atom_2 == nullptr) {
                Logger::Log(LogLevel::Warning, "Cannot find atom object for connection: " + conn_id_label);
                return;
            }
            if (model_number_1 != model_number_2) {
                Logger::Log(
                    LogLevel::Warning,
                    "Skip cross-model connection [" + conn_id_label + "] between model " + std::to_string(model_number_1) + " and model " + std::to_string(model_number_2) + ".");
                return;
            }

            auto bond_object{std::make_unique<BondObject>(atom_1, atom_2)};
            bond_object->SetBondKey(bond_key);
            bond_object->SetBondType(bond_entry.bond_type);
            bond_object->SetBondOrder(bond_entry.bond_order);
            bond_object->SetSpecialBondFlag(atom_1->GetSpecialAtomFlag());
            state.data_block->AddBondObject(std::move(bond_object));
        });
}

void LoadStructureSheetBlock(CifFormatState& state) {
    ParseLoopBlock(
        state,
        "_struct_sheet.",
        [&state](const CifColumnIndexMap& index_map, const std::vector<std::string>& token_list) {
            auto sheet_id_opt{GetTokenOptional(index_map, token_list, {"id"})};
            auto strands_size_opt{GetTokenOptional(index_map, token_list, {"number_strands"})};
            if (!sheet_id_opt.has_value() || !strands_size_opt.has_value())
                return;
            auto strands_size{
                ParseIntOrDefault(
                    *strands_size_opt,
                    0,
                    "_struct_sheet.number_strands",
                    "LoadStructureSheetBlock()")};
            if (strands_size <= 0)
                return;
            state.data_block->AddSheetStrands(*sheet_id_opt, strands_size);
        });

    ParseLoopBlock(
        state,
        "_struct_sheet_range.",
        [&state](const CifColumnIndexMap& index_map, const std::vector<std::string>& token_list) {
            auto sheet_id_opt{GetTokenOptional(index_map, token_list, {"sheet_id"})};
            auto range_id_opt{GetTokenOptional(index_map, token_list, {"id"})};
            auto chain_id_beg_opt{GetTokenOptional(index_map, token_list, {"beg_label_asym_id"})};
            auto residue_id_beg_opt{GetTokenOptional(index_map, token_list, {"beg_label_seq_id"})};
            auto chain_id_end_opt{GetTokenOptional(index_map, token_list, {"end_label_asym_id"})};
            auto residue_id_end_opt{GetTokenOptional(index_map, token_list, {"end_label_seq_id"})};
            if (!sheet_id_opt.has_value() || !range_id_opt.has_value() ||
                !chain_id_beg_opt.has_value() || !residue_id_beg_opt.has_value() ||
                !chain_id_end_opt.has_value() || !residue_id_end_opt.has_value()) {
                return;
            }

            auto residue_id_beg{TryParseInt(*residue_id_beg_opt)};
            auto residue_id_end{TryParseInt(*residue_id_end_opt)};
            if (!residue_id_beg.has_value() || !residue_id_end.has_value()) {
                Logger::Log(
                    LogLevel::Warning,
                    "LoadStructureSheetBlock(): invalid residue range in _struct_sheet_range. "
                    "row is skipped.");
                return;
            }

            auto composite_sheet_id{*sheet_id_opt + *range_id_opt};
            SheetRange range;
            range.chain_id_beg = *chain_id_beg_opt;
            range.seq_id_beg = *residue_id_beg;
            range.chain_id_end = *chain_id_end_opt;
            range.seq_id_end = *residue_id_end;
            state.data_block->AddSheetRange(composite_sheet_id, range);
        });
}

struct AtomAltLocKey {
    int model_number;
    std::string chain_id;
    std::string comp_id;
    std::string sequence_id_token;
    std::string atom_id;

    bool operator==(const AtomAltLocKey& other) const {
        return model_number == other.model_number &&
               chain_id == other.chain_id &&
               comp_id == other.comp_id &&
               sequence_id_token == other.sequence_id_token &&
               atom_id == other.atom_id;
    }
};

struct AtomAltLocKeyHash {
    size_t operator()(const AtomAltLocKey& key) const {
        size_t seed{std::hash<int>{}(key.model_number)};
        seed ^= std::hash<std::string>{}(key.chain_id) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= std::hash<std::string>{}(key.comp_id) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= std::hash<std::string>{}(key.sequence_id_token) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= std::hash<std::string>{}(key.atom_id) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

void LoadAtomSiteBlock(CifFormatState& state) {
    std::unordered_map<AtomAltLocKey, AtomObject*, AtomAltLocKeyHash> altloc_primary_atom_map;
    size_t atom_site_row_count{0};
    size_t atom_site_skip_count{0};
    ParseLoopBlock(
        state,
        "_atom_site.",
        [&state, &altloc_primary_atom_map, &atom_site_row_count, &atom_site_skip_count](
            const CifColumnIndexMap& index_map,
            const std::vector<std::string>& token_list) {
            ++atom_site_row_count;
            const std::string context{
                "LoadAtomSiteBlock[row " + std::to_string(atom_site_row_count) + "]"};

            auto group_type{GetTokenOptional(index_map, token_list, {"group_PDB"}).value_or("ATOM")};
            auto element_type{GetTokenOptional(index_map, token_list, {"type_symbol"}).value_or("UNK")};
            auto comp_id{GetTokenOptional(index_map, token_list, {"label_comp_id", "auth_comp_id"})
                             .value_or("UNK")};
            auto atom_id_opt{GetTokenOptional(index_map, token_list, {"label_atom_id", "auth_atom_id"})};
            auto indicator{GetTokenOptional(index_map, token_list, {"label_alt_id"}).value_or(".")};
            auto sequence_id_token{
                GetTokenOptional(index_map, token_list, {"label_seq_id", "auth_seq_id"}).value_or(".")};
            auto serial_id{
                GetTokenOptional(index_map, token_list, {"id"}).value_or(std::to_string(atom_site_row_count))};
            auto chain_id{
                GetTokenOptional(index_map, token_list, {"label_asym_id", "auth_asym_id"}).value_or("?")};
            auto position_x_str{GetTokenOptional(index_map, token_list, {"Cartn_x"})};
            auto position_y_str{GetTokenOptional(index_map, token_list, {"Cartn_y"})};
            auto position_z_str{GetTokenOptional(index_map, token_list, {"Cartn_z"})};
            auto occupancy_str{GetTokenOptional(index_map, token_list, {"occupancy"})};
            auto temperature_str{GetTokenOptional(index_map, token_list, {"B_iso_or_equiv"})};
            auto model_number_str{GetTokenOptional(index_map, token_list, {"pdbx_PDB_model_num"})};

            if (!atom_id_opt.has_value()) {
                Logger::Log(LogLevel::Warning, context + " Missing atom id, skip this row.");
                ++atom_site_skip_count;
                return;
            }

            if (IsMmCifMissingValue(indicator))
                indicator = ".";
            auto parsed_x{TryParseFloat(position_x_str.value_or("?"))};
            auto parsed_y{TryParseFloat(position_y_str.value_or("?"))};
            auto parsed_z{TryParseFloat(position_z_str.value_or("?"))};
            if (!parsed_x.has_value() || !parsed_y.has_value() || !parsed_z.has_value()) {
                Logger::Log(
                    LogLevel::Warning,
                    context + " Invalid/missing Cartesian coordinates, skip this row.");
                ++atom_site_skip_count;
                return;
            }

            auto occupancy{
                ParseFloatOrDefault(occupancy_str.value_or("?"), 1.0f, "_atom_site.occupancy", context)};
            auto temperature{
                ParseFloatOrDefault(
                    temperature_str.value_or("?"), 0.0f, "_atom_site.B_iso_or_equiv", context)};
            auto model_number_id{
                ParseIntOrDefault(
                    model_number_str.value_or("1"),
                    1,
                    "_atom_site.pdbx_PDB_model_num",
                    context)};
            if (IsMmCifMissingValue(sequence_id_token))
                sequence_id_token = ".";
            auto sequence_id_value{
                ParseIntOrDefault(sequence_id_token, -1, "_atom_site.(label/auth)_seq_id", context)};
            auto serial_id_value{
                ParseIntOrDefault(
                    serial_id,
                    static_cast<int>(atom_site_row_count),
                    "_atom_site.id",
                    context)};

            auto atom_id{*atom_id_opt};
            StringHelper::EraseCharFromString(atom_id, '\"');

            auto element_enum{ChemicalDataHelper::GetElementFromString(element_type)};
            if (element_enum == Element::HYDROGEN)
                return;
            state.observed_component_atom_id_set_map[comp_id].insert(atom_id);

            auto is_special_atom{group_type == "HETATM"};

            if (!state.find_chemical_component_entry) {
                BuildDefaultChemicalComponentEntry(state, comp_id);
            }
            state.data_block->GetComponentKeySystemPtr()->RegisterComponent(comp_id);
            auto component_key{state.data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id)};
            if (!state.data_block->HasChemicalComponentEntry(component_key)) {
                BuildDefaultChemicalComponentEntry(state, comp_id);
                component_key = state.data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id);
            }

            if (!state.find_component_atom_entry) {
                BuildDefaultComponentAtomEntry(state, comp_id, atom_id, element_type);
            }
            state.data_block->GetAtomKeySystemPtr()->RegisterAtom(atom_id);
            auto atom_key{state.data_block->GetAtomKeySystemPtr()->GetAtomKey(atom_id)};

            AtomAltLocKey atom_altloc_key{
                model_number_id, chain_id, comp_id, sequence_id_token, atom_id};
            auto add_atom_as_primary =
                [&](std::unique_ptr<AtomObject> atom_object, const AtomAltLocKey& key) -> AtomObject* {
                auto raw_atom_ptr{atom_object.get()};
                state.data_block->AddAtomObject(
                    model_number_id, std::move(atom_object), sequence_id_token);
                altloc_primary_atom_map[key] = raw_atom_ptr;
                return raw_atom_ptr;
            };

            if (indicator != "." && indicator != "A") {
                auto primary_iter{altloc_primary_atom_map.find(atom_altloc_key)};
                if (primary_iter != altloc_primary_atom_map.end()) {
                    primary_iter->second->AddAlternatePosition(indicator, {*parsed_x, *parsed_y, *parsed_z});
                    primary_iter->second->AddAlternateOccupancy(indicator, occupancy);
                    primary_iter->second->AddAlternateTemperature(indicator, temperature);
                    return;
                }
            }

            auto atom_object{std::make_unique<AtomObject>()};
            atom_object->SetComponentID(comp_id);
            atom_object->SetComponentKey(component_key);
            atom_object->SetAtomID(atom_id);
            atom_object->SetAtomKey(atom_key);
            atom_object->SetElement(element_enum);
            atom_object->SetIndicator(indicator);
            atom_object->SetSequenceID(sequence_id_value);
            atom_object->SetSerialID(serial_id_value);
            atom_object->SetChainID(chain_id);
            atom_object->SetPosition(*parsed_x, *parsed_y, *parsed_z);
            atom_object->SetOccupancy(occupancy);
            atom_object->SetTemperature(temperature);
            atom_object->SetSpecialAtomFlag(is_special_atom);
            state.data_block->SetStructureInfo(atom_object.get());

            if (indicator == "." || indicator == "A") {
                add_atom_as_primary(std::move(atom_object), atom_altloc_key);
                return;
            }

            add_atom_as_primary(std::move(atom_object), atom_altloc_key);
        });
    Logger::Log(
        LogLevel::Info,
        "LoadAtomSiteBlock parsed rows = " + std::to_string(atom_site_row_count) + ", skipped rows = " + std::to_string(atom_site_skip_count) + ".");
}

struct CanonicalAtomPair {
    const AtomObject* atom_a;
    const AtomObject* atom_b;

    bool operator==(const CanonicalAtomPair& other) const {
        return atom_a == other.atom_a && atom_b == other.atom_b;
    }
};

struct CanonicalAtomPairHash {
    size_t operator()(const CanonicalAtomPair& pair) const {
        auto a{reinterpret_cast<std::uintptr_t>(pair.atom_a)};
        auto b{reinterpret_cast<std::uintptr_t>(pair.atom_b)};
        return std::hash<std::uintptr_t>{}(a ^ (b << 1));
    }
};

CanonicalAtomPair BuildCanonicalAtomPair(const AtomObject* atom_1, const AtomObject* atom_2) {
    if (atom_1 < atom_2)
        return {atom_1, atom_2};
    return {atom_2, atom_1};
}

void ConstructBondList(CifFormatState& state) {
    BuildPepetideBondEntry(state);
    BuildPhosphodiesterBondEntry(state);
    BuildObservedPolymerComponentBondEntry(state);
    const auto bond_count_before{state.data_block->GetBondObjectList().size()};
    const auto& atom_object_map{state.data_block->GetAtomObjectMap()};
    auto bond_key_system{state.data_block->GetBondKeySystemPtr()};
    auto component_key_system{state.data_block->GetComponentKeySystemPtr()};
    for (const auto& [model_number, atom_object_list] : atom_object_map) {
        if (atom_object_list.empty())
            continue;
        std::vector<AtomObject*> atom_ptr_list;
        atom_ptr_list.reserve(atom_object_list.size());
        for (const auto& atom : atom_object_list) {
            atom_ptr_list.emplace_back(atom.get());
        }
        auto kd_tree_root{KDTreeAlgorithm<AtomObject>::BuildKDTree(atom_ptr_list, 0)};
        std::unordered_set<CanonicalAtomPair, CanonicalAtomPairHash> processed_bond_pair_set;
        processed_bond_pair_set.reserve(atom_object_list.size() * 2);
        for (const auto& atom : atom_object_list) {
            auto component_id_1{atom->GetComponentID()};
            auto atom_id_1{atom->GetAtomID()};
            auto sequence_id_1{atom->GetSequenceID()};
            auto chain_id_1{atom->GetChainID()};
            auto neighbor_atom_list{
                KDTreeAlgorithm<AtomObject>::RangeSearch(
                    kd_tree_root.get(), atom.get(), kBondSearchingRadius)};

            for (auto neighbor_atom : neighbor_atom_list) {
                if (neighbor_atom == atom.get())
                    continue;
                auto canonical_pair{BuildCanonicalAtomPair(atom.get(), neighbor_atom)};
                if (processed_bond_pair_set.find(canonical_pair) != processed_bond_pair_set.end()) {
                    continue;
                }
                processed_bond_pair_set.insert(canonical_pair);

                auto component_id_2{neighbor_atom->GetComponentID()};
                auto atom_id_2{neighbor_atom->GetAtomID()};
                auto sequence_id_2{neighbor_atom->GetSequenceID()};
                auto chain_id_2{neighbor_atom->GetChainID()};
                if (!bond_key_system->IsRegistedBond(atom_id_1, atom_id_2))
                    continue;
                auto component_key_1{component_key_system->GetComponentKey(component_id_1)};
                auto bond_key{bond_key_system->GetBondKey(atom_id_1, atom_id_2)};
                if (!state.data_block->HasComponentBondEntry(component_key_1, bond_key))
                    continue;

                bool is_in_same_component{component_id_1 == component_id_2};
                bool is_in_same_chain{chain_id_1 == chain_id_2};
                bool is_in_consecutive_sequence{sequence_id_1 + 1 == sequence_id_2};

                bool is_peptide_bond{
                    !is_in_same_component &&
                    is_in_same_chain &&
                    is_in_consecutive_sequence &&
                    (bond_key == static_cast<BondKey>(Link::C_N))};

                bool is_phosphodiester_bond{
                    !is_in_same_component &&
                    is_in_same_chain &&
                    is_in_consecutive_sequence &&
                    (bond_key == static_cast<BondKey>(Link::P_O3p))};

                if (!is_in_same_component &&
                    !is_peptide_bond &&
                    !is_phosphodiester_bond) {
                    continue;
                }

                auto bond_entry{state.data_block->GetComponentBondEntryPtr(component_key_1, bond_key)};
                if (bond_entry == nullptr)
                    continue;

                auto bond_object{std::make_unique<BondObject>(atom.get(), neighbor_atom)};
                bond_object->SetBondKey(bond_key);
                bond_object->SetBondType(bond_entry->bond_type);
                bond_object->SetBondOrder(bond_entry->bond_order);
                bond_object->SetSpecialBondFlag(false);
                state.data_block->AddBondObject(std::move(bond_object));
            }
        }
    }
    Logger::Log(
        LogLevel::Info,
        "Construct " +
            std::to_string(state.data_block->GetBondObjectList().size() - bond_count_before) + " bonds (total " + std::to_string(state.data_block->GetBondObjectList().size()) + ").");
}

void WriteAtomSiteBlockEntry(
    const AtomObject* atom,
    const std::array<float, 3>& position,
    const std::string& alt_id,
    float occupancy,
    float temperature,
    int model_number,
    std::ostream& stream) {
    std::string group_pdb{atom->GetSpecialAtomFlag() ? "HETATM" : "ATOM"};
    std::string type_symbol{ChemicalDataHelper::GetLabel(atom->GetElement())};
    std::string label_atom_id{atom->GetAtomID()};
    std::string label_alt_id{alt_id.empty() ? "." : alt_id};
    std::string label_comp_id{atom->GetComponentID()};
    std::string label_asym_id{atom->GetChainID()};
    std::string label_seq_id{
        (atom->GetSequenceID() == -1) ? "." : std::to_string(atom->GetSequenceID())};

    stream << group_pdb << ' '
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

void WriteAtomSiteBlock(
    const ModelObject& model_object,
    std::ostream& stream,
    int model_par) {
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

    const int model_number{1};
    for (const auto& atom_ptr : model_object.GetAtomList()) {
        const AtomObject* atom{atom_ptr.get()};
        if (ModelAnalysisData::FindLocalEntry(*atom) == nullptr)
            continue;
        const auto & model_entry{ModelAnalysisData::RequireLocalEntry(*atom)};
        auto gaus_estimate{model_entry.GetEstimateMDPDE().GetParameter(model_par)};
        auto position{atom->GetPosition()};
        WriteAtomSiteBlockEntry(
            atom,
            position,
            atom->GetIndicator(),
            atom->GetOccupancy(),
            static_cast<float>(gaus_estimate),
            model_number,
            stream);
    }
    stream << "#\n";
}

} // namespace

CifFormat::CifFormat() : m_read_state{std::make_unique<CifFormatState>()} {
}

CifFormat::~CifFormat() = default;

std::unique_ptr<ModelObject> CifFormat::ReadModel(
    std::istream& stream,
    const std::string& source_name) {
    ParseMmCifDocument(stream, source_name, *m_read_state);
    LoadChemicalComponentBlock(*m_read_state);
    LoadDatabaseBlock(*m_read_state);
    LoadEntityBlock(*m_read_state);
    LoadPdbxData(*m_read_state);
    LoadAtomTypeBlock(*m_read_state);
    LoadStructureConformationBlock(*m_read_state);
    LoadStructureSheetBlock(*m_read_state);

    if (!m_read_state->find_component_bond_entry) {
        Logger::Log(
            LogLevel::Info,
            "No component bond entry found in the CIF file. Build default component bond entries.");
        BuildDefaultComponentBondEntry(*m_read_state);
    }
    LoadAtomSiteBlock(*m_read_state);
    ConstructBondList(*m_read_state);
    LoadStructureConnectionBlock(*m_read_state);
    LogHeaderSummary(*m_read_state->data_block);
    return m_read_state->data_block->TakeModelObject();
}

void CifFormat::WriteModel(
    const ModelObject& model_object,
    std::ostream& stream,
    int model_par) {
    stream << "data_" << model_object.GetKeyTag() << '\n';
    stream << "#\n";
    WriteAtomSiteBlock(model_object, stream, model_par);
}

} // namespace rhbm_gem
