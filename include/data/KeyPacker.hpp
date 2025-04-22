#pragma once

#include <cstdint>
#include <tuple>
#include <type_traits>

enum class Residue : uint16_t;
enum class Element : uint16_t;
enum class Remoteness : uint8_t;
enum class Branch : uint8_t;

static_assert(std::is_same_v<std::underlying_type_t<Residue>, uint16_t>);
static_assert(std::is_same_v<std::underlying_type_t<Element>, uint16_t>);
static_assert(std::is_same_v<std::underlying_type_t<Remoteness>, uint8_t>);
static_assert(std::is_same_v<std::underlying_type_t<Branch>,    uint8_t>);

class KeyPackerResidueClass
{
public:
    // Bits allocation
    // ┌─0…15───Residue───┐┌─16…31──Element───┐┌─32…39─Remoteness─┐┌─40…47─Branch─┐┌─48─Flag─┐
    // | 16 bits          || 16 bits          || 8 bits           || 8 bits       || 1 bit   |
    static uint64_t Pack(
        Residue residue, Element element, Remoteness remoteness, Branch branch, bool flag)
    {
        return static_cast<uint64_t>(residue)
            | (static_cast<uint64_t>(element)      << 16)
            | (static_cast<uint64_t>(remoteness)   << 32)
            | (static_cast<uint64_t>(branch)       << 40)
            | (static_cast<uint64_t>(flag ? 1 : 0) << 48);
    }

    static std::tuple<Residue, Element, Remoteness, Branch, bool> Unpack(uint64_t key)
    {
        constexpr uint64_t mask_16bit{ 0xFFFF };  // 16 bits mask
        constexpr uint64_t mask_8bit { 0xFF };    //  8 bits mask

        auto residue{    static_cast<Residue>(   (key      ) & mask_16bit) };
        auto element{    static_cast<Element>(   (key >> 16) & mask_16bit) };
        auto remoteness{ static_cast<Remoteness>((key >> 32) & mask_8bit ) };
        auto branch{     static_cast<Branch>(    (key >> 40) & mask_8bit ) };
        bool flag{ ((key >> 48    ) & 0x1 ) };

        return { residue, element, remoteness, branch, flag };
    }
};