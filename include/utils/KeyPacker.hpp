#pragma once

#include <cstdint>
#include <tuple>
#include <type_traits>

enum class Residue : uint16_t;
enum class Element : uint16_t;
enum class Remoteness : uint8_t;
enum class Branch : uint8_t;
enum class Structure : uint8_t;

static_assert(std::is_same_v<std::underlying_type_t<Residue>,    uint16_t>);
static_assert(std::is_same_v<std::underlying_type_t<Element>,    uint16_t>);
static_assert(std::is_same_v<std::underlying_type_t<Remoteness>, uint8_t>);
static_assert(std::is_same_v<std::underlying_type_t<Branch>,     uint8_t>);
static_assert(std::is_same_v<std::underlying_type_t<Structure>,  uint8_t>);

struct KeyPackerElementClass
{
    // Bits allocation
    // ┌─0~15─Element─┐┌─16~23─Remo.─┐┌─24─Flag─┐
    // | 16 bits      || 8 bits      || 1 bit   |
    static uint64_t Pack(Element element, Remoteness remoteness, bool flag)
    {
        return static_cast<uint64_t>(element)
            | (static_cast<uint64_t>(remoteness)   << 16)
            | (static_cast<uint64_t>(flag ? 1 : 0) << 24);
    }

    static std::tuple<Element, Remoteness, bool> Unpack(uint64_t key)
    {
        constexpr uint64_t mask_16bit{ 0xFFFF };
        constexpr uint64_t mask_8bit { 0xFF };

        auto element{    static_cast<Element>(   (key      ) & mask_16bit) };
        auto remoteness{ static_cast<Remoteness>((key >> 16) & mask_8bit ) };
        bool flag{       static_cast<bool>(      (key >> 24) & 0x1) };

        return { element, remoteness, flag };
    }
};
struct KeyPackerResidueClass
{
    // Bits allocation
    // ┌─0~15─Residue─┐┌─16~31─Element─┐┌─32~39─Remo.─┐┌─40~47─Branch─┐┌─48─Flag─┐
    // | 16 bits      || 16 bits       || 8 bits      || 8 bits       || 1 bit   |
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
        bool flag{       static_cast<bool>(      (key >> 48) & 0x1 ) };

        return { residue, element, remoteness, branch, flag };
    }
};

struct KeyPackerStructureClass
{
    // Bits allocation
    // ┌─0~7─Struc.─┐┌─8~23─Residue─┐┌─24~39─Element─┐┌─40~47─Remo.─┐┌─48~55─Branch─┐┌─56─Flag─┐
    // | 8 bits     || 16 bits      || 16 bits       || 8 bits      || 8 bits       || 1 bit   |
    static uint64_t Pack(
        Structure structure, Residue residue,
        Element element, Remoteness remoteness, Branch branch, bool flag)
    {
        return static_cast<uint64_t>(structure)
            | (static_cast<uint64_t>(residue)      <<  8)
            | (static_cast<uint64_t>(element)      << 24)
            | (static_cast<uint64_t>(remoteness)   << 40)
            | (static_cast<uint64_t>(branch)       << 48)
            | (static_cast<uint64_t>(flag ? 1 : 0) << 56);
    }

    static std::tuple<Structure, Residue, Element, Remoteness, Branch, bool> Unpack(uint64_t key)
    {
        constexpr uint64_t mask_16bit{ 0xFFFF };  // 16 bits mask
        constexpr uint64_t mask_8bit { 0xFF };    //  8 bits mask

        auto structure{  static_cast<Structure>( (key      ) & mask_8bit ) };
        auto residue{    static_cast<Residue>(   (key >>  8) & mask_16bit) };
        auto element{    static_cast<Element>(   (key >> 24) & mask_16bit) };
        auto remoteness{ static_cast<Remoteness>((key >> 40) & mask_8bit ) };
        auto branch{     static_cast<Branch>(    (key >> 48) & mask_8bit ) };
        bool flag{       static_cast<bool>(      (key >> 56) & 0x1 ) };

        return { structure, residue, element, remoteness, branch, flag };
    }
};
