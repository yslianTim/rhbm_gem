#pragma once

#include <cstdint>
#include <tuple>
#include <type_traits>

#include "GlobalEnumClass.hpp"

static_assert(std::is_same_v<std::underlying_type_t<Residue>,    uint8_t>);
static_assert(std::is_same_v<std::underlying_type_t<Element>,    uint8_t>);
static_assert(std::is_same_v<std::underlying_type_t<Remoteness>, uint8_t>);
static_assert(std::is_same_v<std::underlying_type_t<Branch>,     uint8_t>);
static_assert(std::is_same_v<std::underlying_type_t<Structure>,  uint8_t>);

struct KeyPackerElementClass
{
    // Bits allocation
    // ┌─0~7─Element─┐┌─8~15─Remo.─┐┌─16─Flag─┐
    // | 8 bits      || 8 bits     || 1 bit   |
    static uint64_t Pack(Element element, Remoteness remoteness, bool flag)
    {
        return static_cast<uint64_t>(element)
            | (static_cast<uint64_t>(remoteness)   << 8)
            | (static_cast<uint64_t>(flag ? 1 : 0) << 16);
    }

    static std::tuple<Element, Remoteness, bool> Unpack(uint64_t key)
    {
        constexpr uint64_t mask_8bit { 0xFF };

        auto element{    static_cast<Element>(   (key      ) & mask_8bit) };
        auto remoteness{ static_cast<Remoteness>((key >>  8) & mask_8bit ) };
        bool flag{       static_cast<bool>(      (key >> 16) & 0x1) };

        return { element, remoteness, flag };
    }
};
struct KeyPackerResidueClass
{
    // Bits allocation
    // ┌─0~7─Residue─┐┌─8~15─Element─┐┌─16~23─Remo.─┐┌─24~31─Branch─┐┌─32─Flag─┐
    // | 8 bits      || 8 bits       || 8 bits      || 8 bits       || 1 bit   |
    static uint64_t Pack(
        Residue residue, Element element, Remoteness remoteness, Branch branch, bool flag)
    {
        return static_cast<uint64_t>(residue)
            | (static_cast<uint64_t>(element)      << 8)
            | (static_cast<uint64_t>(remoteness)   << 16)
            | (static_cast<uint64_t>(branch)       << 24)
            | (static_cast<uint64_t>(flag ? 1 : 0) << 32);
    }

    static std::tuple<Residue, Element, Remoteness, Branch, bool> Unpack(uint64_t key)
    {
        constexpr uint64_t mask_8bit { 0xFF };    //  8 bits mask

        auto residue{    static_cast<Residue>(   (key      ) & mask_8bit) };
        auto element{    static_cast<Element>(   (key >>  8) & mask_8bit) };
        auto remoteness{ static_cast<Remoteness>((key >> 16) & mask_8bit ) };
        auto branch{     static_cast<Branch>(    (key >> 24) & mask_8bit ) };
        bool flag{       static_cast<bool>(      (key >> 32) & 0x1 ) };

        return { residue, element, remoteness, branch, flag };
    }
};

struct KeyPackerStructureClass
{
    // Bits allocation
    // ┌─0~7─Struc.─┐┌─8~15─Residue─┐┌─16~23─Element─┐┌─24~31─Remo.─┐┌─32~39─Branch─┐┌─40─Flag─┐
    // | 8 bits     || 8 bits       || 8 bits        || 8 bits      || 8 bits       || 1 bit   |
    static uint64_t Pack(
        Structure structure, Residue residue,
        Element element, Remoteness remoteness, Branch branch, bool flag)
    {
        return static_cast<uint64_t>(structure)
            | (static_cast<uint64_t>(residue)      <<  8)
            | (static_cast<uint64_t>(element)      << 16)
            | (static_cast<uint64_t>(remoteness)   << 24)
            | (static_cast<uint64_t>(branch)       << 32)
            | (static_cast<uint64_t>(flag ? 1 : 0) << 40);
    }

    static std::tuple<Structure, Residue, Element, Remoteness, Branch, bool> Unpack(uint64_t key)
    {
        constexpr uint64_t mask_8bit { 0xFF };    //  8 bits mask

        auto structure{  static_cast<Structure>( (key      ) & mask_8bit) };
        auto residue{    static_cast<Residue>(   (key >>  8) & mask_8bit) };
        auto element{    static_cast<Element>(   (key >> 16) & mask_8bit) };
        auto remoteness{ static_cast<Remoteness>((key >> 24) & mask_8bit) };
        auto branch{     static_cast<Branch>(    (key >> 32) & mask_8bit) };
        bool flag{       static_cast<bool>(      (key >> 40) & 0x1 ) };

        return { structure, residue, element, remoteness, branch, flag };
    }
};
