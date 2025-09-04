#pragma once

#include <cstdint>
#include <tuple>

#include "GlobalEnumClass.hpp"
#include "ComponentKeySystem.hpp"
#include "AtomKeySystem.hpp"

static_assert(std::is_same_v<std::underlying_type_t<Structure>,  uint8_t>);
static_assert(std::is_same_v<ComponentKey, uint8_t>);
static_assert(std::is_same_v<AtomKey, uint32_t>);

struct KeyPackerElementClass
{
    // Bits allocation
    // ┌─0~31─AtomKey─┐┌─32─Flag─┐
    // | 32 bits      || 1 bit   |
    static uint64_t Pack(AtomKey atom_key, bool flag)
    {
        return static_cast<uint64_t>(atom_key)
            | (static_cast<uint64_t>(flag ? 1 : 0) << 32);
    }

    static std::tuple<AtomKey, bool> Unpack(uint64_t key)
    {
        constexpr uint64_t mask_32bit{ 0xFFFFFFFF };
        AtomKey atom_key{ static_cast<AtomKey>((key      ) & mask_32bit) };
        bool flag{        static_cast<bool>(   (key >> 32) & 0x1)        };
        return { atom_key, flag };
    }
};
struct KeyPackerResidueClass
{
    // Bits allocation
    // ┌─0~7─ComponentKey─┐┌─8~39─AtomKey─┐┌─40─Flag─┐
    // | 8 bits           || 32 bits      || 1 bit   |
    static uint64_t Pack(
        ComponentKey component_key, AtomKey atom_key, bool flag)
    {
        return static_cast<uint64_t>(component_key)
            | (static_cast<uint64_t>(atom_key)      <<  8)
            | (static_cast<uint64_t>(flag ? 1 : 0)  << 40);
    }

    static std::tuple<ComponentKey, AtomKey, bool> Unpack(uint64_t key)
    {
        constexpr uint64_t mask_8bit { 0xFF };    //  8 bits mask
        constexpr uint64_t mask_32bit{ 0xFFFFFFFF };
        auto component_key{ static_cast<ComponentKey>( (key      ) & mask_8bit)  };
        auto atom_key{      static_cast<AtomKey>(      (key >> 8 ) & mask_32bit) };
        bool flag{          static_cast<bool>(         (key >> 40) & 0x1)        };
        return { component_key, atom_key, flag };
    }
};

struct KeyPackerStructureClass
{
    // Bits allocation
    // ┌─0~7─Struc.─┐┌─8~15─ComponentKey─┐┌─16~47─AtomKey─┐┌─48─Flag─┐
    // | 8 bits     || 8 bits            || 32 bits       || 1 bit   |
    static uint64_t Pack(
        Structure structure, ComponentKey component_key, AtomKey atom_key, bool flag)
    {
        return static_cast<uint64_t>(structure)
            | (static_cast<uint64_t>(component_key) <<  8)
            | (static_cast<uint64_t>(atom_key)      << 16)
            | (static_cast<uint64_t>(flag ? 1 : 0)  << 48);
    }

    static std::tuple<Structure, ComponentKey, AtomKey, bool> Unpack(uint64_t key)
    {
        constexpr uint64_t mask_8bit { 0xFF };    //  8 bits mask
        constexpr uint64_t mask_32bit{ 0xFFFFFFFF };
        auto structure{ static_cast<Structure>(       (key      ) & mask_8bit)  };
        auto component_key{ static_cast<ComponentKey>((key >> 8)  & mask_8bit)  };
        auto atom_key{ static_cast<AtomKey>(          (key >> 16) & mask_32bit) };
        bool flag{ static_cast<bool>(                 (key >> 48) & 0x1)        };
        return { structure, component_key, atom_key, flag };
    }
};
