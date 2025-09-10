#pragma once

#include <cstdint>
#include <tuple>
#include <string>

#include "GlobalEnumClass.hpp"
#include "ComponentKeySystem.hpp"
#include "AtomKeySystem.hpp"

static_assert(std::is_same_v<std::underlying_type_t<Structure>,  uint8_t>);
static_assert(std::is_same_v<GroupKey, uint64_t>);
static_assert(std::is_same_v<ComponentKey, uint8_t>);
static_assert(std::is_same_v<AtomKey, uint16_t>);

struct KeyPackerSimpleAtomClass
{
    // Bits allocation
    // ┌─0~15─AtomKey─┐┌─16─Flag─┐
    // | 16 bits      || 1 bit   |
    static GroupKey Pack(AtomKey atom_key, bool flag)
    {
        return static_cast<GroupKey>(atom_key)
            | (static_cast<GroupKey>(flag ? 1 : 0) << 16);
    }

    static std::tuple<AtomKey, bool> Unpack(GroupKey key)
    {
        constexpr GroupKey mask_16bit{ 0xFFFF };
        AtomKey atom_key{ static_cast<AtomKey>((key      ) & mask_16bit) };
        bool flag{        static_cast<bool>(   (key >> 16) & 0x1)        };
        return { atom_key, flag };
    }

    static std::string GetKeyString(GroupKey key)
    {
        auto [atom_key, flag]{ Unpack(key) };
        return "<" + std::to_string(static_cast<int>(atom_key)) + ", "
                   + std::to_string(flag) + ">";
    }
};

struct KeyPackerSimpleBondClass
{
    // Bits allocation
    // ┌─0~15─AtomKey1─┐┌─16~31─AtomKey2─┐
    // | 16 bits       || 16 bits.       |
    static GroupKey Pack(AtomKey atom_key_1, AtomKey atom_key_2)
    {
        return static_cast<GroupKey>(atom_key_1)
            | (static_cast<GroupKey>(atom_key_2) << 16);
    }

    static std::tuple<AtomKey, AtomKey> Unpack(GroupKey key)
    {
        constexpr GroupKey mask_16bit{ 0xFFFF };
        AtomKey atom_key_1{ static_cast<AtomKey>((key      ) & mask_16bit) };
        AtomKey atom_key_2{ static_cast<AtomKey>((key >> 16) & mask_16bit) };
        return { atom_key_1, atom_key_2 };
    }

    static std::string GetKeyString(GroupKey key)
    {
        auto [atom_key_1, atom_key_2]{ Unpack(key) };
        return "<" + std::to_string(static_cast<int>(atom_key_1)) + ", "
                   + std::to_string(static_cast<int>(atom_key_2)) + ">";
    }
};

struct KeyPackerComponentAtomClass
{
    // Bits allocation
    // ┌─0~7─ComponentKey─┐┌─8~23─AtomKey─┐
    // | 8 bits           || 16 bits      |
    static GroupKey Pack(ComponentKey component_key, AtomKey atom_key)
    {
        return static_cast<GroupKey>(component_key)
            | (static_cast<GroupKey>(atom_key)      <<  8);
    }

    static std::tuple<ComponentKey, AtomKey> Unpack(GroupKey key)
    {
        constexpr GroupKey mask_8bit { 0xFF };
        constexpr GroupKey mask_16bit{ 0xFFFF };
        auto component_key{ static_cast<ComponentKey>( (key      ) & mask_8bit)  };
        auto atom_key{      static_cast<AtomKey>(      (key >> 8 ) & mask_16bit) };
        return { component_key, atom_key };
    }

    static std::string GetKeyString(GroupKey key)
    {
        auto [component_key, atom_key]{ Unpack(key) };
        return "<" + std::to_string(static_cast<int>(component_key)) + ", "
                   + std::to_string(static_cast<int>(atom_key)) + ">";
    }
};

struct KeyPackerComponentBondClass
{
    // Bits allocation
    // ┌─0~7─ComponentKey1─┐┌─8~15─ComponentKey2─┐┌─16~31─AtomKey1─┐┌─32~47─AtomKey2─┐
    // | 8 bits            || 8 bits             || 16 bits        || 16 bits        |
    static GroupKey Pack(
        ComponentKey component_key_1,
        ComponentKey component_key_2,
        AtomKey atom_key_1,
        AtomKey atom_key_2)
    {
        return static_cast<GroupKey>(component_key_1)
            | (static_cast<GroupKey>(component_key_2) <<  8)
            | (static_cast<GroupKey>(atom_key_1)      << 16)
            | (static_cast<GroupKey>(atom_key_2)      << 32);
    }

    static std::tuple<ComponentKey, ComponentKey, AtomKey, AtomKey> Unpack(GroupKey key)
    {
        constexpr GroupKey mask_8bit { 0xFF };
        constexpr GroupKey mask_16bit{ 0xFFFF };
        auto component_key_1{ static_cast<ComponentKey>( (key      ) & mask_8bit)  };
        auto component_key_2{ static_cast<ComponentKey>( (key >>  8) & mask_8bit)  };
        auto atom_key_1{      static_cast<AtomKey>(      (key >> 16) & mask_16bit) };
        auto atom_key_2{      static_cast<AtomKey>(      (key >> 32) & mask_16bit) };
        return { component_key_1, component_key_2, atom_key_1, atom_key_2 };
    }

    static std::string GetKeyString(GroupKey key)
    {
        auto [component_key_1, component_key_2, atom_key_1, atom_key_2]{ Unpack(key) };
        return "<" + std::to_string(static_cast<int>(component_key_1)) + ", "
                   + std::to_string(static_cast<int>(component_key_2)) + ", "
                   + std::to_string(static_cast<int>(atom_key_1)) + ", "
                   + std::to_string(static_cast<int>(atom_key_2)) + ">";
    }
};

struct KeyPackerStructureAtomClass
{
    // Bits allocation
    // ┌─0~7─Struc.─┐┌─8~15─ComponentKey─┐┌─16~31─AtomKey─┐
    // | 8 bits     || 8 bits            || 16 bits       |
    static GroupKey Pack(Structure structure, ComponentKey component_key, AtomKey atom_key)
    {
        return static_cast<GroupKey>(structure)
            | (static_cast<GroupKey>(component_key) <<  8)
            | (static_cast<GroupKey>(atom_key)      << 16);
    }

    static std::tuple<Structure, ComponentKey, AtomKey> Unpack(GroupKey key)
    {
        constexpr GroupKey mask_8bit { 0xFF };    //  8 bits mask
        constexpr GroupKey mask_16bit{ 0xFFFF };
        auto structure{ static_cast<Structure>(       (key      ) & mask_8bit)  };
        auto component_key{ static_cast<ComponentKey>((key >> 8)  & mask_8bit)  };
        auto atom_key{ static_cast<AtomKey>(          (key >> 16) & mask_16bit) };
        return { structure, component_key, atom_key };
    }

    static std::string GetKeyString(GroupKey key)
    {
        auto [structure, component_key, atom_key]{ Unpack(key) };
        return "<" + std::to_string(static_cast<int>(structure)) + ", "
                   + std::to_string(static_cast<int>(component_key)) + ", "
                   + std::to_string(static_cast<int>(atom_key)) + ">";
    }
};
