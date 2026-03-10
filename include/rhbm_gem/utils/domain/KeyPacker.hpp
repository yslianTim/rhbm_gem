#pragma once

#include <cstdint>
#include <tuple>
#include <string>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/domain/ComponentKeySystem.hpp>
#include <rhbm_gem/utils/domain/AtomKeySystem.hpp>
#include <rhbm_gem/utils/domain/BondKeySystem.hpp>

static_assert(std::is_same_v<std::underlying_type_t<Structure>,  uint8_t>);
static_assert(std::is_same_v<GroupKey, uint64_t>);
static_assert(std::is_same_v<ComponentKey, uint8_t>);
static_assert(std::is_same_v<AtomKey, uint16_t>);
static_assert(std::is_same_v<BondKey, uint16_t>);

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
    // ┌─0~15─BondKey─┐┌─16─Flag─┐
    // | 16 bits      || 1 bit   |
    static GroupKey Pack(BondKey bond_key, bool flag)
    {
        return static_cast<GroupKey>(bond_key)
            | (static_cast<GroupKey>(flag ? 1 : 0) << 16);
    }

    static std::tuple<BondKey, bool> Unpack(GroupKey key)
    {
        constexpr GroupKey mask_16bit{ 0xFFFF };
        BondKey bond_key{ static_cast<BondKey>((key      ) & mask_16bit) };
        bool flag{        static_cast<bool>(   (key >> 16) & 0x1)        };
        return { bond_key, flag };
    }

    static std::string GetKeyString(GroupKey key)
    {
        auto [bond_key, flag]{ Unpack(key) };
        return "<" + std::to_string(static_cast<int>(bond_key)) + ", "
                   + std::to_string(flag) + ">";
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
        auto atom_key{      static_cast<AtomKey>(      (key >>  8) & mask_16bit) };
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
    // ┌─0~7─ComponentKey1─┐┌─8~23─BondKey─┐
    // | 8 bits            || 16 bits      |
    static GroupKey Pack(
        ComponentKey component_key,
        BondKey bond_key)
    {
        return static_cast<GroupKey>(component_key)
            | (static_cast<GroupKey>(bond_key)      <<  8);
    }

    static std::tuple<ComponentKey, BondKey> Unpack(GroupKey key)
    {
        constexpr GroupKey mask_8bit { 0xFF };
        constexpr GroupKey mask_16bit{ 0xFFFF };
        auto component_key{ static_cast<ComponentKey>( (key      ) & mask_8bit)  };
        auto bond_key{        static_cast<BondKey>(    (key >>  8) & mask_16bit) };
        return { component_key, bond_key };
    }

    static std::string GetKeyString(GroupKey key)
    {
        auto [component_key, bond_key]{ Unpack(key) };
        return "<" + std::to_string(static_cast<int>(component_key)) + ", "
                   + std::to_string(static_cast<int>(bond_key)) + ">";
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
