#pragma once

#include <cstdint>
#include <tuple>
#include <type_traits>

#include "GlobalEnumClass.hpp"
#include "ComponentKeySystem.hpp"
#include "AtomKeySystem.hpp"
#include "AtomicInfoHelper.hpp"

static_assert(std::is_same_v<std::underlying_type_t<Residue>,    uint8_t>);
static_assert(std::is_same_v<std::underlying_type_t<Spot>,       uint32_t>);
static_assert(std::is_same_v<std::underlying_type_t<Element>,    uint8_t>);
static_assert(std::is_same_v<std::underlying_type_t<Remoteness>, uint8_t>);
static_assert(std::is_same_v<std::underlying_type_t<Branch>,     uint8_t>);
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

    /*
    static std::tuple<AtomKey, bool> Unpack(uint64_t key)
    {
        constexpr uint64_t mask_32bit{ 0xFFFFFFFF };
        AtomKey atom_key{ static_cast<AtomKey>(key & mask_32bit) };
        bool flag{ static_cast<bool>((key >> 32) & 0x1) };
        return { atom_key, flag };
    }*/

    static std::tuple<Element, Remoteness, bool> Unpack(uint64_t key)
    {
        constexpr uint64_t mask_8bit { 0xFF };
        constexpr uint64_t mask_32bit{ 0xFFFFFFFF };

        AtomKey atom_key{ static_cast<AtomKey>(key & mask_32bit) };
        bool flag{ static_cast<bool>((key >> 32) & 0x1) };
        Element element;
        Remoteness remoteness;
        if (AtomKeySystem::Instance().IsBuildInAtom(atom_key))
        {
            element    = static_cast<Element>(   atom_key        & mask_8bit);
            remoteness = static_cast<Remoteness>((atom_key >> 8) & mask_8bit);
        }
        else
        {
            element    = Element::UNK;
            remoteness = Remoteness::UNK;
        }
        return { element, remoteness, flag };
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

    static std::tuple<Residue, Element, Remoteness, Branch, bool> Unpack(uint64_t key)
    {
        constexpr uint64_t mask_8bit  { 0xFF };    //  8 bits mask
        constexpr uint64_t mask_32bit { 0xFFFFFFFF };

        auto component_key{ static_cast<ComponentKey>( (key      ) & mask_8bit) };
        auto atom_key{      static_cast<AtomKey>(      (key >> 8) & mask_32bit) };
        bool flag{          static_cast<bool>(         (key >> 40) & 0x1) };

        Residue residue{ Residue::UNK };
        auto comp_id{ ComponentKeySystem::Instance().GetComponentId(component_key) };
        residue = AtomicInfoHelper::GetResidueFromString(comp_id);

        Element element{ Element::UNK };
        Remoteness remoteness{ Remoteness::UNK };
        Branch branch{ Branch::UNK };

        if (AtomKeySystem::Instance().IsBuildInAtom(atom_key))
        {
            element    = static_cast<Element>(   (atom_key      ) & mask_8bit);
            remoteness = static_cast<Remoteness>((atom_key >>  8) & mask_8bit);
            branch     = static_cast<Branch>(    (atom_key >> 16) & mask_8bit);
        }
        return { residue, element, remoteness, branch, flag };
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

    static std::tuple<Structure, Residue, Element, Remoteness, Branch, bool> Unpack(uint64_t key)
    {
        constexpr uint64_t mask_8bit { 0xFF };    //  8 bits mask

        auto structure{ static_cast<Structure>( (key      ) & mask_8bit) };
        auto component_key{ static_cast<ComponentKey>((key >> 8)  & mask_8bit) };
        auto atom_key{ static_cast<AtomKey>( (key >> 16) & 0xFFFFFFFF ) };
        bool flag{ static_cast<bool>( (key >> 48) & 0x1) };

        // Resolve residue using ComponentKeySystem. Unknown keys default to UNK.
        Residue residue{ Residue::UNK };
        auto comp_id{ ComponentKeySystem::Instance().GetComponentId(component_key) };
        residue = AtomicInfoHelper::GetResidueFromString(comp_id);

        // Decode atom key to element/remoteness/branch for build-in atoms.
        Element element{ Element::UNK };
        Remoteness remoteness{ Remoteness::UNK };
        Branch branch{ Branch::UNK };
        if (AtomKeySystem::Instance().IsBuildInAtom(atom_key))
        {
            element    = static_cast<Element>(   (atom_key      ) & mask_8bit);
            remoteness = static_cast<Remoteness>((atom_key >>  8) & mask_8bit);
            branch     = static_cast<Branch>(    (atom_key >> 16) & mask_8bit);
        }
        return { structure, residue, element, remoteness, branch, flag };
    }
};
