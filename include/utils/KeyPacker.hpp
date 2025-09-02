#pragma once

#include <cstdint>
#include <tuple>
#include <type_traits>

#include "GlobalEnumClass.hpp"
#include "ComponentKeySystem.hpp"
#include "AtomKeySystem.hpp"
#include "AtomicInfoHelper.hpp"

static_assert(std::is_same_v<std::underlying_type_t<Residue>,    uint8_t>);
static_assert(std::is_same_v<std::underlying_type_t<Element>,    uint8_t>);
static_assert(std::is_same_v<std::underlying_type_t<Remoteness>, uint8_t>);
static_assert(std::is_same_v<std::underlying_type_t<Branch>,     uint8_t>);
static_assert(std::is_same_v<std::underlying_type_t<Structure>,  uint8_t>);
static_assert(std::is_same_v<ComponentKey, uint8_t>);
static_assert(std::is_same_v<AtomKey, uint32_t>);

struct KeyPackerElementClass
{
    // Bits allocation (Type 1)
    // ┌─0~7─Element─┐┌─8~15─Remo.─┐┌─16─Flag─┐
    // | 8 bits      || 8 bits     || 1 bit   |
    static uint64_t Pack(Element element, Remoteness remoteness, bool flag)
    {
        return static_cast<uint64_t>(element)
            | (static_cast<uint64_t>(remoteness)   <<  8)
            | (static_cast<uint64_t>(flag ? 1 : 0) << 16);
    }

    // Bits allocation (Type 2)
    // ┌─0~31─AtomKey─┐┌─32─Flag─┐
    // | 32 bits      || 1 bit   |
    static uint64_t Pack(AtomKey atom_key, bool flag)
    {
        return static_cast<uint64_t>(atom_key)
            | (static_cast<uint64_t>(flag ? 1 : 0) << 32);
    }

    static std::tuple<Element, Remoteness, bool> Unpack(uint64_t key)
    {
        constexpr uint64_t mask_8bit { 0xFF };
        constexpr uint64_t mask_32bit{ 0xFFFFFFFF };

        // Type 2 detection: keys produced by the AtomKey-based Pack use bits
        // beyond bit 16 or the dedicated flag bit at position 32. Keys created
        // by the Element/Remoteness pack only occupy bits 0-16.
        const bool is_type2{
            (key >> 17) != 0 || (key & (static_cast<uint64_t>(1) << 32)) != 0
        };

        if (is_type2)
        {
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

        // Type 1 key
        auto element{    static_cast<Element>(   (key      ) & mask_8bit) };
        auto remoteness{ static_cast<Remoteness>((key >>  8) & mask_8bit) };
        bool flag{       static_cast<bool>(      (key >> 16) & 0x1) };
        return { element, remoteness, flag };
    }
};
struct KeyPackerResidueClass
{
    // Bits allocation (Type 1)
    // ┌─0~7─Residue─┐┌─8~15─Element─┐┌─16~23─Remo.─┐┌─24~31─Branch─┐┌─32─Flag─┐
    // | 8 bits      || 8 bits       || 8 bits      || 8 bits       || 1 bit   |
    static uint64_t Pack(
        Residue residue, Element element, Remoteness remoteness, Branch branch, bool flag)
    {
        return static_cast<uint64_t>(residue)
            | (static_cast<uint64_t>(element)      <<  8)
            | (static_cast<uint64_t>(remoteness)   << 16)
            | (static_cast<uint64_t>(branch)       << 24)
            | (static_cast<uint64_t>(flag ? 1 : 0) << 32);
    }

    // Bits allocation (Type 2)
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
        constexpr uint64_t mask_8bit { 0xFF };    //  8 bits mask

        // Type 2 key has data in bits >= 40
        if ((key >> 40) != 0)
        {
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

        // Type 1 key
        auto residue{    static_cast<Residue>(   (key      ) & mask_8bit) };
        auto element{    static_cast<Element>(   (key >>  8) & mask_8bit) };
        auto remoteness{ static_cast<Remoteness>((key >> 16) & mask_8bit) };
        auto branch{     static_cast<Branch>(    (key >> 24) & mask_8bit) };
        bool flag{       static_cast<bool>(      (key >> 32) & 0x1 ) };
        return { residue, element, remoteness, branch, flag };
    }
};

struct KeyPackerStructureClass
{
    // Bits allocation (Type 1)
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

    // Bits allocation (Type 2)
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

        // Detect type-2 key when any bits beyond 47 are set (including the flag
        // at bit 48). Type-1 keys only use bits up to 40.
        if (key >> 48)
        {
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

        // Type 1 key
        auto structure{  static_cast<Structure>( (key      ) & mask_8bit) };
        auto residue{    static_cast<Residue>(   (key >>  8) & mask_8bit) };
        auto element{    static_cast<Element>(   (key >> 16) & mask_8bit) };
        auto remoteness{ static_cast<Remoteness>((key >> 24) & mask_8bit) };
        auto branch{     static_cast<Branch>(    (key >> 32) & mask_8bit) };
        bool flag{       static_cast<bool>(      (key >> 40) & 0x1 ) };
        return { structure, residue, element, remoteness, branch, flag };
    }
};
