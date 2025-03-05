#pragma once

#include <iostream>
#include <memory>
#include <tuple>
#include <any>

#include "GroupPotentialEntryBase.hpp"
#include "GroupPotentialEntry.hpp"
#include "AtomicInfoHelper.hpp"
#include "AtomObject.hpp"

class GroupPotentialEntryFactory
{
public:
    static std::unique_ptr<GroupPotentialEntryBase> Create(const std::string & class_key)
    {
        if (class_key == "residue_class")
        {
            using Tag = ResidueGroupClassifierTag;
            using GroupKey = typename GroupKeyMapping<Tag>::type;
            return std::make_unique<GroupPotentialEntry<GroupKey>>();
        }
        else if (class_key == "element_class")
        {
            using Tag = ElementGroupClassifierTag;
            using GroupKey = typename GroupKeyMapping<Tag>::type;
            return std::make_unique<GroupPotentialEntry<GroupKey>>();
        }
        throw std::runtime_error("Unknown class key: " + class_key);
    }

    static std::any CreateGroupKeyTuple(const std::string & class_key, const AtomObject * atom)
    {
        if (class_key == "residue_class")
        {
            return std::make_tuple(
                atom->GetResidue(),
                atom->GetElement(),
                atom->GetRemoteness(),
                atom->GetBranch(),
                atom->GetSpecialAtomFlag());
        }
        else if (class_key == "element_class")
        {
            return std::make_tuple(
                atom->GetElement(),
                atom->GetRemoteness(),
                atom->GetSpecialAtomFlag());
        }
        throw std::runtime_error("Unknown class key: " + class_key);
    }
};