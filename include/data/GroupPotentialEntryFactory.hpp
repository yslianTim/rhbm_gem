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
        if (class_key == AtomicInfoHelper::GetResidueClassKey())
        {
            using Tag = ResidueGroupClassifierTag;
            using GroupKey = typename GroupKeyMapping<Tag>::type;
            return std::make_unique<GroupPotentialEntry<GroupKey>>();
        }
        else if (class_key == AtomicInfoHelper::GetElementClassKey())
        {
            using Tag = ElementGroupClassifierTag;
            using GroupKey = typename GroupKeyMapping<Tag>::type;
            return std::make_unique<GroupPotentialEntry<GroupKey>>();
        }
        throw std::runtime_error("Unknown class key: " + class_key);
    }

    static std::any CreateGroupKeyTuple(const std::string & class_key, const AtomObject * atom)
    {
        if (class_key == AtomicInfoHelper::GetResidueClassKey())
        {
            return std::make_tuple(
                static_cast<int>(atom->GetResidue()),
                static_cast<int>(atom->GetElement()),
                static_cast<int>(atom->GetRemoteness()),
                static_cast<int>(atom->GetBranch()),
                atom->GetSpecialAtomFlag());
        }
        else if (class_key == AtomicInfoHelper::GetElementClassKey())
        {
            return std::make_tuple(
                static_cast<int>(atom->GetElement()),
                static_cast<int>(atom->GetRemoteness()),
                atom->GetSpecialAtomFlag());
        }
        throw std::runtime_error("Unknown class key: " + class_key);
    }
};