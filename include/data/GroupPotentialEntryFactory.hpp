#pragma once

#include <iostream>
#include <memory>

#include "GroupPotentialEntryBase.hpp"
#include "GroupPotentialEntry.hpp"
#include "AtomicInfoHelper.hpp"

class GroupPotentialEntryFactory
{
public:
    static std::unique_ptr<GroupPotentialEntryBase> Create(const std::string & classification)
    {
        if (classification == "residue_group")
        {
            using Tag = ResidueGroupClassifierTag;
            using GroupKey = typename GroupKeyMapping<Tag>::type;
            return std::make_unique<GroupPotentialEntry<GroupKey>>();
        }
        else if (classification == "element_group")
        {
            using Tag = ElementGroupClassifierTag;
            using GroupKey = typename GroupKeyMapping<Tag>::type;
            return std::make_unique<GroupPotentialEntry<GroupKey>>();
        }
        throw std::runtime_error("Unknown classification: " + classification);
    }
};