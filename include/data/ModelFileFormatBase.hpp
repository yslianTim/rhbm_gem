#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <any>

class AtomObject;

class ModelFileFormatBase
{
public:
    virtual ~ModelFileFormatBase() = default;
    virtual void LoadHeader(const std::string & filename) = 0;
    virtual void PrintHeader(void) const = 0;
    virtual void LoadDataArray(const std::string & filename) = 0;
    virtual void BuildAtomObject(std::any atom_info, bool is_special_atom) = 0;
    virtual std::vector<std::unique_ptr<AtomObject>> GetAtomObjectList(void) = 0;
    virtual std::string GetPdbID(void) const = 0;
    virtual std::string GetEmdID(void) const = 0;
    virtual double GetResolution(void) const = 0;
    virtual std::string GetResolutionMethod(void) const = 0;

};