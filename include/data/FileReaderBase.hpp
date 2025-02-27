#pragma once

#include <iostream>
#include <memory>
#include <string>

class FileReaderBase
{
public:
    virtual ~FileReaderBase() = default;
    virtual void Read(void) = 0;
};