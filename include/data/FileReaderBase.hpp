#pragma once

class FileReaderBase
{
public:
    virtual ~FileReaderBase() = default;
    virtual void Read(void) = 0;
};