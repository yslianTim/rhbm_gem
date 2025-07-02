#pragma once

class FileReaderBase
{
public:
    virtual ~FileReaderBase() = default;
    virtual void Read(void) = 0;
    virtual bool IsSuccessfullyRead(void) const = 0;
};
