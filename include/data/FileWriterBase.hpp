#pragma once

class FileWriterBase
{
public:
    virtual ~FileWriterBase() = default;
    virtual void Write(void) = 0;
};
