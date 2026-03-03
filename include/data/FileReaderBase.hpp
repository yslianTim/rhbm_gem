#pragma once

namespace rhbm_gem {

class FileReaderBase
{
public:
    virtual ~FileReaderBase() = default;
    virtual void Read(void) = 0;
    virtual bool IsSuccessfullyRead(void) const = 0;
};

} // namespace rhbm_gem
