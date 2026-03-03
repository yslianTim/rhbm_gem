#pragma once

namespace rhbm_gem {

class FileReaderBase
{
public:
    virtual ~FileReaderBase() = default;
    virtual void Read() = 0;
    virtual bool IsSuccessfullyRead() const = 0;
};

} // namespace rhbm_gem
