#pragma once

namespace rhbm_gem {

class FileWriterBase
{
public:
    virtual ~FileWriterBase() = default;
    virtual void Write() = 0;
};

} // namespace rhbm_gem
