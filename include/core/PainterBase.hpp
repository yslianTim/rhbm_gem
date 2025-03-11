#pragma once

#include <iostream>
#include <string>

class PainterBase
{
public:
    virtual ~PainterBase() = default;
    virtual void SetFolder(std::string folder_path) = 0;
    virtual void Painting(void) = 0;

};