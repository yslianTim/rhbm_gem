#include <iostream>

#include "BuiltInCommandCatalog.hpp"

int main()
{
    for (const auto & descriptor : rhbm_gem::BuiltInCommandCatalog())
    {
        if (!rhbm_gem::IsPythonPublic(descriptor.binding_exposure)) continue;
        if (descriptor.python_binding_name.empty()) continue;
        std::cout << descriptor.python_binding_name << '\n';
    }
    return 0;
}
