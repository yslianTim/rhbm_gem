#include <iostream>

#include "BuiltInCommandCatalogInternal.hpp"

int main()
{
    for (const auto & descriptor : rhbm_gem::BuiltInCommandCatalog())
    {
        if (descriptor.python_binding_name.empty())
        {
            std::cerr << "Missing python_binding_name for built-in command: "
                      << descriptor.name << '\n';
            return 1;
        }
        std::cout << descriptor.python_binding_name << '\n';
    }
    return 0;
}
