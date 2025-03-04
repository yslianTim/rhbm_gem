#include <iostream>
#include <memory>
#include "DataObjectManager.hpp"
#include "PotentialAnalysisVisitor.hpp"
#include "AtomSelector.hpp"
#include "SphereSampler.hpp"

using std::cout;
using std::endl;

int main(int argc, char* argv[])
{
    for (int i = 0; i < argc; i++)
    {
        std::cout << argv[i] << std::endl;
    }

    auto data_manager{ std::make_unique<DataObjectManager>("database.sqlite") };
    data_manager->ProcessFile("/Users/yslian/data/6z6u.cif","model_1");
    data_manager->ProcessFile("/Users/yslian/data/emd_11103_unsharpened.map","map_1");
    data_manager->PrintDataObjectInfo("model_1");

    auto atom_selector{ std::make_shared<AtomSelector>() };
    atom_selector->PickChainID("A,AA,BA");
    atom_selector->PickIndicator(".,A");

    auto sphere_sampler{ std::make_shared<SphereSampler>() };
    sphere_sampler->SetSamplingSize(100);
    sphere_sampler->SetDistanceRangeMinimum(0.0);
    sphere_sampler->SetDistanceRangeMaximum(1.5);

    auto analyzer{ std::make_unique<PotentialAnalysisVisitor>(atom_selector, sphere_sampler) };
    analyzer->SetMapObjectKeyTag("map_1");
    analyzer->SetModelObjectKeyTag("model_1");
    analyzer->SetFitRange(0.0, 1.0);
    analyzer->SetAlphaR(0.1);
    analyzer->SetAlphaG(0.2);

    data_manager->Accept(analyzer.get());

    data_manager->SaveDataObject("model_1");

    return 0;
}