// Standalone Release arithmetic replay of the frozen coverage distance loop.
// Compile with the experiment's compiler and -O3 -std=gnu++20, adding -Isrc.
#include "core/command/detail/SimulationGeometry.hpp"
#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using Position = std::array<double, 3>;
// Original ObservationMatchedPrediction.cpp::SquareDistance, unchanged at ecf55f45.
__attribute__((noinline)) double LegacySquare(const Position & a, const Position & b)
{
    double out{}; for (std::size_t k=0;k<3;++k) out += (a[k]-b[k])*(a[k]-b[k]); return out;
}
std::vector<std::vector<double>> Read(const char * path)
{
    std::ifstream in(path); if(!in) throw std::runtime_error("Missing input table.");
    std::string line; std::getline(in,line); std::vector<std::vector<double>> rows;
    while(std::getline(in,line))
    {
        std::istringstream stream(line); std::string cell; auto & row=rows.emplace_back();
        while(std::getline(stream,cell,',')) row.push_back(std::stod(cell));
    }
    return rows;
}
int main(int argc, char ** argv)
{
    if(argc!=4) return 2;
    const auto voxels=Read(argv[1]),atoms=Read(argv[2]),contributors=Read(argv[3]);
    std::size_t index{},memberships{},membership_difference{},square_difference{},csr_difference{};
    for(std::size_t row=0;row<voxels.size();++row)
    {
        const Position p{voxels[row][2],voxels[row][3],voxels[row][4]};
        for(std::size_t atom=0;atom<atoms.size();++atom)
        {
            const Position a{atoms[atom][0],atoms[atom][1],atoms[atom][2]};
            const double legacy=LegacySquare(p,a),current=rhbm_gem::core::simulation::SupportSquare(p,a);
            membership_difference+=(legacy<=6.25)!=(current<=6.25);
            square_difference+=legacy!=current;
            if(legacy<=6.25)
            {
                ++memberships;
                if(index>=contributors.size() || contributors[index][0]!=static_cast<double>(row) ||
                    contributors[index][1]!=static_cast<double>(atom) || contributors[index][2]!=current) ++csr_difference;
                ++index;
            }
        }
    }
    csr_difference+=index!=contributors.size();
    std::cout<<"{\"rows\":"<<voxels.size()<<",\"memberships\":"<<memberships
        <<",\"membership_differences\":"<<membership_difference<<",\"square_differences\":"<<square_difference
        <<",\"csr_differences\":"<<csr_difference<<"}\n";
    return membership_difference || square_difference || csr_difference ? 1 : 0;
}
