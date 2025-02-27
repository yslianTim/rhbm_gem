#include <vector>
#include <tuple>

class DataObjectBase;

struct AtomicPotentialEntry
{
    std::vector<std::tuple<float, float>> distance_and_map_value_list;
    std::tuple<double, double> gaus_estimate_posterior, gaus_variance_posterior;
    std::tuple<double, double> gaus_estimate_mdpde;
    std::tuple<double, double> gaus_estimate_ols;
    bool outlier_tag;
    double statistical_distance;
};

struct GroupPotentialEntry
{
    int classification_type;
    std::tuple<double, double> gaus_estimate_prior, gaus_variance_prior;
    std::tuple<double, double> gaus_estimate_mdpde;
    std::tuple<double, double> gaus_estimate_mean;
};