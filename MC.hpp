#include <vector>

class MC {

private:
    std::vector<double> centroid;
    int lastUpdate;
    double weight;

public:
    MC(std::vector<double> centroid, int lastUpdate, double weight);
    MC(std::vector<double> centroid, int lastUpdate);
    std::vector<double> getCentroid();
    double getWeight();
    double distance(MC mc);
    double distance(std::vector<double> &x);
    void merge(MC mc, int t, double lambda, double r);
    void fade(int t, double lambda);

};


