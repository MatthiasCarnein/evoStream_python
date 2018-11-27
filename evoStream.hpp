#include <iostream>
#include <vector>
#include <random>

#include "MC.hpp"


class EvoStream {

private:
  double r;
  double lambda;
  int tgap;
  unsigned int k;
  double crossoverRate;
  double mutationRate;
  unsigned int populationSize;
  unsigned int initializeAfter;
  int reclusterGenerations;
  double omega;
  int t;
  int init;
  int upToDate;
  double delay;
  int initTime;
  std::vector<MC> micro;
  std::vector<std::vector< std::vector<double> > > macro; // [solution][cluster][centre]
  std::vector<double> macroFitness;
  std::mt19937 rng;


  void updateWeights();
  void removeMicroCluster(int i);
  void insert(std::vector<double> &distances, MC mc);
  void evolution();
  void calculateFitness();
  double fitness(std::vector< std::vector<double> > &centres);
  void initialize();
  std::vector<std::vector< std::vector<double> > > selection();
  std::vector<std::vector< std::vector<double> > > recombination(std::vector<std::vector< std::vector<double> > > &individuals);
  std::vector<std::vector< std::vector<double> > > mutation(std::vector<std::vector< std::vector<double> > > &individuals);
  std::vector<std::vector< std::vector<double> > > recombinationGAclust(std::vector<std::vector< std::vector<double> > > &individuals);
  std::vector<std::vector< std::vector<double> > > recombinationPESAII(std::vector<std::vector< std::vector<double> > > &individuals);
  std::vector<std::vector< std::vector<double> > > mutationGAclust(std::vector<std::vector< std::vector<double> > > &individuals);
  std::vector<std::vector< std::vector<double> > > mutationPESAII(std::vector<std::vector< std::vector<double> > > &individuals);
  
  int ndimensions();
  int sampleProportionally(std::vector<double> &data);
  std::vector<double> getDistanceVector(MC mc, std::vector<MC> &cluster);

public:   
  EvoStream(double r, double lambda, int tgap, unsigned int k, double crossoverRate, double mutationRate, int populationSize, unsigned int initializeAfter, int reclusterGenerations);
  std::vector< std::vector<double> > get_microclusters();
  std::vector<double> get_microweights();
  std::vector< std::vector<double> > get_macroclusters();
  std::vector<double> get_macroweights();
  std::vector<int> microToMacro();
  void recluster(int generations);
  std::vector<int> getAssignment(std::vector< std::vector<double> > &centres);
  double getMaxFitness();
  void cluster(std::vector<double> &data);
  void cleanup();

};

