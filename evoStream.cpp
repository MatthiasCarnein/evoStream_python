#include <algorithm>
#include <ctime>
#include <limits>


#include "evoStream.hpp"


#define VERBOSE 0


EvoStream::EvoStream(double r, double lambda, int tgap, unsigned int k, double crossoverRate, double mutationRate, int populationSize, unsigned int initializeAfter, int reclusterGenerations){
	this->r=r;
	this->lambda=lambda;
	this->tgap=tgap;
	this->reclusterGenerations = reclusterGenerations;
	this->k=k;
	this->crossoverRate=crossoverRate;
	this->mutationRate=mutationRate;
	this->populationSize=populationSize;
	this->initializeAfter=initializeAfter;
	this->init=0;

	this->macroFitness = std::vector<double>(this->populationSize,0);
	this->omega = pow(2, (-1*lambda * tgap));
	this->t=0;
	this->upToDate=0;

	// random_device not working on all systems, fall back to old C style seeding
	// std::random_device rd; 
	// this->rng = std::mt19937(rd());  
	this->rng = std::mt19937(time(NULL));
 }




////////////////////////// Interfacer


std::vector< std::vector<double> > EvoStream::get_microclusters(){
	int d = this->ndimensions();
	std::vector< std::vector<double> > x(this->micro.size(),std::vector<double>(d));
	for(unsigned int i=0; i<this->micro.size();i++){
	  std::vector<double> mc = this->micro[i].getCentroid();
	  for(unsigned int j=0; j<mc.size();j++){
		x[i][j]=mc[j];
	  }
	}
	return(x);
}


std::vector<double> EvoStream::get_microweights(){
	std::vector<double> x(this->micro.size());
	for(unsigned int i=0; i<this->micro.size();i++){
	  x[i]=this->micro[i].getWeight();
	}
	return(x);
}



std::vector< std::vector<double> > EvoStream::get_macroclusters(){

	if(!this->init){
		return(std::vector<std::vector<double> >(0));
	}

	if(this->reclusterGenerations!=0 && this->upToDate==0){
	  this->recluster(reclusterGenerations);
	  this->upToDate=1;
	}

	// find max fitness
	int maxIdx=-1;
	double max=-1*std::numeric_limits<double>::max();

	for(unsigned int i=0; i<this->macroFitness.size(); i++){
	  if(this->macroFitness[i]>max){
		max=this->macroFitness[i];
		maxIdx = i;
	  }
	}
	return(this->macro[maxIdx]);
}



std::vector<double> EvoStream::get_macroweights(){

	if(!this->init){
		return(std::vector<double>(0));
	}

	if(reclusterGenerations!=0 && this->upToDate==0){
	  this->recluster(reclusterGenerations);
	  this->upToDate=1;
	}

	std::vector<int> clusterAssignment = this->microToMacro();
	std::vector<double> microWeights = this->get_microweights();
	std::vector<double> macroWeights(this->k);

	// for every cluster
	for(unsigned int i=0; i<clusterAssignment.size(); i++){
	  macroWeights[clusterAssignment[i]] += microWeights[i];
	}

	return(macroWeights);
}


std::vector<int> EvoStream::microToMacro(){

	if(!this->init){
		return(std::vector<int>(0));
	}

	std::vector<std::vector<double> > centres = this->get_macroclusters();
	std::vector<int> assignment = this->getAssignment(centres);
	return(assignment);
}







////////////////////////////// Online component



void EvoStream::cluster(std::vector<double> &data){

	this->upToDate=0;

	// increment time
	this->t++;

	// create temporary mc
	MC mc(data, this->t);

	// get distance to all other mcs
	std::vector<double> distances = this->getDistanceVector(mc, this->micro);

	// insert into closest or init new mc
	this->insert(distances, mc);

	// cleanup every tgap
	if(this->t % this->tgap == 0){
	  this->cleanup();
	}

	if(!this->init && this->micro.size()==this->initializeAfter){
	  this->initialize();
	}

}






void EvoStream::cleanup(){

	#if VERBOSE >= 2
	std::cout << "Cleanup" << std::endl;
	#endif
	this->updateWeights();

	// merge close mcs
	for(int i = this->micro.size()-1; i >=0; i--){
	  for(int j = i-1; j>=0; j--){
		double d = this->micro[i].distance(this->micro[j]); // calc distance
		if(d <= r){
	#if VERBOSE >= 2
		  std::cout << "Merging cluster" << i << " into " << j << std::endl;
	#endif
		  this->micro[j].merge(this->micro[i], this->t, this->lambda, this->r); //merge mc into the other
		  this->removeMicroCluster(i); // remove the old mc
		  break;
		}
	  }
	}
}

void EvoStream::updateWeights(){

	// fade Clusters
	for(int i=micro.size()-1; i>=0; i--){
	  this->micro[i].fade(this->t, this->lambda);
	  // remove insufficient weight
	  if( this->micro[i].getWeight() <= this->omega){
		this->removeMicroCluster(i);
	  }
	}
}



void EvoStream::removeMicroCluster(int i){
	#if VERBOSE >= 2
	std::cout << "Remove Cluster " << i << std::endl;
	#endif
	this->micro.erase(this->micro.begin()+i);
}



void EvoStream::insert(std::vector<double> &distances, MC mc){

	// merge if within radius
	bool merged = 0;
	for(unsigned int i=0; i < micro.size(); i++){
	  if(distances[i] <= r){
	#if VERBOSE >= 2
		std::cout << "Merge observation " << t << " into Micro Cluster " << i  << std::endl;
	#endif
		this->micro[i].merge(mc, this->t, this->lambda, r);
		merged=1;
	  }
	}

	// if none close enough create new
	if(!merged){
	#if VERBOSE >= 2
	  std::cout << "Use observation " << t <<" to create new Micro Cluster" << std::endl;
	#endif
	  this->micro.push_back(mc);
	}
}

















//////////////////////////////////// Offline component




void EvoStream::evolution(){
	if(!this->init) return;

	this->calculateFitness();

	// perform evolutionary step
	std::vector<std::vector< std::vector<double> > > selected = this->selection(); // select parents
	std::vector<std::vector< std::vector<double> > > offsprings = this->recombination(selected); // recombine parents
	std::vector<std::vector< std::vector<double> > > mutants = this->mutation(offsprings); // mutate offspring




	for(unsigned int i=0; i<mutants.size(); i++){
	  double fit = this->fitness(mutants[i]); // evaluate new solution

	#if VERBOSE >= 1
	  std::cout << "Fitness of child: " << fit << std::endl;
	#endif

	  // replace weakest solution so far
	  int minIdx = 0;
	  double min = std::numeric_limits<double>::max();
	  for(unsigned int j=0; j<this->macroFitness.size(); j++){
		if(this->macroFitness[j] < min){
		  min = this->macroFitness[j];
		  minIdx = j;
		}
	  }

	  if(min < fit){
	#if VERBOSE >= 1
		std::cout << "Use new solution to replace solution " << minIdx << " with fitness " << this->macroFitness[minIdx] << std::endl;
	#endif

		this->macro[minIdx] = mutants[i];
		this->macroFitness[minIdx] = fit;
	  }
	}


	#if VERBOSE >= 1
	std::cout << "Highest Fitness is  " <<  this->getMaxFitness() << std::endl;
	#endif

}



void EvoStream::recluster(int generations){
	if(!this->init) return;

	// fixed number of generations
	for(int i=0; i<generations; i++){
	#if VERBOSE >= 1
	  std::cout << "Generation " << i << std::endl;
	#endif
	  this->evolution();
	}
}




void EvoStream::calculateFitness(){
	// for every solution
	for(unsigned int i=0; i<this->macro.size(); i++){
	  this->macroFitness[i] = this->fitness(this->macro[i]);
	}
}



double EvoStream::fitness(std::vector< std::vector<double> > &centres){
	double result=0.0;

	std::vector<int> assignment = this->getAssignment(centres);

	// calc distance to centers
	for(unsigned int i=0; i<assignment.size(); i++){
	  result += pow(this->micro[i].distance(centres[assignment[i]]),2) * this->micro[i].getWeight();
	}
	return(1/result);
}



void EvoStream::initialize(){
	#if VERBOSE >= 2
	std::cout << "Initialize Solutions" << std::endl;
	#endif
	this->initTime = this->t;

	// for entire population
	for(unsigned int i=0; i<this->populationSize; i++){

	  //init
	  std::vector< std::vector<double> > temp = std::vector< std::vector<double> >(this->k, std::vector<double>(this->ndimensions()));
	  macro.push_back(temp);

	  // choose k centers
	  std::vector<int> choose(this->micro.size());
	  for(unsigned int j=0; j<this->micro.size(); j++){
		choose[j]=j;
	  }

	  std::shuffle(choose.begin(), choose.end(), this->rng);

	  for(unsigned int j=0; j<this->k; j++){
		int l = j % this->micro.size();
		macro[i][j] = this->micro[choose[l]].getCentroid();
	  }
	}

	this->init=1;
}


std::vector<std::vector< std::vector<double> > > EvoStream::selection(){

	std::vector<std::vector< std::vector<double> > > individuals;
	individuals.reserve(2);

	// Select proportionally to fitness
	double sum=0.0;
	std::vector<double> probability(this->macroFitness.size());
	for(unsigned int i=0; i<this->macroFitness.size(); i++){
	  sum += this->macroFitness[i];
	  probability[i] = this->macroFitness[i];
	}

	// sample two parents
	std::vector<int> selected(2);
	for(int i=0; i<2; i++){
	std::uniform_real_distribution<double> unif(0.0,1.0);
	  double rand = unif(this->rng);

	  for(unsigned int j=0; j<probability.size(); j++){

		double val = probability[j]/sum;

		if(val > rand){
		  selected[i]=j;
		  sum -= probability[j];
		  probability.erase(probability.begin()+j);
		  break;
		}
		rand -= val;
	  }
	}
	// since we deleted one element
	if(selected[1] >= selected[0]){
	  selected[1]++;
	}

	#if VERBOSE >= 1
	std::cout << "Select " << selected[0] << " and " << selected[1] << " with fitness " << this->macroFitness[selected[0]] << " and " << this->macroFitness[selected[1]] << " as parents" << std::endl;
	#endif

	individuals.push_back(this->macro[selected[0]]);
	individuals.push_back(this->macro[selected[1]]);

	return(individuals);
}


// wrapper for easier change of recombination algorithm
std::vector<std::vector< std::vector<double> > > EvoStream::recombination(std::vector<std::vector< std::vector<double> > > &individuals){
	return(recombinationGAclust(individuals));
}

// wrapper for easier change of mutation algorithm
std::vector<std::vector< std::vector<double> > > EvoStream::mutation(std::vector<std::vector< std::vector<double> > > &individuals){
	return(mutationGAclust(individuals));
}


// GA Clustering recombination approach
std::vector<std::vector< std::vector<double> > > EvoStream::recombinationGAclust(std::vector<std::vector< std::vector<double> > > &individuals){

	std::uniform_real_distribution<double> unif(0.0,1.0);

	if(unif(this->rng) < this->crossoverRate){
	  unsigned int nrow =  individuals[0].size();
	  unsigned int ncol = individuals[0][0].size();
	  unsigned int size = nrow * ncol;
	  int crossoverPoint = unif(this->rng)*(size-1); // range 1 - length-1
	#if VERBOSE >= 2
	  std::cout << "CrossoverPoint: " << crossoverPoint << std::endl;
	#endif

	  int pos=0;
	  for(unsigned int row=0; row < nrow; row++){
		for(unsigned int col=0; col<ncol; col++){
		  if(pos>crossoverPoint){
			double temp = individuals[0][row][col];
			individuals[0][row][col] = individuals[1][row][col];
			individuals[1][row][col] = temp;
		  }
		  pos++;
		}
	  }
	}

	return(individuals);
}


// PESAII reclustering approach
std::vector<std::vector< std::vector<double> > > EvoStream::recombinationPESAII(std::vector<std::vector< std::vector<double> > > &individuals){

	std::uniform_real_distribution<double> unif(0.0,1.0);


	for(unsigned int i=0; i<individuals.size(); i++){
	  unsigned int nrow =  individuals[i].size();
	  unsigned int ncol = individuals[i][0].size();
	  for(unsigned int row=0; row < nrow; row++){
		for(unsigned int col=0; col < ncol; col++){
		  if(unif(this->rng) < this->crossoverRate){
			double temp = individuals[0][row][col];
			individuals[0][row][col] = individuals[1][row][col];
			individuals[1][row][col] = temp;
		  }
		}
	  }
	}

	return(individuals);
}



// GA Clustering mutation approach
std::vector<std::vector< std::vector<double> > > EvoStream::mutationGAclust(std::vector<std::vector< std::vector<double> > > &individuals){

	std::uniform_real_distribution<double> unif(0.0,1.0);

	// for both individuals
	for(unsigned int i=0; i<individuals.size(); i++){
	  unsigned int nrow =  individuals[i].size();
	  unsigned int ncol = individuals[i][0].size();
	  for(unsigned int row=0; row < nrow; row++){
		for(unsigned int col=0; col<ncol; col++){
		  if(unif(this->rng) < this->mutationRate){
			double val=0.0;

			if(individuals[i][row][col]!=0){
			  val = 2 * unif(this->rng) * individuals[i][row][col];
			} else{
			  val = 2 * unif(this->rng);
			}

			// + or - with equal change
			if(unif(this->rng)>0.5){
			  individuals[i][row][col] += val;
			} else{
			  individuals[i][row][col] -= val;
			}
		  }
		}
	  }
	}

	return(individuals);
}

// PESAII mutation approach
std::vector<std::vector< std::vector<double> > > EvoStream::mutationPESAII(std::vector<std::vector< std::vector<double> > > &individuals){

	std::uniform_real_distribution<double> unif(0.0,1.0);
	std::normal_distribution<double> norm(0,.3); 

	// for both individuals
	for(unsigned int i=0; i<individuals.size(); i++){
	  unsigned int nrow =  individuals[i].size();
	  unsigned int ncol = individuals[i][0].size();
	  for(unsigned int row=0; row < nrow; row++){
		for(unsigned int col=0; col<ncol; col++){
		  if(unif(this->rng) < this->mutationRate){
			individuals[i][row][col] += norm(this->rng);
		  }
		}
	  }
	}

	return(individuals);
}






///////////////////////// Helper

std::vector<int> EvoStream::getAssignment(std::vector< std::vector<double> > &centres){
	std::vector<int> assignment(this->micro.size());
	// for every cluster
	for(unsigned int i=0; i<this->micro.size(); i++){
	  double min=std::numeric_limits<double>::max();
	  // find closest centre
	  for(unsigned int j=0; j<centres.size(); j++){
		double dist = this->micro[i].distance(centres[j]);
		if(dist < min){
		  min=dist;
		  assignment[i]=j;
		}
	  }
	}
	return(assignment);
}



double EvoStream::getMaxFitness(){
	// find max fitness
	if(!this->init){
	  return(0);
	}

	double max=-1*std::numeric_limits<double>::max();
	for(unsigned int i=0; i<this->macroFitness.size(); i++){
	  if(this->macroFitness[i]>max){
		max=this->macroFitness[i];
	  }
	}
	return(max);
}



int EvoStream::ndimensions(){
	if(!micro.size()){
	  return(0);
	} else{
	  return(this->micro[0].getCentroid().size());
	}
}



int EvoStream::sampleProportionally(std::vector<double> &data){

	std::uniform_real_distribution<double> unif(0.0,1.0);

	unsigned int j=0;

	// sample proportionally to distance
	double cumsum=0.0;
	for(j=0; j<data.size(); j++){
	  cumsum += data[j];
	  if(cumsum >= unif(this->rng)){
		return(j);
	  }
	}
	return(-1);
}



std::vector<double> EvoStream::getDistanceVector(MC mc, std::vector<MC> &cluster){

	std::vector<double> distances(cluster.size());
	for(unsigned int i=0; i < cluster.size(); i++){
	  distances[i] = mc.distance(cluster[i]);
	}
	return(distances);
}