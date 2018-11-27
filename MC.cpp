#include <math.h>
#include "MC.hpp"

MC::MC(std::vector<double> centroid, int lastUpdate, double weight) {
    this->centroid=centroid;
    this->lastUpdate=lastUpdate;
    this->weight=weight;
}


MC::MC(std::vector<double> centroid, int lastUpdate) {
    this->centroid=centroid;
    this->lastUpdate=lastUpdate;
    this->weight=1;
}

double MC::getWeight(){
	return(this->weight);
}

std::vector<double> MC::getCentroid(){
	return(this->centroid);
}

void MC::merge(MC mc, int t, double lambda, double r) {
	mc.fade(t, lambda);
	this->fade(t, lambda);

	// update statistics
	this->weight += mc.weight;

	// competetive learning
	double d = this->distance(mc);

	std::vector<double> mcCentroid = mc.getCentroid();
	for(unsigned int i=0; i<this->centroid.size(); i++){
	   this->centroid[i] += exp(-pow(d/r*3.0, 2.0) /2.0) * (mcCentroid[i]-this->centroid[i]);
	}
}


void MC::fade(int t, double lambda){
	// apply fading
	this->weight *= pow(2,(-lambda * (t-this->lastUpdate)));
	// update time
	this-> lastUpdate = t;
}

double MC::distance(MC mc){
	std::vector<double> thisCentre = this->getCentroid();
	std::vector<double> mcCentre = mc.getCentroid();
	double sum = 0.0;
	for(unsigned int i=0; i<thisCentre.size(); i++){
	  sum += pow(thisCentre[i] - mcCentre[i], 2);
	}
	return(sqrt(sum));
}

double MC::distance(std::vector<double> &x){
	std::vector<double> thisCentre = this->getCentroid();
	double sum = 0.0;
	for(unsigned int i=0; i<thisCentre.size(); i++){
	  sum += pow(thisCentre[i] - x[i], 2);
	}
	return(sqrt(sum));
}