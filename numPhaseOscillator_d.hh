#include <iostream>
#include <cmath>
#include <fstream>
#include <math.h> 
#include <algorithm> //need std::copy to copy arrays
#include <math.h>
#include <complex>
#include <assert.h>
#include <omp.h>
#include <string>
#include <cstdlib> //used in command line input
#include <stdio.h>
#include <sys/time.h>
#include <cstring>

using namespace std;

/* To Do List:

g++ numPhaseOscillator_d.hh numPhaseOscillator_d.cpp -fopenmp  -fexceptions -O3
*/



class numPhaseOscillator {

public:
	static const int N=10000; //number of oscillators
	double Beta; //non-odd coupling part
	double epsilon;
	double gamma; //spread of the freq dist
	double Rhope; //pass this to try and generate a system with this equil R value
	double K0; //coupling constant
	double phase_dist[N]; //store the equilbrium phase dist
	double wi[N]; //store the frequency dist of each node
	double initial[N];
	double collectivefreq;
	int generatefreq; //whether to load the freq or to generate them randomly 1=generate, 0=load

	//Numerical Parameters
	double h; //this the the timestep in the numerical integration alg



	/*Methods for this class */

	numPhaseOscillator(double Rh); //pass it a hopefull R value
	double Z_individual(double phi); //the individual PRC

	void initvalues(void); //init the frequency array, and initial array
	inline void derv(double *y, double t, double *d); //compute the derivative, write the answer on the array d, given array y #
	double heaviside(double x); //#
	std::complex<double> findOrderParameter(double *results); //#
	void findPRCPoint(double phase, std::complex<double> *results); //given a starting point, the period the strength and the start time of pertub find the phase shift of the mean phase
	void numDeltaZero(double *phi, double *phishifted); //given a phase dist implement a delta shift

	void findPRC(void); //find the collective PRC and output the results to a file
	
	void addarraysMe(double *a, double *b, double *r, double mul); //r=a+b*mul #
	void RK4sum(double *k1, double *k2, double *k3, double *k4, double *y); //compute the last sum in RK4, overwrite y #
	void integrateSystem(double *initial, double tend, double *results); //integrate the system from the initial condition from t=0 to tend, write the results into the results 2DIM array #
	void integrateSystemLast(double *initial, double tend, double *results); //integrate the system from the initial condition and only report the last state of the system at tend
	void integrateTransients(double tend); //reset initial conditions to final state
	void integrateSystemOP(double *initial, double tend, std::complex<double> *OP); //find OP as a function of time
	void findPhaseDist(void); //integrate the system to find the equil phase dist then normalize to have mean phase=0



	void writePhaseData(void);
	double bestFitSlope(double *x, double *y, int size);
	void printMeanPhase(void); //given a time series produce a file with time series of order parameter

}; //class numCollective
