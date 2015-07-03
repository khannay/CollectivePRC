
#include "numPhaseOscillator_d.hh"  

using namespace std;

/* This file gives all the function definitions for the numPhaseOscillator class without using gsl */

numPhaseOscillator::numPhaseOscillator(double Rh) {
	Rhope=Rh;
	Beta=0.5;
	epsilon=0.1;
	gamma=0.5;
	generatefreq=0;

	h=0.1; //numerical timestep
	initvalues();
	findPhaseDist(); //find the phase dist
		


}//constructor


double numPhaseOscillator::Z_individual(double phi) {


       double twoPi = 2.0 * M_PI;
       double modphi=phi - twoPi * floor( phi / twoPi );
       double heavfactor;

       if (modphi < M_PI) {
	       heavfactor=0.0;
       } else {
	       heavfactor=1.0;
       }//else


       double x=heavfactor*(-sin(2*phi)-sin(2*phi)*cos(2*phi) );

       return(x);

}//Z_indiv





void numPhaseOscillator::initvalues(void) {

	K0=2.0*gamma/((1-pow(Rhope,2))*cos(Beta)); //hopefully this will give an equil R value close
	std::cout << "Coupling Constant  " << K0 << endl;


	for (int i=0; i< N; i++) {
		initial[i]=0.0;
	}//for

	//load the freq data
	ifstream input("Frequency_Data.txt",  ios::in);
    
	for (int i = 0; i < N; i++) {
		input >> wi[i];
		//std::cout << wi[i] << endl;
	}
	input.close();

}//init values



double numPhaseOscillator::heaviside(double x) {
	if ( x < 0.0) {
		return(0.0);
	} else {
		return(1.0);
	}
}//heaviside






inline void numPhaseOscillator::derv(double *y, double t, double *d) {
	/* This will only be used after the delta function is applied to integrate and find Delta_R */

	std::complex<double> Z(0, 0); //declare the order parameter as a complex number;
	std::complex<double> Inum(0,1); //define i as Inum

	for (int i=0; i< N; i++) {
		Z+=std::exp(Inum*y[i]);
	}
	double normf=(double) 1.0/N;
	Z=Z*normf;

	double R=std::abs(Z);
	double Psi=std::arg(Z);

	for (int i=0; i< N; i++) {
		d[i]=wi[i]+K0*R*sin(Psi-y[i]+Beta);
		//std::cout << d[i] << "\t";
	}//for
	//std::cout << "\n";

	//this should now run in O(2*N) time instead of O(N^2) a significant improvement

}//derv


void numPhaseOscillator::numDeltaZero(double *phi, double *phishifted) {

	for (int i=0; i< N; i++) {
		phishifted[i]=phi[i]+epsilon*Z_individual(phi[i]); //dumbest possible way to shift the oscillators
	}//for

}//numDeltaZero

	
	
void numPhaseOscillator::integrateSystem(double *initial, double tend, double *results) {
	/*implement a RK4 algorithm to compute the diff equations store the results in the passed two dimensional array results, each row is a time point */

	double d[N]; //declare the d array, to be used in passing to derv
	double y[N]; //state of the system
	double r[N]; //storage array used in the for loop
	double k1[N];
	double k2[N];
	double k3[N];
	double k4[N];
	std::copy(initial+0, initial+N, y); //copy the initial array into y
	double time=0.0;
	int timepoints= (int) tend/h; //number of timepoints (integer)
	for (int i=0; i < timepoints; i++) {

		for (int k=0; k < N; k++) {
			results[k+N*i]=y[k];
		}//copy for loop


		derv(y,time,k1);

		addarraysMe(y, k1, r, 0.5*h);
		derv(r,time+h/2.0,k2);

		addarraysMe(y,k2,r,0.5*h);
		derv(r,time+h/2.0,k3);

		addarraysMe(y, k3,r,h);
		derv(r,time+h,k4);

		RK4sum(k1,k2,k3,k4,y); //overwrite y with the RK4 sum

		time+=h;
	}//for timepoints
	

}//integrateSystem


void numPhaseOscillator::integrateSystemOP(double *initial, double tend, std::complex<double> *OP) {
	/*implement a RK4 algorithm to compute the diff equations store the results in the passed two dimensional array results, each row is a time point */

	double d[N]; //declare the d array, to be used in passing to derv
	double y[N]; //state of the system
	double r[N]; //storage array used in the for loop
	double k1[N];
	double k2[N];
	double k3[N];
	double k4[N];
	std::copy(initial+0, initial+N, y); //copy the initial array into y
	double time=0.0;
	int timepoints= (int) tend/h; //number of timepoints (integer)
	for (int i=0; i < timepoints; i++) {

		OP[i]=findOrderParameter(y);


		derv(y,time,k1);

		addarraysMe(y, k1, r, 0.5*h);
		derv(r,time+h/2.0,k2);

		addarraysMe(y,k2,r,0.5*h);
		derv(r,time+h/2.0,k3);

		addarraysMe(y, k3,r,h);
		derv(r,time+h,k4);

		RK4sum(k1,k2,k3,k4,y); //overwrite y with the RK4 sum

		time+=h;
	}//for timepoints
	

}//integrateSystemOP








void numPhaseOscillator::integrateSystemLast(double *initial, double tend, double *results) {

		/*implement a RK4 algorithm to compute the diff equations and only give back the final results after the full integration */

	double d[N]; //declare the d array, to be used in passing to derv
	double y[N]; //state of the system
	double r[N]; //storage array used in the for loop
	double k1[N];
	double k2[N];
	double k3[N];
	double k4[N];
	std::copy(initial+0, initial+N, y); //copy the initial array into y
	double time=0.0;
	int timepoints= (int) tend/h; //number of timepoints (integer)

	for (int i=0; i < timepoints; i++) {


		derv(y,time,k1);
		

		addarraysMe(y, k1, r, 0.5*h);//r=y+0.5*k1*h
                derv(r,time+h/2.0,k2); //get k2

		addarraysMe(y,k2,r,0.5*h); //r=y+0.5*h*k2
                derv(r,time+h/2.0,k3); //get k3

		addarraysMe(y, k3,r,h); //r=y+k3*h
		derv(r,time+h,k4);//get k4

		RK4sum(k1,k2,k3,k4,y); //overwrite y with the RK4 sum

		time+=h;
	}//for timepoints

	for (int i=0; i < N; i++) {
		results[i]=y[i]; //write the last y value to the results and only return that
	}//for

}//integrateSystemLast





void numPhaseOscillator::addarraysMe(double *a, double *b, double *r, double mul) {
	for (int i=0; i < N; i++) {
		r[i]=a[i]+mul*b[i];
	}//for

}//addArraysMe




void numPhaseOscillator::RK4sum(double *k1, double *k2, double *k3, double *k4, double *y) {
	for (int i=0; i < N; i++) {
		y[i]=y[i]+h/6.0*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i]);
	}//for
}//RK4sum




std::complex<double> numPhaseOscillator::findOrderParameter(double *results) {
	/*Find the complex order parameter given a vector of the phases of each oscillator N length array*/


	std::complex<double> orderP(0, 0); //declare the order parameter as a complex number;
	std::complex<double> Inum(0,1); //define i as Inum

	for (int i=0; i< N; i++) {
		orderP+=std::exp(Inum*results[i]);
	}
	double normf=(double) 1.0/N;
	orderP=orderP*normf;
	//R=std::abs(orderP);
	return(orderP);

}//findOrderParameter





void numPhaseOscillator::integrateTransients(double tend) {
	/*Integrate the system for a long time, then reset initial to the last state of the system */
	int timepoints= (int) tend/h;
	double results[N]; //init the results array
	integrateSystemLast(initial, tend, results);
	for (int i=0; i < N; i++) {
		initial[i]=results[i]; //reset the initial conditions to the final state
		//std::cout << initial[i] << endl;
	}//for
}//integrateTransients




void numPhaseOscillator::findPRCPoint(double phase, std::complex<double> *results) {

	double phase_dist_tmp[N];
	double delta0phi[N];
	double phase_dist_final[N];
	double phase_dist_unpre[N];

	for (int i=0; i< N; i++) {
		phase_dist_tmp[i]=phase_dist[i]+phase; //shift the mean phase of equil dist
	}//for

	numDeltaZero(phase_dist_tmp,delta0phi); //shift the system using the deltazero function

	double itime=1000.0;

	integrateSystemLast(delta0phi,itime,phase_dist_final); //integrate the system to get Delta_R results are stored in last
	integrateSystemLast(phase_dist_tmp, itime, phase_dist_unpre);



	std::complex<double> Zbar=findOrderParameter(delta0phi); //init shift

	std::complex<double> Zfinal=findOrderParameter(phase_dist_final);

	std::complex<double> Zfinalun=findOrderParameter(phase_dist_unpre);

	results[0]=Zbar;
	results[1]=Zfinal; 
	results[2]=Zfinalun;



}//findPRCPoint








void numPhaseOscillator::findPRC(void) {
	
	std::complex<double> normalZ=findOrderParameter(phase_dist);
	std::complex<double> Ival (0,1);

	//now compute the phase shift for diff time points in the oscillation
	int numSamples=100; //number of points to divide the period into
	double phases[numSamples]; //array to record the times or perturbations
	double phaseCoherence[numSamples];
	double meanPhase[numSamples];
	std::complex<double> Z0[numSamples];
	std::complex<double> Zbar[numSamples];
	std::complex<double> Zfinal[numSamples];
	std::complex<double> Zfinalun[numSamples];
	std::complex<double> Zans[3];

	double startPhase=0.0;
	double increment=2.0*M_PI/numSamples;
	for (int i=0; i < numSamples; i++) {
		phases[i]=startPhase;
		startPhase+=increment;
	}//for

#pragma omp parallel for
	for (int i=0; i < numSamples; i++) {
		Z0[i]=std::abs(normalZ)*exp(Ival*phases[i]);
		findPRCPoint(phases[i], Zans);
		Zfinal[i]=Zans[1];
		Zbar[i]=Zans[0];
		Zfinalun[i]=Zans[2];
	}//for

	/*Create a Report from this data */

	ofstream myfile;
	std::ostringstream sstream;
	sstream << Rhope;
	std::string varAsString = sstream.str();
	std::string filename1="PhaseData_light_"+varAsString+".txt"; //use std::to_string to add things to file name
	myfile.open(filename1.c_str(), ios::out);
	myfile << "Report from a run of Delta Infinity Finding Program. 0.5sinx-cos(x)\n";
	myfile << "Number of oscillators: " << N << endl;
	myfile << "Phase Coherence in locked state: " << std::abs(normalZ) << endl;
	myfile << "Coupling Term has Beta= " << Beta << endl;
	myfile << "Now the data: " << endl;
	myfile << "Z_0\t" << "Zbar\t" << "Zfinal\t" << endl;
	for (int i=0; i < numSamples; i++) {
		std::complex<double> s=Zfinal[i]/Zfinalun[i];
		myfile << Z0[i] << "\t" << Zbar[i] << "\t" << s << endl;
	}//for
	myfile.close();
	
}//findPRC




void numPhaseOscillator::findPhaseDist(void) {

	integrateTransients(3000); //integrate for a long time to get the equil phase dist


	std::complex<double> Z=findOrderParameter(initial); //find the OP of the final state
	std::cout << "The Equilbrium R value is: " << std::abs(Z) <<endl;
	
	std::complex<double> Ival (0,1); //declare the number i as Ival


	for (int i=0; i<N; i++) { 
		std::complex<double> tmp=exp(Ival*initial[i])/Z; 
		phase_dist[i]=std::arg(tmp); //find the angle for that complex number
	}//for


	
	//Try and make the collective frequency actually be zero
	double timeend=10.0;
	int timepoints= (int) timeend/h;
	double *results=new double[N*timepoints]; //one dimensional array to store all the results in (instead of a 2 dim array)
	integrateSystem(phase_dist,timeend,results); 

	//extract time series 
	double xtimes[timepoints];
	double xt=0.0;
	double meanPhase[timepoints];
	int counter=0;
	for (int i=0; i < N*timepoints; i+=N) {
		double tmpPhase[N];
		for (int k=0; k < N; k++) {
			tmpPhase[k]=results[k+i];
		}//for N
		meanPhase[counter]=std::arg(findOrderParameter(tmpPhase));
		xtimes[counter]=xt;
		//if (meanPhase[counter] < 0.0) { meanPhase[counter]+=2*M_PI; }
		xt=xt+h;
		counter+=1;
	}//for timepoints

	
	//convert the meanPhase data to absolute phase


	//find the collective frequency for the phase dist
	double slope=bestFitSlope(xtimes, meanPhase, timepoints);
	std::cout << "Slope estimate is: " << slope << endl;

	//adjust the system so the collective system will have collective freq approx zero
	for (int i=0; i < N; i++) {
		wi[i]=wi[i]-slope;
	}//for


	/* NOW WE MEASURE THE ERROR AND STORE THAT IN COLLECTIVE FREQUENCY */


	integrateSystem(phase_dist,timeend,results); 

	//extract time series 

	counter=0;
	for (int i=0; i < N*timepoints; i+=N) {
		double tmpPhase[N];
		for (int k=0; k < N; k++) {
			tmpPhase[k]=results[k+i];
		}//for N
		meanPhase[counter]=std::arg(findOrderParameter(tmpPhase));
		xtimes[counter]=xt;
		xt=xt+h;
		counter+=1;
	}//for timepoints

	
	//convert the meanPhase data to absolute phase


	//find the collective frequency for the phase dist
	double collectivefreq=bestFitSlope(xtimes, meanPhase, timepoints);
	std::cout << "Collective Freq is: " << collectivefreq << endl;







	/* Write the Phase Distribution to a File */
	ofstream myfile;
	std::ostringstream sstream;
	sstream << Rhope;
	std::string varAsString = sstream.str();
	std::string filename1="PhaseDistribution" +varAsString +".txt"; //use std::to_string to add things to file name
	myfile.open(filename1.c_str(), ios::out);
	for (int i=0; i < N; i++) {
		myfile << phase_dist[i] << "\t";
	}
	myfile << "\n";
	myfile.close();



}//find phase dist







void numPhaseOscillator::writePhaseData(void)  {
	/*write the phase oscillators dynamics to a file */

	double phase_dist_tmp[N];
	double phase=M_PI/2.0;

	for (int i=0; i< N; i++) {
		phase_dist_tmp[i]=phase_dist[i]+phase; //shift the mean phase of equil dist
	}//for

	//Integrate for a bit before the pertub is applied
	double timeendb=10.0;
	int timepointsb= (int) timeendb/h;
	double *resultsb=new double[N*timepointsb]; //one dimensional array to store all the results in (instead of a 2 dim array)
	integrateSystem(phase_dist_tmp,timeendb,resultsb); //results are written to 1d array results



	//apply the perturbation
	double shifted_pd[N];
	numDeltaZero(phase_dist_tmp, shifted_pd);




	//Integrate the system
	double timeend=20.0;
	int timepoints= (int) timeend/h;
	double *results=new double[N*timepoints]; //one dimensional array to store all the results in (instead of a 2 dim array)
	integrateSystem(shifted_pd,timeend,results); //results are written to 1d array results

	//Now write the data to a file

		 
	 ofstream myfile;
	 myfile.open ("phaseDynamics.txt");
	 //write the data from before the perturbation
	 int counter=1;
	 for (int i=0; i < N*timepointsb; i++) {
		 myfile << resultsb[i] << "\t";
		 counter+=1;
		 if (counter >= N) {
			 myfile << "\n";
			 counter=0;
		 }//if
	 }//for

	 //write the data from after the perturbation
	 counter=1;
	 for (int i=0; i < N*timepoints; i++) {
		 myfile << results[i] << "\t";
		 counter+=1;
		 if (counter >= N) {
			 myfile << "\n";
			 counter=0;
		 }//if
	 }//for
	 myfile.close();


}//write the phase data





void numPhaseOscillator::printMeanPhase(void) {

	double phase_dist_tmp[N];
	double phase=M_PI/2.0;

	for (int i=0; i< N; i++) {
		phase_dist_tmp[i]=phase_dist[i]+phase; //shift the mean phase of equil dist
	}//for


	//apply the perturbation
	double shifted_pd[N];
	numDeltaZero(phase_dist_tmp, shifted_pd);


	//Integrate the system
	double timeend=500.0;
	int timepoints= (int) timeend/h;
	std::complex<double> OP[timepoints]; //one dim array of the OP
	integrateSystemOP(shifted_pd,timeend,OP); //results are written to 1d array results

	


	//Now write the data to a file
	 ofstream myfile;
	 myfile.open ("meanPhase.txt");
	 double timec=0.0;
	 int counter=1;
	 for (int i=0; i < timepoints; i++) {
		 myfile << timec << "\t";
		 myfile << OP[i] << "\n";
		 timec+=0.1;
		 
	 }//for


	
	 myfile.close();


}//printMeanPhase









double numPhaseOscillator::bestFitSlope(double *x, double *y, int size) {
	
	double sumx=0.0;
	double sumy=0.0;
	double sumxy=0.0;
	double sumx2=0.0;

	for (int i=0; i < size; i++) {
		sumx+=x[i];
		sumy+=y[i];
		sumxy+=x[i]*y[i];
		sumx2+=x[i]*x[i];
	}//for

	double slope=(sumxy-(1.0/size)*sumx*sumy)/(sumx2-(1.0/size)*sumx*sumx);

	return(slope);

}//linearregression




int main(void) {

	
	double Rinit=0.50;
	while (Rinit <=0.95) {
		numPhaseOscillator a(Rinit);
		a.findPRC();
		Rinit+=0.05;
	}//while
	
	
	

	//numPhaseOscillator a(0.8);
	//a.findPRC();





}//main (Test)
