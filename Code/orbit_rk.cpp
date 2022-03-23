/*

Orbital motion
Runge-Kutta 

Raphael Van Laer

10-23-2021.

*/


#include <chrono>
#include <iostream>
#include <cmath>
#include <fstream>


using namespace std;
using namespace std::chrono;


//Global Variables:
double PI = acos(-1.0);
double GM = 4 * PI * PI;


int rk4( double *q, double t, double tau, 
    int (*derfunc)(double *q, double t, double *par, double *derive),   // This passes a function into another function
    double *par){

	// Somewhat general RK4 update for two dimensional problem

	double *k1, *k2, *k3, *k4;
	double *qold, *qtmp;

	k1 = new double[4];
	k2 = new double[4];
	k3 = new double[4];
	k4 = new double[4];
	qold = new double[4];
	qtmp = new double[4];


	// can't do: qold = q
	for(int i = 0; i < 4; i++){
		qold[i] = q[i];
		qtmp[i] = q[i];
	}

	double t12;

	// Define k1:
	int c1 = derfunc(q, t, par, k1);
	t12 = t + tau / 2.0;

	for(int i = 0; i < 4; i++){
		qtmp[i] = q[i] + tau * k1[i] / 2.0;
	}

	// Define k2:
	int c2 = derfunc(qtmp, t12, par, k2);

	for(int i = 0; i < 4; i++){
		qtmp[i] = q[i] + tau * k2[i] / 2.0;
	}

	// Define k3:
	int c3 = derfunc(qtmp, t12, par, k3);

	for(int i = 0; i < 4; i++){
		qtmp[i] = q[i] + tau * k3[i];
	}

	// Define k4:
	int c4 = derfunc(qtmp, t + tau, par, k4);

	for(int i = 0; i < 4; i++){
		q[i] = qold[i] + (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) * tau/6;
	}


	// Deallocate memory for this array

	delete[] k1;
	delete[] k2;
	delete[] k3;
	delete[] k4;
	delete[] qold;
	delete[] qtmp;


	return c1 + c2 + c3 + c4;
}


int gravity(double *q, double t, double *par, double *deriv){
	// This is my "derfunc" for GM/r^2
	// par: array of parameters; par is a 1-d array
	// with par[0] = GM

	double x, y, r;

	x = q[0];
	y = q[1];
	r = sqrt(x*x + y*y);

	deriv[0] = q[2];						// Derivative of q[0] --> vx --> q[2]
	deriv[1] = q[3];
	deriv[2] = -(par[0]/pow(r,3)) * x;				// derivative of q[2] --> a_x
	deriv[3] = -(par[0]/pow(r,3)) * y;

	return 0;
}


int main(int argc, char *argv[]){

	auto start = high_resolution_clock::now();

	double tau;
	double t = 0.0;
	double x0, y0, vx0, vy0;
	double *q, *q0;							// pointer declaration
	double *par;

	q = new double[4];						// new: new array, double is type, 2 is number of element
	q0 = new double[4];
	par = new double[1];
	par[0] = GM;
	
	if(argc < 6){ 

		cout << "Usage: " << argv[0] << " x0 y0 vx0 vy0 tau" << endl;
		return -1;
	}
	
	else{

		q0[0] = atof(argv[1]);				
		q0[1] = atof(argv[2]);				
		q0[2] = atof(argv[3]);
		q0[3] = atof(argv[4]);
		tau = atof(argv[5]);
	}

	// Initialize arrays
	for(int i = 0; i <4; i++){
		q[i] = q0[i];
	}

	double energy;						// energy / mass

	energy = 0.5 * (q[2]*q[2] + q[3]*q[3]) - GM/sqrt(q[0]*q[0]+q[1]*q[1]);

	double AnglMoment;

	AnglMoment = -1.0 * q[1]*q[2] + q[0]*q[3];

	ofstream fout ; 				
	fout.open("orbit_rk.dat"); 			

	fout << t << "\t" << q[0] << "\t" << q[1];
	fout << "\t" << q[2] << "\t" << q[3] << "\t" << energy << "\t" << AnglMoment << endl;	


//	int start = 0;
	int nmax = 10 / tau;
	int count = 0;

	double period, t0 = 0.0;
								// I want to measure the period 10 times.

	while( count < nmax ){

		rk4( q, t, tau, gravity, par);

		t = t + tau;			

		energy = 0.5 * (q[2]*q[2] + q[3]*q[3]) - GM/sqrt(q[0]*q[0]+q[1]*q[1]);

		AnglMoment = -1.0 * q[1]*q[2] + q[0]*q[3];

		// PRINTING OUT

		fout << t << "\t" << q[0] << "\t" << q[1];
		fout << "\t" << q[2] << "\t" << q[3] << "\t" << energy << "\t" << AnglMoment << endl;	
	
		count++;
	}

	fout.close();

	FILE *gnuplot = popen("gnuplot", "w");			// opens a pipe
	fprintf(gnuplot, "set out 'orbit_rk.ps'\n ");		// Opens a window
	fprintf(gnuplot, "set term post land\n ");		// using a poster file in landscape
	fprintf(gnuplot, "set size square\n");
	fprintf(gnuplot, "set title 'Circular, Orbital Motion in RK4'\n");
	fprintf(gnuplot, "set xlabel 'x (m)'\n");
	fprintf(gnuplot, "set ylabel 'y (m)'\n");
	fprintf(gnuplot, "plot 'orbit_rk.dat' u 2:3 w l title'orbit'\n");
	fprintf(gnuplot, "set title 'E/m vs time for Circular Orbit in RK4'\n");
	fprintf(gnuplot, "set xlabel 't'\n");
	fprintf(gnuplot, "set ylabel 'E/m'\n");
	fprintf(gnuplot, "plot 'orbit_rk.dat' u 1:6 w l\n"); 
								// plot the data from the file using column 2 and 3 with lines
	pclose(gnuplot);					// closes the plotting

	system( (char *) "ps2pdf orbit_rk.ps");		// transforms the plotting into a pdf

	delete[] q;
	delete[] q0;

	auto stop = high_resolution_clock::now();

	auto duration = duration_cast<microseconds>(stop - start);

	cout << "Time taken by Algorithm = " << duration.count() << " microseconds." << endl;

	return 0;

}	