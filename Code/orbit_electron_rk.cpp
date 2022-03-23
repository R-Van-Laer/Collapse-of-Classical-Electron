/*

Collapse of Classical Electron 
Runge-Kutta 4 Method

Raphael Van Laer

10-23-2021.

*/


#include <iostream>
#include <cmath>
#include <fstream>


using namespace std;


//Global Variables:
double PI = acos(-1.0);
double GM = 4 * PI * PI;
double c = 3.0;
double a0 = 2.8;

int rk4( double *q, double t, double tau, 
    int (*derfunc)(double *q, double t, double *par, double *derive),
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

	double r;

	r = q[0];

	deriv[0] = q[2];
	deriv[1] = q[3];
	deriv[2] = -32.0/9.0 * (par[0]*par[0]*par[0]*par[0])/(r*r*r*r*r) * par[1]*par[1] - par[0]/(r*r) * par[1]*par[1];
	deriv[3] = -8.0/3.0 * sqrt(par[0]/r) * (par[0]*par[0])/(r*r*r) * par[1]*par[1];

	return 0;
}


int main(int argc, char *argv[]){

	double tau;
	double t = 0.0;
	double *q, *q0;							// pointer declaration
	double *par;

	q = new double[4];						// new: new array, double is type, 2 is number of element
	q0 = new double[4];
	par = new double[2];
	par[0] = a0;
	par[1] = c;
	
	if(argc < 4){ 

		cout << "Usage: " << argv[0] << " r0 Theta0 tau" << endl;
		cout << "The radius will be in fm, meaning 1 = 10^-15 m ." << endl;
		cout << "The angle will be in radians." << endl;
		cout << "The time will be in RVLs (invented second variant), meaning 1 = 10^-23 s." << endl;
		return -1;
	}
	
	else{

		q0[0] = atof(argv[1]);				
		q0[1] = atof(argv[2]);				
		q0[2] = 0.0;
		q0[3] = sqrt(par[0]/q0[0]) * par[1];
		tau = atof(argv[3]);
	}

	for(int i = 0; i < 4; i++){
		q[i] = q0[i];
	}

	ofstream fout ; 				
	fout.open("orbit_electron_rk.dat"); 			

	fout << t << "\t" << q[0] << "\t" << q[1];
	fout << "\t" << q[2] << "\t" << q[3] << endl;	

	while( t < 100000000){

		rk4( q, t, tau, gravity, par);

		t = t + tau;			

		if(q[0] < 0){
			break;
		}

		fout << t << "\t" << q[0] << "\t" << q[1];
		fout << "\t" << q[2] << "\t" << q[3] << endl;	
	
	}

	fout.close();

//	cout << t << endl;

	FILE *gnuplot = popen("gnuplot", "w");		
	fprintf(gnuplot, "set out 'orbit_electron_rk.ps'\n ");
	fprintf(gnuplot, "set polar\n");
	fprintf(gnuplot, "set term post land\n ");
	fprintf(gnuplot, "set title 'The Collapse of Classical Electron in RK4'\n");
	fprintf(gnuplot, "set size square\n");
	fprintf(gnuplot, "set xlabel 'x (fm)'\n");
	fprintf(gnuplot, "set ylabel 'y (fm)'\n");
	fprintf(gnuplot, "plot 'orbit_electron_rk.dat' u 3:2 w l title'orbit'\n");

	pclose(gnuplot);					

	system( (char *) "ps2pdf orbit_electron_rk.ps");	

	delete[] q;
	delete[] q0;

	return 0;

}	