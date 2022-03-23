/*

Orbital motion
Velocity Verlet 

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


int main(int argc, char *argv[]){

	auto start = high_resolution_clock::now();

	double tau;
	double t = 0.0;
	double x0, y0, vx0, vy0;
	
	if(argc < 6){ 

		cout << "Usage: " << argv[0] << " x0 y0 vx0 vy0 tau" << endl;
		return -1;
	}
	
	else{

		x0 = atof(argv[1]);				
		y0 = atof(argv[2]);				
		vx0 = atof(argv[3]);
		vy0 = atof(argv[4]);
		tau = atof(argv[5]);
	}


	double r, rold;
	double x, xold, vx2, vx1, vxold, ax, axold;
	double y, yold, vy2, vy1, vyold, ay, ayold;

	xold = x0;
	x = x0;
	vxold = vx0;
	vx1 = vx0;
	vx2 = vx0;

	yold = y0;
	y = y0;
	vyold = vy0;
	vy1 = vy0;
	vy2 = vy0;

	r = sqrt(x*x + y*y);
	rold = sqrt(xold*xold + yold*yold);
	

	double energy;

	energy = 0.5 * (vx2*vx2 + vy2*vy2) - GM/r;


	ofstream fout ; 				
	fout.open("orbit_vl.dat"); 			

	fout << t << "\t" << x << "\t" << y;
	fout << "\t" << vx2 << "\t" << vy2 << "\t" << energy << endl;	


	int nmax = 10 / tau;
	int count = 0;


	while( count < nmax ){

		axold = -(GM/pow(rold,3)) * xold;
		vx1 = vxold + axold * tau / 2.0;
		x = xold + vx1 * tau;
		ax = -(GM/pow(r, 3)) * x;
		vx2 = vx1 + ax * tau / 2.0;

		ayold = -(GM/pow(rold,3)) * yold;
		vy1 = vyold + ayold * tau / 2.0;
		y = yold + vy1 * tau;
		ay = -(GM/pow(r,3)) * y;
		vy2 = vy1 + ay * tau / 2.0;

		xold = x;
		vxold = vx2;
		axold = ax;

		yold = y;
		vyold = vy2;
		ayold = ay;

		energy = 0.5 * (vx2*vx2 + vy2*vy2) - GM/r;

		t = t + tau;			


		// PRINTING OUT

		fout << t << "\t" << x << "\t" << y;
		fout << "\t" << vx2 << "\t" << vy2 << "\t" << energy << endl;	
	
		count++;
	}

	fout.close();

	FILE *gnuplot = popen("gnuplot", "w");			// opens a pipe
	fprintf(gnuplot, "set out 'orbit_vl.ps'\n ");		// Opens a window
	fprintf(gnuplot, "set term post land\n ");		// using a poster file in landscape
	fprintf(gnuplot, "set size square\n");
	fprintf(gnuplot, "set title 'Circular, Orbital Motion in Velocity Verlet'\n");
	fprintf(gnuplot, "set xlabel 'x (m)'\n");
	fprintf(gnuplot, "set ylabel 'y (m)'\n");
	fprintf(gnuplot, "plot 'orbit_vl.dat' u 2:3 w l title'orbit'\n");
	fprintf(gnuplot, "set title 'E/m vs time for Circular Orbit in Velocity Verlet'\n");
	fprintf(gnuplot, "set xlabel 't'\n");
	fprintf(gnuplot, "set ylabel 'E/m'\n");
	fprintf(gnuplot, "plot 'orbit_vl.dat' u 1:6 w l\n"); 
								// plot the data from the file using column 2 and 3 with lines
	pclose(gnuplot);					// closes the plotting

	system( (char *) "ps2pdf orbit_vl.ps");		// transforms the plotting into a pdf

	auto stop = high_resolution_clock::now();

	auto duration = duration_cast<microseconds>(stop - start);

	cout << "Time taken by Algorithm = " << duration.count() << " microseconds." << endl;

	return 0;

}	