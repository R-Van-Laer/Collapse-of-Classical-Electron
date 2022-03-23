/*

Collapse of Classical Electron
Velocity Verlet 

Raphael Van Laer

11-6-2021.

*/


#include <iostream>
#include <cmath>
#include <fstream>


using namespace std;


//Global Variables:
double PI = acos(-1.0);
double GM = 4 * PI * PI;
double c =3.0;
double a0 = 2.8;


int main(int argc, char *argv[]){

	double tau;
	double t = 0.0;
	double r0, theta0, vr0, vtheta0;
	
	if(argc < 4){ 

		cout << "Usage: " << argv[0] << " r0 theta0 tau" << endl;
		cout << "The radius will be in fm, meaning 1 = 10^-15 m." << endl;
		cout << "The angle will be in radians." << endl;
		cout << "The time will be in RVLs (invented second variant), meaning 1 = 10^-23 s." << endl;
		return -1;
	}
	
	else{

		r0 = atof(argv[1]);				
		theta0 = atof(argv[2]);				
		vr0 = 0.0;
		vtheta0 = sqrt(a0/r0) * c;
		tau = atof(argv[3]);
	}


	double r, rold, vr2, vr1, vrold, ar, arold;
	double theta, thetaold, vtheta2, vtheta1, vthetaold, atheta, athetaold;

	rold = r0;
	r = r0;
	vrold = vr0;
	vr1 = vr0;
	vr2 = vr0;

	thetaold = theta0;
	theta = theta0;
	vthetaold = vtheta0;
	vtheta1 = vtheta0;
	vtheta2 = vtheta0;


	ofstream fout ; 				
	fout.open("orbit_electron_vl.dat"); 			

	fout << t << "\t" << r << "\t" << theta;
	fout << "\t" << vr2 << "\t" << vtheta2 << endl;	


	while( t < 10000000 ){

		arold = -32.0 / 9.0 * (a0*a0*a0*a0)/(rold*rold*rold*rold*rold) * c*c - a0/(rold*rold) * c*c;
		vr1 = vrold + arold * tau / 2.0;
		r = rold + vr1 * tau;
		ar = -32.0 / 9.0 * (a0*a0*a0*a0)/(r*r*r*r*r) * c*c - a0/(r*r) * c*c;
		vr2 = vr1 + ar * tau / 2.0;

		athetaold = -8.0/3.0 * sqrt(a0/rold) * (a0*a0)/(rold*rold*rold) * c*c;
		vtheta1 = vthetaold + athetaold * tau / 2.0;
		theta = thetaold + vtheta1 * tau;
		atheta = -8.0/3.0 * sqrt(a0/r) * (a0*a0)/(r*r*r) * c*c;
		vtheta2 = vtheta1 + atheta * tau / 2.0;

		rold = r;
		vrold = vr2;
		arold = ar;

		thetaold = theta;
		vthetaold = vtheta2;
		athetaold = atheta;

		t = t + tau;			

		if(r < 0){
			break;
		}

		fout << t << "\t" << r << "\t" << theta;
		fout << "\t" << vr2 << "\t" << vtheta2 << endl;	
	
	}

	fout.close();

	FILE *gnuplot = popen("gnuplot", "w");		
	fprintf(gnuplot, "set out 'orbit_electron_vl.ps'\n ");
	fprintf(gnuplot, "set polar\n");
	fprintf(gnuplot, "set term post land\n ");
	fprintf(gnuplot, "set size square\n");
	fprintf(gnuplot, "set title 'The Collapse of CLassical Electron in Velocity Verlet'\n");
	fprintf(gnuplot, "set xlabel 'x (fm)'\n");
	fprintf(gnuplot, "set ylabel 'y (fm)'\n");
	fprintf(gnuplot, "plot 'orbit_electron_vl.dat' u 2:3 w l title'orbit'\n");

	pclose(gnuplot);					

	system( (char *) "ps2pdf orbit_electron_vl.ps");	

	return 0;

}	