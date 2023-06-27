// *********************************************************
// *** Simulations for Kuramoto Oscillators ****************
// *** Chimeras set up. Look at Laing et. al 2009 **********
// ***  g++ -O3 -o ring chimera_on_a_ring.cpp mt.cpp *
// *********************************************************


#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>		// std::ostringstream
#include <iomanip>		//setw, setfill
#include <stdlib.h>     // strtod
#include "stdio.h"
#include "mt.h"
#include "string.h"
#include "helfer.h"      //needed for matrix&vector class
#include <time.h>
#include <vector>
using namespace std;

//***********************************************************
template <typename T>  string nameconv(T param)
{
    ostringstream convert;   // stream used for the conversion
    convert << param;      // insert the textual representation of'Number' in the characters in the stream
	string paramstr = convert.str(); // set 'paramstr' to the contents of the stream
    convert.str("");
    return paramstr;
}
//***********************************************************
// Polar-Method should be slightly faster:
// It generates a pair of normal random variables
// This has the additional advantage that you use both random numbers
template <class T>
void noise_polar (T& mt, double *x1, double *x2)
{
   double u, v, q, p;

   do 
   {
      u = 2.0 * mt.random() - 1;
      v = 2.0 * mt.random() - 1;
      q  = u * u + v * v;
   } 
   
   while (q == 0.0 || q >= 1.0);

   p = sqrt(-2 * log(q) / q);
   *x1 = u * p;
   *x2 = v * p;
}

// ********************************************************
// Coupling term: all to all coupled
double all_to_all_coupling(int J, VecDoub& theta, int j, double K)
{
	double coup = 0;
    for (int z = 0; z < J; z++)
    {
		coup += sin(theta[z]-(theta[j]));
	}
    return (K/J)*coup;
}

// ********************************************************
// Coupling term: non-local coupling like Laing et. al 2009
double non_local_coupling(int J, VecDoub& theta, int j, double beta, double K)
{
	double coup = 0;
    for (int z = 0; z < J; z++)
    {
		coup += (1.+K*cos(2*M_PI*abs(j-z)/J))*cos(theta[j]-theta[z] - beta);
	}
    return coup/J;
}

/* Coupling term: non-local coupling like Zakharova et. al
 * double non_local_coupling (int J, int P, VecDoub theta, double act, int j) 
{
	double coup = 0.;
	
	 for (int z = j - P; z <= j + P; z++)
		coup += sin(theta[z]-(theta[j]+act));
     return coup;
	
}*/

void KRM(double TC, VecDoub om, int j, double K, int J, VecDoub& theta, VecDoub& ftheta, double beta)
{
    ftheta[j] = om[j] - non_local_coupling(J, theta, j, beta, K);
    ftheta[j] *= TC;
}

void E(VecDoub& theta, VecDoub& ftheta, int j, double dt)
{
    theta[j] += dt*ftheta[j]; 
}



// ********************************************************
int main (int argc, char **argv)
{
	time_t start;
	time(&start);

	/* Definitions */
    int J = 500;
    double A = 0.95;
    double beta = 0.2;
    double omega = 2.8; // fixed omega
    double TC = 1; //fixed TC. For user input use: strtod(argv[1], 0);
    int tfinal =5e4;
	double dt = 0.01;
	int N = int(tfinal/dt);		
	VecDoub theta(J);
	VecDoub ftheta(J);
	double t, c, c_2, p, p_2, K;
    VecDoub marcador(J);
    VecDoub rotations(J);
    VecDoub om(J,omega);
    VecDoub theta0(J);
    
    /* Theta intial conditions */
    string nestring = string("initial_theta_N_500_chimera.txt");
	char *nechar= new char[nestring.length()+1];
	strcpy(nechar,nestring.c_str());
	ifstream fthe;
    fthe.open (nechar, ios::in);
	
    VecDoub chimera(N);
    for (int col = 0; col < N; col++)
        fthe >> chimera[col];
    
    theta = chimera;
    
    string nostring = string("theta_A_")+nameconv(A)+string("_N_")+nameconv(J)+("_beta_0.2_TC_")+nameconv(TC)+string("_ic_chimera_om_")+nameconv(omega)+string("_tfinal_")+nameconv(tfinal)+string(".txt");
    char *nochar= new char[nostring.length()+1];
    strcpy(nochar,nostring.c_str());
    ofstream fitxer1;
    fitxer1.open (nochar, ios::trunc);
    
    for (int i = 0; i < N; i++)
    {
        t = i*dt;
        if (i%10 == 0)
        {

            fitxer1 << t << "\t";
        }
        for (int j = 0; j < J; j++)
        {
            KRM(TC,om, j, A, J, theta, ftheta, beta);
            E(theta, ftheta, j, dt);
            if (abs(theta[j]) > M_PI)
            {
                if (theta[j] > 0)
                    theta[j] -= 2*M_PI;

                else if (theta[j] < 0)
                    theta[j] += 2*M_PI;
            }
            

                
            if (i%10 == 0)
            {

                fitxer1 << theta[j] << "\t";
            }
        }
        if (i%10 == 0)
        {

            fitxer1 << endl;
        }



    


}
	
	time_t fin;
	time(&fin);
	cout << fin-start << endl; 
}
