// *********************************************************
// *** Simulations for Kuramoto Oscillators ****************
// *** Two pop set up. Look at Panaggio et. al 2016 ********
// *** Using Euler's Method   ******************************
// ***  g++ -O3 -o example chimera_in_small_pop.cpp mt.cpp *
// *********************************************************

#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>        // std::ostringstream
#include <iomanip>        //setw, setfill
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

/* Coupling term: non-local coupling like Zakharova et. al*/
double non_local_coupling_zak (int J, int P, VecDoub theta, double act, int j)
{
    double coup = 0.;
    
     for (int z = j - P; z <= j + P; z++)
        coup += sin(theta[z]-(theta[j]+act));
     return coup;
    
}

// ********************************************************
// Coupling terms: between nodes and between populations
double pop_coupling_theta(int J, VecDoub& theta, VecDoub& phi, int j, double beta, double A)
{
    double coup = 0;
    for (int z = 0; z < J; z++)
    {
        coup += ((1.+ A)/(2*J))*cos(theta[j]-theta[z] - beta) + ((1.- A)/(2*J))*cos(theta[j]-phi[z] - beta);
    }
    return coup;
}

double pop_coupling_phi(int J, VecDoub& theta, VecDoub& phi, int j, double beta, double A)
{
    double coup = 0;
    for (int z = 0; z < J; z++)
    {
        coup += ((1.+ A)/(2*J))*cos(phi[j]-phi[z] - beta) + ((1.- A)/(2*J))*cos(phi[j]-theta[z] - beta);
    }
    return coup;
}

void KRM(VecDoub om, int j, double A, int J, VecDoub& theta, VecDoub& ftheta, VecDoub& phi, VecDoub& fphi, double beta, double TC)
{
    ftheta[j] = om[j] - pop_coupling_theta(J, theta, phi, j, beta, A);
    ftheta[j] *= TC;
    fphi[j] = om[j] - pop_coupling_phi(J, theta, phi, j, beta, A);
    fphi[j] *= TC;
}

void E(VecDoub& theta, VecDoub& ftheta, VecDoub& phi, VecDoub& fphi, int j, double dt)
{
    theta[j] += dt*ftheta[j];
    phi[j] += dt*fphi[j];
}



// ******************************************
int main (int argc, char **argv)
{
    time_t start;
    time(&start);


    /* CHOOSE BETWEEN THE TWO SYSTEMS J = 3 OR J = 25*/
    int J = 3;
    double num = 0; //if J = 3, num = 7
    
    /*
    int J = 25;
    double num = 2; //if J = 25, num = 2*/
   

    /* Definitions */
    double A = 0.1;
    double beta = 0.025;
    double TC = 1; //time constant. It is used in order to slower it down or speed it up. 
    int tfinal = 5e4;
    double dt = 0.01;
    int N = int(tfinal/dt);
    VecDoub theta(J);
    VecDoub ftheta(J);
    VecDoub phi(J);
    VecDoub fphi(J);
    double t;
    VecDoub theta0(J);
    VecDoub phi0(J);
    double omega = 2.6;
    VecDoub om(J, omega);
    int l;
    int step = (1000/dt);

    
    double seed = strtod(argv[1], 0);
    MersenneTwister mt(seed);
    double epsilon = strtod(argv[2], 0);
    MersenneTwister mt_2(seed+23);
    double D, D_2;
    
    /* Theta and Phi intial conditions */
    string nestring = nameconv(num)+string("_initial_conditions.txt");
    char *nechar= new char[nestring.length()+1];
    strcpy(nechar,nestring.c_str());
    ifstream fthe;
    fthe.open (nechar, ios::in);
    
    int doubleJ = 2*J;
    VecDoub chimera(doubleJ);
    for (int col = 0; col < doubleJ; col++)
        fthe >> chimera[col];

    for (int col = 0; col < J; col++)
        theta0[col] = chimera[col];
    
    for (int col = J; col < doubleJ; col++)
        phi0[col-J] = chimera[col];
        
    theta = theta0;
    phi = phi0;
    
    /* Perturbation */
    int length_pert = 44;
    VecDoub t_pert(length_pert);
    
    string nelstring = string("random_kick_time_A_0.1_N_3_beta_0.025_seed_4_epsilon_0.3_omega_2.6_skip_nt_500000.txt");
    char *nelchar= new char[nelstring.length()+1];
    strcpy(nelchar,nelstring.c_str());
    ifstream fthel;
    fthel.open (nelchar, ios::in);
    
    for (int col = 0; col < length_pert; col++)
        fthel >> t_pert[col];

   
    
    int transients = int(5000/dt);
    
    string nostring = string("theta_A_")+nameconv(A)+string("_N_")+nameconv(J)+("_beta_0.025_seed_")+nameconv(seed)+string("_epsilon_")+nameconv(epsilon)+string("_omega_")+nameconv(omega)+string("_skip_nt_")+nameconv(transients)+string("_1_kick.txt");
    char *nochar= new char[nostring.length()+1];
    strcpy(nochar,nostring.c_str());
    ofstream fitxer1;
    fitxer1.open (nochar, ios::trunc);
    
    string cstring = string("phi_A_")+nameconv(A)+string("_N_")+nameconv(J)+("_beta_0.025_seed_")+nameconv(seed)+string("_epsilon_")+nameconv(epsilon)+string("_omega_")+nameconv(omega)+string("_skip_nt_")+nameconv(transients)+string("_1_kick.txt");
    char *cchar= new char[cstring.length()+1];
    strcpy(cchar,cstring.c_str());
    ofstream fitxer2;
    fitxer2.open (cchar, ios::trunc);

    string costring = string("random_kick_time_A_")+nameconv(A)+string("_N_")+nameconv(J)+("_beta_0.025_seed_")+nameconv(seed)+string("_epsilon_")+nameconv(epsilon)+string("_omega_")+nameconv(omega)+string("_skip_nt_")+nameconv(transients)+string("_1_kick.txt");
    char *cochar= new char[costring.length()+1];
    strcpy(cochar,costring.c_str());
    ofstream fitxer3;
    fitxer3.open (cochar, ios::trunc);

    int p = 0;
    int ij = 0;
    for (int i = 0; i < N; i++)
    {
        t = i*dt;
        if (i%10 == 0)
        {
            fitxer1 << t << "\t";
            fitxer2 << t << "\t";
        }

        for (int j = 0; j < J; j++)
        {

            KRM(om, j, A, J, theta, ftheta, phi, fphi, beta, TC);
            E(theta, ftheta, phi, fphi, j, dt);
            
            if (abs(theta[j]) > M_PI)
            {
                if (theta[j] > 0)
                    theta[j] -= 2*M_PI;

                else if (theta[j] < 0)
                    theta[j] += 2*M_PI;
            }
            
            if (abs(phi[j]) > M_PI)
            {
                if (phi[j] > 0)
                    phi[j] -= 2*M_PI;

                else if (phi[j] < 0)
                    phi[j] += 2*M_PI;
            }
        
            
            if (t==t_perturb[ij])
            {

            if (i>transients && i%(step+l)==0)
            {
                if(j == 0)
                {
                    if (p%2 == 0)
                        noise_polar(mt_2,&D,&D_2);
                    else
                        D = D_2;
                    p +=1;
                    fitxer3 << t << "\t";
                }

                theta[j] += theta[j] + epsilon*D;
                phi[j] += phi[j] + epsilon*D;

            }
            
            if (i%10 == 0)
            {
                fitxer1 << theta[j] << "\t";
                fitxer2 << phi[j] << "\t";
            }
            
        }


        
        
        if (i%10 == 0)
        {
            fitxer1 << endl;
            fitxer2 << endl;
        }
        
    }


    time_t fin;
    time(&fin);
    cout << fin-start << endl;

    
}
