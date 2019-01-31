#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "bckPotentialsAsymp.hpp"
#include "interp_1d.hpp"

using namespace Rcpp;
using namespace std;
using namespace boost::numeric::odeint ;

// The stiff algorith only accepts these data types
typedef boost::numeric::ublas::vector< long double > state_type;
typedef boost::numeric::ublas::matrix< long double > matrix_type;

struct VQCD
{
    long double x;
    long double tau0;
    long double W0;
    long double V0;
    long double lambda0;
    VQCD():x(1.0), tau0(1.0), W0(12.0/11), V0(12.0), lambda0(8.0 * M_PI * M_PI){}
    VQCD(long double xi, long double ti, long double w0, long double v0, long double l0):x(xi),tau0(ti), W0(w0), V0(v0), lambda0(l0){}
    void operator()(const state_type &X , state_type &dXdt , long double A)
    {
        long double e2A = exp(2.0 * A);
        long double kNr = k(X[1], x) ;
        long double dkNr = dk(X[1],x);
        long double vf = Vf(X[1], X[2], x, W0) ;
        long double dvfdl = dVfdl(X[1], X[2], x, W0) ;
        long double dlogvfdt = dLogVfdt(X[2], x, W0) ;
        long double dlogvfdl = dLogVfdl(X[1], x, W0) ;
        long double G = sqrt(1.0 + kNr * pow(X[4], 2.0) / ( e2A * pow(X[0], 2.0 ) ) ) ;
        // 1st Aeom in Mathematica Notebook 
        dXdt[0] = - X[0] + (4.0 / 9.0) * X[0] * pow(X[3], 2.0) / pow(X[1], 2.0) + x * e2A * pow(X[0], 3.0) * kNr * vf * pow(X[4], 2.0) * G / ( 6.0 * ( e2A * pow( X[0], 2.0 ) + kNr * pow(X[4], 2.0) ) ) ;
        dXdt[1] = X[3] ;
        dXdt[2] = X[4] ;
        // 2nd Aeom in Mathematica Notebook
        dXdt[3] = - (3.0 / 8) * e2A * pow( X[0] * X[1], 2.0 ) * dVg(X[1], V0, lambda0) + 9.0 * pow(X[1], 2.0) / X[3] \
        - 3 * e2A * pow(X[0] * X[1], 2.0) * Vg(X[1], V0, lambda0) / ( 4.0 * X[3] ) - 5.0 * X[3] + pow( X[3], 2.0) / X[1] \
        + ( 4.0 / 9) * pow(X[3], 3.0) / pow(X[1], 2.0) + 3 * e2A * x * pow( X[0] * X[1], 2.0) * vf / (4.0 * X[3] * G) \
        + x * kNr * vf * X[3] * pow( X[4], 2.0) / ( 6.0 * G) + 3.0 * x * vf * pow(X[1] * X[4], 2.0) * dkNr / ( 16.0 * G) \
        + 3 * e2A * x * pow(X[0] * X[1], 2.0) * dvfdl / ( 8.0 * G) + 3.0 * x * kNr * pow(X[1] * X[4], 2.0) * dvfdl / ( 8.0 * G) ;
        // 3rd Aeom in Mathematica Notebook
        dXdt[4] = - 4.0 * X[4] + ( 4.0 / 9 ) * pow( X[3] / X[1] , 2.0) * X[4] - 4.0 * kNr * pow(X[4], 3.0) / ( e2A * X[0] * X[0]) \
        + e2A * x * pow( X[0], 2.0 ) * kNr * vf * pow( X[4], 3.0) * G / ( 6 * ( e2A * pow(X[0],2.0) + kNr * pow(X[4], 2.0) ) )\
        - X[3] * X[4] * dkNr / kNr - X[3] * pow(X[4], 3.0) * dkNr / (2.0 * e2A * pow( X[0], 2.0) ) \
        + e2A * pow(X[0], 2.0) * dlogvfdt / kNr + pow(X[4], 2.0) * dlogvfdt - X[3] * X[4] * dlogvfdl \
        - kNr * X[3] * pow(X[4], 3.0) * dlogvfdl / ( e2A * pow(X[0], 2.0) ) ;
    }
};

// [[Rcpp::export]]
List solveHVQCD(long double xi = 1, long double ti = 1, long double W0 = 12/11.0, long double V0 = 12, long double lambda0 = 8 * M_PI * M_PI)
{
    // Computes dr/dA, tau, lambda, dtau/dA and dlambda/dA given x and tau0
    // x - long double. Physically it means x = N_f / N_c when N_f, N_c -> inf but with fixed coefficient
    // tau0 - > parameter related with the exponential behaviour of the tachyon in the IR
    // Returns a list with quantitites that depend on A, mq and zIR
    
    // Create a vectors containing the values of the fields
    vector< long double > Z, AA, dZ, L, T;
    // Boundary conditions in the IR
    long double zIR = log(70.0 / ti) / CI(xi, W0, V0, lambda0) ;
    long double aIR = AIR(zIR, V0, lambda0) ;
    AA.push_back(aIR) ;
    long double lambdair = lambdaIR(zIR, lambda0) ;
    L.push_back(lambdair);
    long double tauir = 70.0 ;
    T.push_back(tauir) ;
    // dA/dz at zIR
    long double daIR =  173.0 / ( 1728.0 * pow(zIR, 3.0)) + 0.5 /zIR - 2.0 * zIR  ;
    // dtau/dz at zIR
    long double dtauir = CI(xi, W0, V0, lambda0) * tauir ;
    // dlambda/dz at zIR. Expression got from eom2 in Mathematica notebook.
    long double dlambdair = sqrt(1.5) * lambdair * sqrt( 6 * pow(daIR, 2.0) + exp( 2.0 * aIR) * xi * Vf(lambdair,tauir, xi, W0) \
         / (2.0 * sqrt(1+ dtauir * dtauir * k(lambdair, xi) / exp(2 * aIR ) ) ) - 0.5 * exp( 2.0 * aIR) * Vg(lambdair, V0, lambda0) );
    long double dzIR = 1.0 / daIR ;
    dZ.push_back(dzIR) ;
    dlambdair = dlambdair / daIR;
    dtauir = dtauir / daIR;
    // Define the type of the state. We have X = {dz, lambda, tau, dlambda, dtau}
    state_type X (5);
    X <<= dzIR, lambdair, tauir, dlambdair, dtauir;
    // Now compute the starting value of dX
    cout << "Solving HVQCD for x = " << xi << ", tau0 = " << ti << ", W0 = " << W0 << ", V0 = " << V0 << ", lambda0 = " << lambda0 << endl;
    long double Amax = 100.0;
    long double h = 0.1 ;
    dense_output_runge_kutta< controlled_runge_kutta< runge_kutta_dopri5< state_type > > > stepper;
    stepper.initialize( X , aIR , h );
    long double A = aIR;
    while ( A < Amax )
    {
        stepper.do_step( VQCD(xi, ti, W0, V0, lambda0) ) ; 
        X = stepper.current_state();
        A = stepper.current_time();
        //cout << A << '\t' <<  X(0) << '\t' <<  X(1) << '\t' <<  X(2) << endl;
        AA.push_back(A) ;
        dZ.push_back(X(0)) ;
        L.push_back(X(1)) ;
        T.push_back(X(2)) ;
    }
    double mq = X(4) / X(0) ;
    Spline_Interp<long double> dzfun = Spline_Interp<long double>(AA, dZ);
    int n = AA.size() ;
    long double zmin = zIR + dzfun.integrate(AA[ n - 1]) ;
    for(int i = 0; i < n; i++)
    {
        Z.push_back( zIR + dzfun.integrate(AA[i]) - zmin ) ;
    }
    // We need to reverse the lists because later
    // we will use them to compute the spectra of vector mesons.
    reverse(AA.begin(),AA.end()) ;
    reverse(Z.begin(),Z.end()) ;
    reverse(L.begin(),L.end()) ;
    reverse(T.begin(),T.end()) ;
    // Return A, Z, L(A), T(A), zIR and mq
    return List::create(Named("z") = Z, Named("A") = AA, Named("lambda") = L, Named("tau") = T, Named("zIR") = zIR - zmin, Named("mq") = mq );
} ;