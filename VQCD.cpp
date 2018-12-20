#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <math.h>
#include <boost/array.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include <boost/numeric/odeint/integrate/max_step_checker.hpp>

using namespace std;
using namespace boost::numeric::odeint ;
namespace phoenix = boost::phoenix;

// The stiff algorith only accepts these data types
typedef boost::numeric::ublas::vector< double > state_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;

struct output_observer
{
    vector< vector< double > > &data;
    output_observer( vector< vector < double > > &fieldsData ) : data( fieldsData ) { }
    void operator()( const state_type &X , double A ) const
    {
        // Add A to vector of As
        data[0].push_back(A);
        // Add dz/dA to the vector of dzs
        data[1].push_back(X[0]);
        // Add lambda to the vector of lambdas
        data[2].push_back(X[1]);
        // Add tau to the vector of taus
        data[3].push_back(X[2]);
        // Add dlambda to the vector of dlambdas
        data[4].push_back(X[3]);
        // Add dtau to the vector of dtaus
        data[5].push_back(X[4]);
        //cout << A << '\t' << X[0] << '\t' <<  X[1] << '\t'<< X[2] << '\t' << X[3] << '\t' << X[4] << endl;
    }
};

// Definition of constants need to solve ODE
const double V0 = 12.0 ;
const double V1 = 11.0 / ( 27 * M_PI * M_PI);
const double V2 = 4619.0 / ( 46656.0 * pow(M_PI,4.0) );
const double W0 = 12.0/11;
const double l0 = 8 * pow(M_PI, 2.0);
// Definition of functions needed to solve ODE
double W1 (double x) 
{
    return (24 + (11 - 2 * x) * W0) / ( 27.0 * pow(M_PI,2.0) * W0 ) ; 
}
double W2 (double x) 
{
    return (24 * (857 - 46 * x) + (4619 - 1714 * x + 92 * pow(x, 2.0)) * W0) / (46656 * pow(M_PI,4.0) * W0) ;
}
double a0 (double x) 
{
    return (12 - x * W0) / 8.0 ;
}
double a1 (double x) 
{
    return (115 - 16 * x) / (216.0 * pow(M_PI, 2.0)) ;
}
double Vg (double l) 
{
    return V0 * (1+ V1 * l + V2 * pow(l,2.0) * sqrt(1 + log(1.0 + l / l0))/pow(1.0 + l / l0,2.0/3)) ;
}
double dVg (double l)
{
    double ans = V1 + V2 * pow(l, 2.0) / ( 2.0 * pow( 1 + l / l0, 5.0/3) * l0 * sqrt( 1 + log( 1 + l / l0 )));
    ans += 2.0 * V2 * l * sqrt( 1 + log( 1 + l / l0 )) / pow( 1 + l / l0, 2.0/3);
    ans += - 2.0 * V2 * pow(l, 2.0) * sqrt( 1 + log( 1 + l / l0 )) / ( 3.0 * pow( 1 + l / l0, 5.0/3) * l0 );
    return V0 * ans;
}
double d2Vg( double l)
{
    double ans = V0 * V2 ;
    ans = ans * (37 * l * l + 120.0 * l * l0 + 72 * l0 *l0 + 2.0 * log(1.0 + l / l0) * (31 * l * l + 84.0 * l * l0  + 72.0 * l0 * l0 + 4.0 * (2 * l * l + 6 * l * l0 + 9 * l0 * l0) * log(1.0 + l / l0)))  ;
    ans = ans / (36.0 * l0 * l0 * pow(1.0  + l / l0, 8.0/3.0) * pow(1 + log(1.0 + l / l0),1.5)) ;
    return ans;
}
double Vf0 (double l, double x) 
{
    return W0 * ( 1 + W1(x) * l + W2(x) * pow(l, 2.0)) ;
}
double k (double l, double x) 
{
    return 1.0 / pow(1 + 3 * a1(x) * l / 4.0, 4.0 / 3) ;
}
double dk(double l, double x)
{
    return - a1(x) / pow( 1 + 3 * l * a1(x) / 4.0, 7.0 /3 );
}
double d2k(double l, double x)
{
    return 7 * a1(x) * a1(x) / (4.0 * pow(1 + 0.75 * l * a1(x),10.0/3)) ;
}
double Vf (double l, double t, double x) 
{
    return exp( log( Vf0(l,x) ) - a0(x) * pow(t, 2.0) );
}
double dVfdl(double l, double t, double x)
{
    return exp( log( W0 * ( W1(x) + 2.0 * l * W2(x) )) - a0(x) * t * t);
}
double d2Vfdl(double l, double t, double x)
{
    return 2.0 * W0 * W2(x) * exp(- a0(x) * t * t);
}
double d2Vfdl_Vf(double l, double t, double x)
{
    return 2.0 * W0 * W2(x) / Vf0(l,x) ;
}
double d2Vfdldt(double l, double t, double x)
{
    return W0 * (W1(x) + 2 * l * W2(x)) * (-2.0 * a0(x) * t) * exp(-a0(x) * t * t);
}
double d2Vfdldt_Vf(double l, double t, double x)
{
    return W0 * (W1(x) + 2 * l * W2(x)) * (-2.0 * a0(x) * t) / Vf0(l,x) ;
}
double dVfdt(double l, double t, double x)
{
    return  (- 2 * a0(x) * t) * exp( log( Vf0(l,x) ) - a0(x) * pow(t, 2.0) ) ;
}
double d2Vfdt(double l, double t, double x)
{
    return 2 * W0 * a0(x) * exp(-a0(x)*t*t) * (-1 + 2 * t * t * a0(x)) * (1 + l * W1(x) + l * l * W2(x) );
}
double d2Vfdt_Vf(double l, double t, double x)
{
    return 2 * a0(x) * (-1 + 2 * t * t * a0(x)) ;
}
double dLogVfdt (double t, double x)
{
    return - 2.0 * t * a0(x) ;
}
double dLogVfdl(double l, double x)
{
    return ( W1(x) + 2.0 * l * W2(x) ) / ( 1.0 + l * W1(x) + W2(x) * l * l) ;
}

// Definition of the IR assymptotics
double AIR( double z)
{
    return 13.0 / 8 + log(27 * pow(6, 0.25) / sqrt(4619) ) - (173.0 / 3456 ) / pow(z, 2.0) - pow(z, 2.0) + 0.5 * log(z) ;
}

double lambdaIR (double z)
{
    return exp( log(8 * pow(M_PI, 2.0)) - (39.0 / 16) - ( 151.0 / 2304) / pow(z, 2.0) + 1.5 * pow(z, 2.0)) ;
}

double CI ( double x)
{
    return 81 * pow(3.0, 5.0 / 6) * pow( 115 - 16 * x, 4.0/3) * (11.0 - x) / ( 812944 * pow(2.0, 1.0 / 6)) ;
}

double tauIR ( double z, double t0, double x)
{
    return t0 * exp( CI(x) * z) ;
}

struct VQCD
{
    double x;
    double tau0;
    VQCD():x(1.0), tau0(1.0){}
    VQCD(double xi, double ti):x(xi),tau0(ti){}
    void operator()(const state_type &X , state_type &dXdt , double A)
    {
        double e2A = exp(2.0 * A);
        double kNr = k(X[1], x) ;
        double dkNr = dk(X[1],x);
        double vf = Vf(X[1], X[2], x) ;
        double dvfdl = dVfdl(X[1], X[2], x) ;
        double dlogvfdt = dLogVfdt(X[2], x) ;
        double dlogvfdl = dLogVfdl(X[1], x) ;
        double G = sqrt(1.0 + kNr * pow(X[4], 2.0) / ( e2A * pow(X[0], 2.0 ) ) ) ;
        // 1st Aeom in Mathematica Notebook 
        dXdt[0] = - X[0] + (4.0 / 9.0) * X[0] * pow(X[3], 2.0) / pow(X[1], 2.0) + x * e2A * pow(X[0], 3.0) * kNr * vf * pow(X[4], 2.0) * G / ( 6.0 * ( e2A * pow( X[0], 2.0 ) + kNr * pow(X[4], 2.0) ) ) ;
        dXdt[1] = X[3] ;
        dXdt[2] = X[4] ;
        // 2nd Aeom in Mathematica Notebook
        dXdt[3] = - (3.0 / 8) * e2A * pow( X[0] * X[1], 2.0 ) * dVg(X[1]) + 9.0 * pow(X[1], 2.0) / X[3] \
        - 3 * e2A * pow(X[0] * X[1], 2.0) * Vg(X[1]) / ( 4.0 * X[3] ) - 5.0 * X[3] + pow( X[3], 2.0) / X[1] \
        + ( 4.0 / 9) * pow(X[3], 3.0) / pow(X[1], 2.0) + 3 * e2A * x * pow( X[0] * X[1], 2.0) * vf / (4.0 * X[3] * G) \
        + x * kNr * vf * X[3] * pow( X[4], 2.0) / ( 6.0 * G) + 3.0 * x * vf * pow(X[1] * X[4], 2.0) * dkNr / ( 16.0 * G) \
        + 3 * e2A * x * pow(X[0] * X[1], 2.0) * dvfdl / ( 8.0 * G) + 3.0 * x * kNr * pow(X[1] * X[4], 2.0) * dvfdl / ( 8.0 * G) ;
        // 3rd Aeom in Mathematica Notebook
        dXdt[4] = - 4.0 * X[4] + ( 4.0 / 9 ) * pow( X[1] / X[3] , 2.0) * X[4] - 4.0 * kNr * pow(X[4], 3.0) / ( e2A * X[0] * X[0]) \
        + e2A * x * pow( X[0], 2.0 ) * kNr * vf * pow( X[4], 3.0) * G / ( 6 * ( e2A * pow(X[0],2.0) + kNr * pow(X[4], 2.0) ) )\
        - X[3] * X[4] * dkNr / kNr - X[3] * pow(X[4], 3.0) * dkNr / (2.0 * e2A * pow( X[0], 2.0) ) \
        + e2A * pow(X[0], 2.0) * dlogvfdt / kNr + pow(X[4], 2.0) * dlogvfdt - X[3] * X[4] * dlogvfdl \
        - kNr * X[3] * pow(X[4], 3.0) * dlogvfdl / ( e2A * pow(X[0], 2.0) ) ;
    }
};

struct jacobian
{
    double x ;
    jacobian( double xi):x(xi){}
    void operator()( const state_type & X, matrix_type &J , const double &A , state_type &dfdt )
    {
        double X1 = X[0] ;
        double X2 = X[1] ;
        double X3 = X[2] ;
        double X4 = X[3] ;
        double X5 = X[4] ;
        double e2A = exp(2.0 * A);
        double e4A = exp(4.0 * A);
        double kNr = k(X2, x) ;
        double dkNr = dk(X2,x);
        double d2kNr = d2k(X2,x) ;
        double vg = Vg(X2) ;
        double dvg = dVg(X2);
        double d2vg = d2Vg(X2) ;
        double vf = Vf(X2, X3, x) ;
        double dvfdl = dVfdl(X2, X3, x) ;
        double d2vfdl = d2Vfdl(X2, X3, x) ;
        double d2vfdl_vf = d2Vfdl_Vf(X2, X3, x) ;
        double d2vfdldt = d2Vfdldt(X2, X3, x) ;
        double d2vfdldt_vf = d2Vfdldt_Vf(X2, X3, x) ;
        double dvfdt = dVfdt(X2, X3, x) ;
        double d2vfdt = d2Vfdt(X2, X3, x);
        double d2vfdt_vf = d2Vfdt_Vf(X2,X3, x) ;
        double dlogvfdt = dLogVfdt(X3, x) ;
        double dlogvfdl = dLogVfdl(X2, x) ;
        double G = sqrt(1.0 + kNr * pow(X5, 2.0) / ( e2A * pow(X1, 2.0 ) ) ) ;
        double denom1 = e2A * pow(X1,2.0) + pow(X5,2.0) * kNr;
        J( 0 , 0 ) = - e2A * pow(X1,2.0) / denom1 + (4.0 / 9) * e2A * pow(X1 * X4,2.0) / (pow(X2,2.0) * denom1) \
        - pow(X5, 2.0) * kNr / denom1 + (4.0/9) * pow(X4*X5,2.0) * kNr / (pow(X2,2.0) * denom1) \
        + e2A * x * pow(X1*X5,2.0) * kNr * vf / (6.0 * denom1 * G) + x * pow(X5,4.0) * pow(kNr,2.0) * vf / (3.0 * denom1 * G);
        J( 0 , 1 ) = - (8.0 / 9) * e2A * pow(X1,3.0) * pow(X4,2.0) / (pow(X2,3.0) * denom1) \
        - (8.0 / 9) * X1 * pow(X4 * X5, 2.0) * kNr / ( pow(X2,3.0) * denom1) \
        + e2A * x * pow(X1,3.0) * pow(X5,2.0) * vf * dkNr / (6 * denom1 * G) \
        + x * X1 * pow(X5,4.0) * kNr * vf * dkNr / (12.0 * denom1 * G) \
        + e2A * x * pow(X1,3.0) * pow(X5,2.0) * kNr * dvfdl / (6.0 * denom1 * G) \
        + x * X1 * pow(X5, 4.0) * kNr * kNr * dvfdl / (6.0 * denom1 * G) ;
        J( 0 , 2 ) = e2A * x * pow(X1,3.0) * pow(X5,2.0) * kNr * G * dvfdt / (6.0 * denom1);
        J( 0 , 3 ) = (8.0 / 9) * X1 * X4 / pow(X2,2.0);
        J( 0 , 4 ) = e2A * x * pow(X1,3.0) * X5 * kNr * vf / (3.0 * denom1 * G) \
        + x * X1 * pow(X5,3.0) * pow(kNr,2.0) * vf / (6.0 * denom1 * G);
        J( 1 , 0 ) = 0.0;
        J( 1 , 1 ) = 0.0;
        J( 1 , 2 ) = 0.0;
        J( 1 , 3 ) = 1.0;
        J( 1 , 4 ) = 0.0;
        J( 2 , 0 ) = 0.0;
        J( 2 , 1 ) = 0.0;
        J( 2 , 2 ) = 0.0;
        J( 2 , 3 ) = 0.0;
        J( 2 , 4 ) = 1.0;
        J( 3 , 0 ) = 1.5 * e4A * x * pow(X1,3.0) * pow(X2,2.0) * vf / (X4 * denom1 * G) \
        + (9.0 / 4) * e2A * x * X1 * pow(X2 * X5,2.0) * kNr * vf / (X4 * denom1 * G) \
        + x * X4 * pow(X5, 4.0) * pow(kNr,2.0) * vf / (6 * X1 * denom1 * G) \
        - 1.5 * e4A * pow(X1,3.0) * pow(X2,2.0) * vg / (X4 * denom1) \
        - 1.5 * e2A * X1 * pow(X2*X5,2.0) * kNr * vg / (X4 * denom1) \
        - 0.75 * e4A * pow(X1,3.0) * pow(X2,2.0) * dvg / denom1 - 0.75 * e2A * X1 * pow(X2*X5,2.0) * kNr * dvg / denom1 \
        + (3.0 / 16) * x * pow(X2,2.0) * pow(X5, 4.0) * kNr * vf * dkNr / ( X1 * denom1 * G) \
        + 0.75 * e4A * x * pow(X1,3.0) * pow(X2, 2.0) * dvfdl / (denom1 * G) \
        + (9.0 / 8) * e2A * x * X1 * pow(X2 * X5,2.0) * kNr * dvfdl / (denom1 * G) \
        + (3.0 / 8) * x * pow(X2 * kNr,2.0) * pow(X5, 4.0) * dvfdl / (X1 * denom1 * G);
        J( 3 , 1 ) = 18.0 * e4A * pow(X1,4.0) * X2 / ( X4 * pow(denom1,2.0)) - e4A * pow(X1,4.0) * pow(X4/X2,2.0) / (pow(denom1,2.0)) \
        - (8.0 / 9) * e4A * pow(X1, 4.0) * pow(X4 / X2, 3.0) / pow(denom1,2.0) \
        + 36 * e2A * pow(X1 * X5, 2.0) * X2 * kNr / (X4 * pow(denom1,2.0)) \
        - 2 * e2A * pow(X1 * X4 * X5,2.0) * kNr / pow(X2 * denom1,2.0) \
        - (16.0/9) * e2A * pow(X1 * X5,2.0) * pow(X4 / X2,3.0) * kNr / pow(denom1,2.0) \
        + 18 * X2 * pow(X5,4.0) * pow(kNr,2.0) / (X4 * pow(denom1,2.0)) - pow(X4 * kNr / (X2 * denom1),2.0) * pow(X5, 4.0) \
        - (8.0/9) * pow(X4/X2,3.0) * pow(X5,4.0) * pow(kNr / denom1,2.0) \
        + 1.5 * e4A * e2A * x * pow(X1,6.0) * X2 * G * vf / (X4 * pow(denom1,2.0)) \
        + 1.5 * e4A * x * pow(X1,4.0) * X2 * pow(X5,2.0) * kNr * G * vf / (X4 * pow(denom1, 2.0)) \
        - 1.5 * e2A * e4A * pow(X1,6.0) * X2 * vg / (X4 * pow(denom1, 2.0)) \
        - 3 * e4A * pow(X1,4.0) * X2 * pow(X5,2.0) * kNr * vg / (X4 * pow(denom1,2.0)) \
        -1.5 * e2A * pow(X1 * kNr, 2.0) * X2 * pow(X5,4.0) * vg / (X4 * pow(denom1,2.0)) \
        - 0.75 * e2A * e4A * pow(X1,6.0) * X2 * dvg / pow(denom1,2.0) \
        - 0.75 * e2A * e4A * pow(X1, 6.0) * pow(X2,2.0) * dvg / (X4 * pow(denom1,2.0)) \
        - 1.5 * e4A * pow(X1,4.0) * X2 * pow(X5, 2.0) * kNr * dvg / pow(denom1, 2.0) \
        - 1.5 * e4A * pow(X1,4.0) * pow(X2 * X5, 2.0) * kNr * dvg / (X4 * pow(denom1,2.0)) \
        - 0.75 * e2A * pow(X1 * kNr,2.0) * X2 * pow(X5,4.0) * dvg / pow(denom1,2.0) \
        - 0.75 * e2A * pow(X1 * X2 * kNr,2.0) * pow(X5, 4.0) * dvg / (X4 * pow(denom1,2.0)) \
        - (3.0 / 8 ) * e2A * e4A * pow(X1, 6.0) * pow(X2,2.0) * d2vg / pow(denom1,2.0) \
        - 0.75 * e4A * pow(X1,4.0) * pow(X2 * X5,2.0) * kNr * d2vg / pow(denom1,2.0) \
        - (3.0 / 8) * e2A * pow(X1 * X2 * kNr, 2.0) * pow(X5, 4.0) * d2vg / pow(denom1,2.0) \
        + (3.0 / 8) * e4A * x * pow(X1,4.0) * X2 * pow(X5,2.0) * G * vf * dkNr / pow(denom1, 2.0) \
        - (3.0 / 8) * e4A * x * pow(X1,4.0) * pow(X2 * X5,2.0) * G * vf * dkNr / (X4 * pow(denom1,2.0)) \
        + e4A * x * pow(X1,4.0) * X4 * pow(X5,2.0) * G * vf * dkNr / (6 * pow(denom1,2.0)) \
        + (3.0 / 8) * e2A * x * pow(X1,2.0) * X2 * pow(X5,4.0) * kNr * G * vf * dkNr / pow(denom1,2.0) \
        + e2A * x * pow(X1,2.0) * X4 * pow(X5, 4.0) * kNr * G * vf * dkNr / (12.0 * pow(denom1,2.0)) \
        - (3.0 / 32) * x * pow(X1 * X2 * dkNr, 2.0) * pow(X5,4.0) * G * vf / pow(denom1, 2.0) \
        + (3.0 / 16) * e4A * x * pow(X1,4.0) * pow(X2 * X5,2.0 ) * G * vf * d2kNr / pow(denom1,2.0) \
        + (3.0 / 16) * e2A * x * pow(X1*X2,2.0) * pow(X5,4.0) * kNr * G * vf * d2kNr / pow(denom1,2.0) \
        + 0.75 * x * e2A * e4A * pow(X1,6.0) * X2 * G * dvfdl / pow(denom1, 2.0) \
        + 0.75 * x * e2A * e4A * pow(X1, 6.0) * pow(X2, 2.0) * G * dvfdl / (X4 * pow(denom1,2.0)) \
        + 1.5 * e4A * x * pow(X1,4.0) * X2 * pow(X5,2.0) * kNr * G * dvfdl / pow(denom1, 2.0) \
        + 0.75 * e4A * x * pow(X1, 4.0) * pow(X2 * X5, 2.0) * kNr * G * dvfdl / (X4 * pow(denom1, 2.0)) \
        + e4A * x * pow(X1,4.0) * X4 * pow(X5,2.0) * kNr * G * dvfdl / (6.0 * pow(denom1,2.0)) \
        + 0.75 * e2A * x * pow(X1 * kNr,2.0) * X2 * pow(X5,4.0) * G * dvfdl / pow(denom1,2.0) \
        + e2A * x * pow(X1*kNr,2.0) * X4 * pow(X5,4.0) * G * dvfdl / (6.0 * pow(denom1,2.0)) \
        + (3.0/8) * e4A * x * pow(X1,4.0) * pow(X2 * X5, 2.0) * G * dkNr * dvfdl / pow(denom1,2.0) \
        + (3.0 / 8) * e2A * x * pow(X1 * X2,2.0) * pow(X5, 4.0) * kNr * G * dkNr * dvfdl / pow(denom1,2.0) \
        + (3.0 / 8) * e2A * e4A * x * pow(X1,6.0) * pow(X2,2.0) * G * d2vfdl / pow(denom1,2.0) \
        + 0.75 * e4A * x * pow(X1,4.0) * pow(X2 * X5,2.0) * kNr * G * d2vfdl / pow(denom1,2.0) \
        + (3.0 / 8) * e2A * x * pow(X1 * X2 * kNr,2.0) * pow(X5, 4.0) * G * d2vfdl / pow(denom1,2.0) ; 
        J( 3 , 2 ) = 0.75 * e2A * x * pow(X1 * X2,2.0) * dvfdt / (X4 * G) \
        + x * X4 * pow(X5,2.0) * kNr * dvfdt / (6.0 * G) \
        + (3.0 / 16) * x * pow(X2 * X5,2.0) * dkNr * dvfdt / G \
        + (3.0 / 8.0) * e2A * x * pow(X1 * X2,2.0) * d2vfdldt / G \
        + (3.0 / 8.0) * x * pow(X2 * X5, 2.0) * kNr * d2vfdldt / G;
        J( 3 , 3 ) = - 5 - 9.0 * pow(X2/X4,2.0) + (4.0 / 3) * pow(X4/X2,2.0) \
        - 0.75 * e2A * x * pow(X1 * X2 / X4, 2.0) * vf / G + x * pow(X5,2.0) * kNr * vf / (6.0 * G) \
        + 0.75 * e2A * pow(X1 * X2 / X4,2.0) * vg ;
        J( 3 , 4 ) = -0.75 * e2A * x * pow(X1 * X2,2.0) * X5 * kNr * vf / ( X4 * denom1 * G) \
        + e2A * x * pow(X1,2.0) * X4 * X5 * kNr * vf / (3.0 * denom1 * G) \
        + x * X4 * pow(X5, 3.0) * pow(kNr,2.0) * vf / (6.0 * denom1 * G) \
        + (3.0 / 8) * e2A * x * pow(X1 * X2,2.0) * X5 * vf * dkNr / (denom1 * G) \
        + (3.0 / 16) * x * pow(X2,2.0) * pow(X5, 3.0) * kNr * vf * dkNr / (denom1 * G) \
        + (3.0 / 8) * e2A * x * pow(X1 * X2,2.0) * X5 * kNr * dvfdl / (denom1 * G) \
        + (3.0 / 8) * x * pow(X2,2.0) * pow(X5,3.0) * pow(kNr,2.0) * dvfdl / (denom1 * G);
        J( 4 , 0 ) = 8.0 * e2A * X1 * pow(X5,3.0) * kNr / pow(denom1,2.0) \
        + 16.0 * pow(X5,5.0) * pow(kNr,2.0) / (X1 * pow(denom1,2.0)) \
        + 8 * pow(X5, 7.0) * pow(kNr / X1, 3.0) / (e2A * pow(denom1, 2.0)) \
        + e2A * x * X1 * pow(X5,5.0) * pow(kNr,2.0) * G * vf / (6.0 * pow(denom1,2.0)) \
        + X4 * pow(X5 / X1,3.0) * dkNr / e2A + 2 * X1 * e2A * dlogvfdt / kNr \
        + 2 * e2A * X1 * X4 * pow(X5, 3.0) * kNr * dlogvfdl / pow(denom1,2.0) \
        + 4 * X4 * pow(X5,5.0) * pow(kNr,2.0) * dlogvfdl / (X1 * pow(denom1,2.0)) \
        + 2.0 * X4 * pow(X5,7.0) * pow(kNr / X1 ,3.0) * dlogvfdl / (e2A * pow(denom1, 2.0) );
        J(4,1) = - (8.0 / 9.0) * pow(X4,2.0) * X5 / pow(X2,3.0) - 4.0 * pow(X5,3.0) * dkNr / (e2A * pow(X1,2.0)) \
        - e2A * x * pow(X1 / denom1, 2.0) * pow(X5,5.0) * kNr * G * vf * dkNr / 12.0 \
        + e2A * x * pow(X1, 2.0) * pow(X5,3.0) * G * vf * dkNr / (6.0 * denom1) + X4 * X5 * pow(dkNr/kNr,2.0)\
        - X4 * pow(X5,3.0) * d2kNr / (2.0 * e2A * pow(X1,2.0)) - X4 * X5 * d2kNr / kNr \
        - e2A * pow(X1,2.0) * dkNr * dlogvfdt / pow(kNr,2.0) + e2A * x * pow(X1,2.0) * pow(X5,3.0) * kNr * G * dvfdl / (6 * denom1) \
        - X4 * pow(X5,3.0) * dkNr * dlogvfdl / (e2A * pow(X1,2.0)) \
        - pow(X5, 2.0) * dlogvfdl * dlogvfdt - e2A * pow(X1, 2.0) * dlogvfdl * dlogvfdt / kNr \
        +  X4 * X5 * pow(dlogvfdl,2.0) + X4 * pow(X5, 3.0) * kNr * pow(dlogvfdl/X1,2.0) / e2A + pow(X5,2.0) * d2vfdldt_vf \
        + e2A * pow(X1, 2.0) * d2vfdldt_vf / kNr - X4 * X5 * d2vfdl_vf - X4 * pow(X5,3.0) * kNr * d2vfdl_vf / (e2A * X1 * X1);
        J(4,2) = e2A * x * pow(X1,2.0) * pow(X5,3.0) * kNr *G * dvfdt / (6.0 * denom1) - pow(X5,2.0) * pow(dlogvfdt,2.0) \
        - e2A * pow(dlogvfdt * X1, 2.0) / kNr + pow(X5,2.0) * d2vfdt_vf + e2A * pow(X1,2.0) * d2vfdt_vf / kNr \
        + X4 * X5 * dlogvfdl * dlogvfdt + X4 * pow(X5, 3.0) * kNr * dlogvfdl * dlogvfdt / (e2A * X1 * X1) \
        - X4 * X5 * d2vfdldt_vf - X4 * pow(X5,3.0) * kNr * d2vfdldt_vf / (e2A * X1 * X1 ) ;
        J( 4 , 3 ) = (8.0 / 9) * X4 * X5 / pow(X2,2.0) - pow(X5,3.0) * dkNr / (2.0 * e2A * pow(X1,2.0)) \
        - X5 * dkNr / kNr - X5 * dlogvfdl - pow(X5,3.0) * kNr * dlogvfdl / (e2A * pow(X1,2.0));
        J( 4 , 4 ) = - 4 + (4.0 / 9) * pow(X4 / X2,2.0) - 12 * pow(X5 / X1,2.0) * kNr / e2A \
        - e2A * x * pow(X1 * kNr / denom1,2.0) * pow(X5,4.0) * G * vf / 6.0 \
        + e2A * x * pow(X1 * X5,2.0) * kNr * G * vf / (2.0 * denom1) - 1.5 * X4 * pow(X5/X1,2.0) * dkNr / e2A \
        - X4 * dkNr / kNr + 2.0 * X5 * dlogvfdt - X4 * dlogvfdl - 3 * X4 * pow(X5 / X1,2.0) * kNr * dlogvfdl / e2A;
        dfdt[0] = x * X1 * pow(X5,4.0) * pow(kNr, 2.0) * vf / ( 6 * denom1 * G ) ;
        dfdt[1] = 0.0;
        dfdt[2] = 0.0;
        dfdt[3] = 3.0 * e4A * x * pow(X1,4.0) * pow(X2,2.0) * vf / ( 2 * X4 * denom1 * G) \
        + (9.0 / 4) * x * e2A * pow(X1 * X2 * X5,2.0) * kNr * vf / (X4 * denom1 * G ) \
        + x * X4 * pow( X5 ,4.0) * pow(kNr,2.0) * vf / ( 6 * denom1 * G ) \
        - 1.5 * e4A * pow(X1,4.0) * pow(X2,2.0) * vg / ( X4 * denom1 ) \
        - 1.5 * e2A * pow( X1 * X2 * X5,2.0) * kNr * vg / ( X4 * denom1 ) \
        - 0.75 * e4A * pow(X1, 4.0) * pow(X2,2.0) * dvg / denom1 \
        - 0.75 * e2A * pow(X1 * X2 * X5, 2.0) * kNr * dvg / denom1 \
        + (3.0 / 16) * x * pow(X2,2.0) * pow(X5,4.0) * kNr * vf * dkNr / (denom1 * G) \
        + 0.75 * e4A * x * pow(X1,4.0) * pow(X2,2.0) * dvfdl / (denom1 * G) \
        + (9.0 / 8) * e2A * x * pow(X1 * X2 * X5,2.0) * kNr * dvfdl / (denom1 * G) \
        + (3.0 / 8) * x * pow(X2,2.0) * pow(X5, 4) * kNr * kNr * dvfdl / (denom1 * G)  ;
        dfdt[4] = 8.0 * e2A * pow(X1, 4.0) * pow(X5,3.0) * kNr / pow(e2A * pow(X1,3.0) + X1 * pow(X5,2.0) * kNr,2.0) \
        + 16 * pow(X1 * kNr,2.0) * pow(X5, 5.0) / pow(e2A * pow(X1,3.0) + X1 * pow(X5,2.0) * kNr,2.0) \
        + (8.0/e2A) * pow(X5,7.0) * pow(kNr,3.0) / pow(e2A * pow(X1,3.0) + X1 * pow(X5,2.0) * kNr,2.0) \
        - e2A * x * pow(X1, 2.0) * pow(X5, 5.0) * pow(kNr, 2.0) * G * vf / (6.0 * pow(e2A * pow(X1,2.0) + pow(X5,2.0) * kNr,2.0)) \
        + e2A * x * pow(X1, 4.0) * pow(X5, 5.0) * pow(kNr,2.0) * G * vf / (3.0 * pow(e2A * pow(X1,3.0) + X1 * pow(X5,2.0) * kNr,2.0)) \
        + X4 * pow(X5,3.0) * dkNr / (e2A * pow(X1,2.0)) + 2 * e2A * pow(X1,2.0) * dlogvfdt / kNr \
        + 2 * e2A * pow(X1, 4.0) * X4 * pow(X5,3.0) * kNr * dlogvfdl / pow(e2A * pow(X1,3.0) + X1 * pow(X5,2.0) * kNr,2.0) \
        + 4.0 * pow(X1 * kNr,2.0) * X4 * pow(X5,5.0) * dlogvfdl / pow(e2A * pow(X1,3.0) + X1 * pow(X5,2.0) * kNr,2.0) \
        + (2.0 / e2A ) * X4 * pow(X5, 7.0) * pow(kNr,3.0) * dlogvfdl / pow(e2A * pow(X1,3.0) + X1 * pow(X5,2.0) * kNr,2.0) ;
        
    }
};

// Solve VHQCD and get the background fields
void solveVHQCD(double xi, double ti)
{
    // Computes dr/dA, tau, lambda, dtau/dA and dlambda/dA given x and tau0
    // x - double. Physically it means x = N_f / N_c when N_f, N_c -> inf but with fixed coefficient
    // tau0 - > parameter related with the exponential behaviour of the tachyon in the IR
    // Does not return anything. Just creates a file with all the relevant data.
    
    // Create file where data will be saved
    ofstream myfile;
    myfile.open ("VQCD.txt");
    myfile << "A" << '\t' << "dz" << '\t' << "lambda" << '\t' << "tau" << '\t' << "dlambda" << '\t' << "dtau"<< endl;
    // Boundary conditions in the IR
    double zIR = log(70.0 / ti) / CI(xi) ;
    double aIR = AIR(zIR) ;
    double lambdair = lambdaIR(zIR) ;
    double tauir = 70.0 ;
    // dA/dz at zIR
    double daIR =  173.0 / ( 1728.0 * pow(zIR, 3.0)) + 0.5 /zIR - 2.0 * zIR  ;
    // dtau/dz at zIR
    double dtauir = CI(xi) * tauir ;
    // dlambda/dz at zIR. Expression got from eom2 in Mathematica notebook.
    double dlambdair = sqrt(1.5) * lambdair * sqrt( 6 * pow(daIR, 2.0) + exp( 2.0 * aIR) * xi * Vf(lambdair,tauir, xi) \
         / (2.0 * sqrt(1+ dtauir * dtauir * k(lambdair, xi) / exp(2 * aIR ) ) ) - 0.5 * exp( 2.0 * aIR) * Vg(lambdair) );
    double dzIR = 1.0 / daIR ;
    dlambdair = dlambdair / daIR;
    dtauir = dtauir / daIR;

    // Define the type of the state. We have X = {dz, lambda, tau, dlambda, dtau}
    state_type X (5);
    X <<= dzIR, lambdair, tauir, dlambdair, dtauir;

    cout << "Solving VHQCD for x = " << xi << ", tau0 = " << ti << endl;
    double AUV = 5.0;
    double h = 0.1 ;
    
    auto stepper = make_dense_output( 1.0e-6 , 1.0e-6 , rosenbrock4< double >() );
    stepper.initialize( X , aIR , h );
    double relDiff = 1e6 ;
    double tau1 = tauir ;
    while ( relDiff > 1e-6 )
    {
        myfile << stepper.current_time() << '\t' << X(0) << '\t' << X(1) << '\t' << X(2) << '\t' << X(3) << '\t' << X(4) << endl;
        stepper.do_step( make_pair( VQCD(xi, ti) , jacobian(xi) ) ) ;
        X = stepper.current_state();
        double tau2 = X(2);
        relDiff = abs(tau2 - tau1) / tau1 ;
        tau1 = tau2 ;
    }
    myfile.close();
} ;


int main(int argc, char * argv[])
{
    if (argc != 3) cout << "Invalid number of parameters. Insert the values of x and tau0." << endl;
    // Get the parameters parsed through CMD
    double x = stod(argv[1]) ;
    double tau0 = stod(argv[2]);
    solveVHQCD(x, tau0);

    return 0;
}