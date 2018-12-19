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
typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;

struct output_observer
{
    vector< vector< double > > &data;
    output_observer( vector< vector < double > > &fieldsData ) : data( fieldsData ) { }
    void operator()( const vector_type &X , double A ) const
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
    void operator()(const vector_type &X , vector_type &dXdt , double A)
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
    void operator()( const vector_type & X, matrix_type &J , const double &A , vector_type &dfdt )
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
        /*cout << "Jacobian : " << endl;
        cout << J(0,0) << '\t' << J(0,1) << '\t' << J(0,2) << '\t' << J(0,3) << '\t' << J(0,4) << endl;
        cout << J(1,0) << '\t' << J(1,1) << '\t' << J(1,2) << '\t' << J(1,3) << '\t' << J(1,4) << endl;
        cout << J(2,0) << '\t' << J(2,1) << '\t' << J(2,2) << '\t' << J(2,3) << '\t' << J(2,4) << endl;
        cout << J(3,0) << '\t' << J(3,1) << '\t' << J(3,2) << '\t' << J(3,3) << '\t' << J(3,4) << endl;
        cout << J(4,0) << '\t' << J(4,1) << '\t' << J(4,2) << '\t' << J(4,3) << '\t' << J(4,4) << endl;
        cout << "dfdt : " << endl;
        cout << dfdt[0] << '\t' << dfdt[1] << '\t' << dfdt[2] << '\t' << dfdt[3] << '\t' << dfdt[4] << endl;*/
    }
};

// Solve VHQCD and get the background fields
vector<vector<double> > solveVHQCD(double xi, double ti)
{
    vector<vector<double> > bckFields ;
    vector < double > As ;
    vector <double> dzs ;
    vector <double> ls ;
    vector <double> taus ;
    vector <double> dls;
    vector <double> dtaus;
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
    vector_type X (5);
    X <<= dzIR, lambdair, tauir, dlambdair, dtauir;
    As.push_back(aIR);
    dzs.push_back(dzIR);
    ls.push_back(lambdair);
    taus.push_back(tauir);
    dls.push_back(dlambdair);
    dtaus.push_back(dtauir);
    bckFields.push_back(As);
    bckFields.push_back(dzs);
    bckFields.push_back(ls);
    bckFields.push_back(taus);
    bckFields.push_back(dls);
    bckFields.push_back(dtaus);
    cout << "Solving VHQCD for x = " << xi << ", tau0 = " << ti << endl;
    double AUV = 2.5;
    double h = 0.01 ;
    
    /*typedef rosenbrock4< double > stepper_type;
    stepper_type stepper;
    typedef stepper_type::state_type state_type;
    typedef stepper_type::value_type stepper_value_type;
    typedef stepper_type::deriv_type deriv_type;
    typedef stepper_type::time_type time_type;
    state_type X(5) , Xerr(5);
    X(0) = dzIR; X(1) = lambdair; X(2) = tauir; X(3) = dlambdair; X(4) = dtauir;
    stepper.do_step( make_pair( VQCD(xi, ti) , jacobian(xi) ) , X , aIR , 0.01 , Xerr );
    cout << X(0) << '\t' << X(1) << '\t' << X(2) << '\t' << X(3) << '\t' << X(4) << endl;
    stepper.do_step( make_pair( VQCD(xi, ti) , jacobian(xi) ) , X , aIR , 0.01 );
    cout << X(0) << '\t' << X(1) << '\t' << X(2) << '\t' << X(3) << '\t' << X(4) << endl;*/
    size_t num_of_steps = integrate_adaptive( make_dense_output< rosenbrock4< double > >( 1.0e-6 , 1.0e-6 ) ,
            make_pair( VQCD(xi, ti) , jacobian(xi) ), X , aIR , AUV , h, output_observer(bckFields) );
    clog << num_of_steps << endl;
    return bckFields;
} ;


int main()
{
    vector<vector<double> > solution = solveVHQCD(1.0,1.0);
    ofstream myfile;
    myfile.open ("VQCD.txt");
    myfile << "A" << '\t' << "dz" << '\t' << "lambda" << '\t' << "tau" << '\t' << "dlambda" << '\t' << "dtau"<< endl;
    int n = solution[0].size() ;
    for(int i = 0 ; i < n ; i++)
    {
        myfile << solution[0][i] << '\t' << solution[1][i] << '\t' << solution[2][i] << '\t' << solution[3][i] << '\t' << solution[4][i] << endl;
    }
    myfile.close();
    return 0;
}