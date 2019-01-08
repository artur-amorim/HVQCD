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
#include "backgroundPotentials.cpp"
#include "backgroundFieldsAsymptotics.cpp"

using namespace std;
using namespace boost::numeric::odeint ;
namespace phoenix = boost::phoenix;

// The stiff algorith only accepts these data types
typedef boost::numeric::ublas::vector< long double > state_type;
typedef boost::numeric::ublas::matrix< long double > matrix_type;

struct output_observer
{
    vector< vector< long double > > &data;
    output_observer( vector< vector < long double > > &fieldsData ) : data( fieldsData ) { }
    void operator()( const state_type &X , long double A ) const
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

struct VQCD
{
    long double x;
    long double tau0;
    VQCD():x(1.0), tau0(1.0){}
    VQCD(long double xi, long double ti):x(xi),tau0(ti){}
    void operator()(const state_type &X , state_type &dXdt , long double A)
    {
        long double e2A = exp(2.0 * A);
        long double kNr = k(X[1], x) ;
        long double dkNr = dk(X[1],x);
        long double vf = Vf(X[1], X[2], x) ;
        long double dvfdl = dVfdl(X[1], X[2], x) ;
        long double dlogvfdt = dLogVfdt(X[2], x) ;
        long double dlogvfdl = dLogVfdl(X[1], x) ;
        long double G = sqrt(1.0 + kNr * pow(X[4], 2.0) / ( e2A * pow(X[0], 2.0 ) ) ) ;
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
    long double x ;
    jacobian( long double xi):x(xi){}
    void operator()( const state_type & X, matrix_type &J , const long double &A , state_type &dfdt )
    {
        long double X1 = X[0] ;
        long double X2 = X[1] ;
        long double X3 = X[2] ;
        long double X4 = X[3] ;
        long double X5 = X[4] ;
        long double e2A = exp(2.0 * A);
        long double e4A = exp(4.0 * A);
        long double kNr = k(X2, x) ;
        long double dkNr = dk(X2,x);
        long double d2kNr = d2k(X2,x) ;
        long double vg = Vg(X2) ;
        long double dvg = dVg(X2);
        long double d2vg = d2Vg(X2) ;
        long double vf = Vf(X2, X3, x) ;
        long double dvfdl = dVfdl(X2, X3, x) ;
        long double d2vfdl = d2Vfdl(X2, X3, x) ;
        long double d2vfdl_vf = d2Vfdl_Vf(X2, X3, x) ;
        long double d2vfdldt = d2Vfdldt(X2, X3, x) ;
        long double d2vfdldt_vf = d2Vfdldt_Vf(X2, X3, x) ;
        long double dvfdt = dVfdt(X2, X3, x) ;
        long double d2vfdt = d2Vfdt(X2, X3, x);
        long double d2vfdt_vf = d2Vfdt_Vf(X2,X3, x) ;
        long double dlogvfdt = dLogVfdt(X3, x) ;
        long double dlogvfdl = dLogVfdl(X2, x) ;
        long double G = sqrt(1.0 + kNr * pow(X5, 2.0) / ( e2A * pow(X1, 2.0 ) ) ) ;
        long double denom1 = e2A * pow(X1,2.0) + pow(X5,2.0) * kNr;
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
void solveVHQCD(long double xi, long double ti)
{
    // Computes dr/dA, tau, lambda, dtau/dA and dlambda/dA given x and tau0
    // x - long double. Physically it means x = N_f / N_c when N_f, N_c -> inf but with fixed coefficient
    // tau0 - > parameter related with the exponential behaviour of the tachyon in the IR
    // Does not return anything. Just creates a file with all the relevant data.
    
    // Create file where data will be saved
    ofstream myfile;
    myfile.open ("VQCD.txt");
    myfile << "A" << '\t' << "dz" << '\t' << "lambda" << '\t' << "tau" << '\t' << "dlambda" << '\t' << "dtau"<< endl;
    // Boundary conditions in the IR
    long double zIR = log(70.0 / ti) / CI(xi) ;
    long double aIR = AIR(zIR) ;
    long double lambdair = lambdaIR(zIR) ;
    long double tauir = 70.0 ;
    // dA/dz at zIR
    long double daIR =  173.0 / ( 1728.0 * pow(zIR, 3.0)) + 0.5 /zIR - 2.0 * zIR  ;
    // dtau/dz at zIR
    long double dtauir = CI(xi) * tauir ;
    // dlambda/dz at zIR. Expression got from eom2 in Mathematica notebook.
    long double dlambdair = sqrt(1.5) * lambdair * sqrt( 6 * pow(daIR, 2.0) + exp( 2.0 * aIR) * xi * Vf(lambdair,tauir, xi) \
         / (2.0 * sqrt(1+ dtauir * dtauir * k(lambdair, xi) / exp(2 * aIR ) ) ) - 0.5 * exp( 2.0 * aIR) * Vg(lambdair) );
    long double dzIR = 1.0 / daIR ;
    dlambdair = dlambdair / daIR;
    dtauir = dtauir / daIR;

    // Define the type of the state. We have X = {dz, lambda, tau, dlambda, dtau}
    state_type X (5);
    X <<= dzIR, lambdair, tauir, dlambdair, dtauir;

    cout << "Solving VHQCD for x = " << xi << ", tau0 = " << ti << endl;
    long double Amax = 5.0;
    long double h = 0.1 ;
    
    auto stepper = make_dense_output( 1.0e-6 , 1.0e-6 , rosenbrock4< long double >() );
    stepper.initialize( X , aIR , h );
    /*dense_output_runge_kutta< controlled_runge_kutta< runge_kutta_dopri5< state_type > > > stepper;
    stepper.initialize( X , aIR , h );*/
    long double relDiff = 1e6 ;
    long double tau1 = tauir ;
    long double A;
    while ( relDiff > 1e-6 /*A < Amax*/ )
    {
        A = stepper.current_time();
        myfile << A << '\t' << X(0) << '\t' << X(1) << '\t' << X(2) << '\t' << X(3) << '\t' << X(4) << endl;
        cout << A << '\t' << X(0) << '\t' << X(1) << '\t' << X(2) << '\t' << X(3) << '\t' << X(4) << endl;
        stepper.do_step( make_pair( VQCD(xi, ti) , jacobian(xi) ) /*VQCD(xi, ti)*/ ) ;
        X = stepper.current_state();
        long double tau2 = X(2);
        relDiff = abs(tau2 - tau1) / tau1 ;
        tau1 = tau2 ;
    }
    /*// Let's initialize a new stepper
    dense_output_runge_kutta< controlled_runge_kutta< runge_kutta_dopri5< state_type > > > stepper1;
    stepper1.initialize( X , A , h );
    while ( A < Amax )
    {
        cout << A << '\t' << X(4) << endl ;
        stepper1.do_step( VQCD(xi, ti) );
        X = stepper1.current_state();
        A = stepper1.current_time();
    }*/
    myfile.close();
} ;

int main(int argc, char * argv[])
{
    if (argc != 3) cout << "Invalid number of parameters. Insert the values of x and tau0." << endl;
    // Get the parameters parsed through CMD
    long double x = stod(argv[1]) ;
    long double tau0 = stod(argv[2]);
    solveVHQCD(x, tau0);

    return 0;
}