#ifndef __BCK_POTENTIALS_ASYMP_HPP
#define __BCK_POTENTIALS_ASYMP_HPP

#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

// Definition of constants need to solve ODE
const long double V0 = 12.0 ;
const long double V1 = 11.0 / ( 27 * M_PI * M_PI);
const long double V2 = 4619.0 / ( 46656.0 * pow(M_PI,4.0) );
const long double W0 = 12.0/11;
const long double l0 = 8 * pow(M_PI, 2.0);
// Definition of functions needed to solve ODE
long double W1 (long double x) 
{
    return (24 + (11 - 2 * x) * W0) / ( 27.0 * pow(M_PI,2.0) * W0 ) ; 
}
long double W2 (long double x) 
{
    return (24 * (857 - 46 * x) + (4619 - 1714 * x + 92 * pow(x, 2.0)) * W0) / (46656 * pow(M_PI,4.0) * W0) ;
}
long double a0 (long double x) 
{
    return (12 - x * W0) / 8.0 ;
}
long double a1 (long double x) 
{
    return (115 - 16 * x) / (216.0 * pow(M_PI, 2.0)) ;
}
long double Vg (long double l) 
{
    return V0 * (1+ V1 * l + V2 * pow(l,2.0) * sqrt(1 + log(1.0 + l / l0))/pow(1.0 + l / l0,2.0/3)) ;
}
long double dVg (long double l)
{
    long double ans = V1 + V2 * pow(l, 2.0) / ( 2.0 * pow( 1 + l / l0, 5.0/3) * l0 * sqrt( 1 + log( 1 + l / l0 )));
    ans += 2.0 * V2 * l * sqrt( 1 + log( 1 + l / l0 )) / pow( 1 + l / l0, 2.0/3);
    ans += - 2.0 * V2 * pow(l, 2.0) * sqrt( 1 + log( 1 + l / l0 )) / ( 3.0 * pow( 1 + l / l0, 5.0/3) * l0 );
    return V0 * ans;
}
long double d2Vg( long double l)
{
    long double ans = V0 * V2 ;
    ans = ans * (37 * l * l + 120.0 * l * l0 + 72 * l0 *l0 + 2.0 * log(1.0 + l / l0) * (31 * l * l + 84.0 * l * l0  + 72.0 * l0 * l0) + 8.0 * (2 * l * l + 6 * l * l0 + 9 * l0 * l0) * pow(log(1.0 + l / l0),2.0))  ;
    ans = ans / (36.0 * l0 * l0 * pow(1.0  + l / l0, 8.0/3.0) * pow(1 + log(1.0 + l / l0),1.5)) ;
    return ans;
}
long double Vf0 (long double l, long double x) 
{
    return W0 * ( 1 + W1(x) * l + W2(x) * pow(l, 2.0)) ;
}
long double k (long double l, long double x) 
{
    return 1.0 / pow(1 + 3 * a1(x) * l / 4.0, 4.0 / 3) ;
}
long double dk(long double l, long double x)
{
    return - a1(x) / pow( 1 + 3 * l * a1(x) / 4.0, 7.0 /3 );
}
long double d2k(long double l, long double x)
{
    return 7 * a1(x) * a1(x) / (4.0 * pow(1 + 0.75 * l * a1(x),10.0/3)) ;
}
long double Vf (long double l, long double t, long double x) 
{
    return exp( log( Vf0(l,x) ) - a0(x) * t * t );
}
long double dVfdl(long double l, long double t, long double x)
{
    return exp( log( W0 * ( W1(x) + 2.0 * l * W2(x) )) - a0(x) * t * t);
}
long double d2Vfdl(long double l, long double t, long double x)
{
    return 2.0 * W0 * W2(x) * exp(- a0(x) * t * t);
}
long double d2Vfdl_Vf(long double l, long double t, long double x)
{
    return 2.0 * W0 * W2(x) / Vf0(l,x) ;
}
long double d2Vfdldt(long double l, long double t, long double x)
{
    return W0 * (W1(x) + 2 * l * W2(x)) * (-2.0 * a0(x) * t) * exp(-a0(x) * t * t);
}
long double d2Vfdldt_Vf(long double l, long double t, long double x)
{
    return W0 * (W1(x) + 2 * l * W2(x)) * (-2.0 * a0(x) * t) / Vf0(l,x) ;
}
long double dVfdt(long double l, long double t, long double x)
{
    return  (- 2 * a0(x) * t) * exp( log( Vf0(l,x) ) - a0(x) * pow(t, 2.0) ) ;
}
long double d2Vfdt(long double l, long double t, long double x)
{
    return 2 * W0 * a0(x) * exp(-a0(x)*t*t) * (-1 + 2 * t * t * a0(x)) * (1 + l * W1(x) + l * l * W2(x) );
}
long double d2Vfdt_Vf(long double l, long double t, long double x)
{
    return 2 * a0(x) * (-1 + 2 * t * t * a0(x)) ;
}
long double dLogVfdt (long double t, long double x)
{
    return - 2.0 * t * a0(x) ;
}
long double dLogVfdl(long double l, long double x)
{
    return ( W1(x) + 2.0 * l * W2(x) ) / ( 1.0 + l * W1(x) + W2(x) * l * l) ;
}

// Definition of the IR asymptotics
long double AIR( long double z)
{
    return 13.0 / 8 + log(27 * pow(6, 0.25) / sqrt(4619) ) - (173.0 / 3456 ) / pow(z, 2.0) - pow(z, 2.0) + 0.5 * log(z) ;
}

long double lambdaIR (long double z)
{
    return exp( log(8 * pow(M_PI, 2.0)) - (39.0 / 16) - ( 151.0 / 2304) / pow(z, 2.0) + 1.5 * pow(z, 2.0)) ;
}

long double CI ( long double x)
{
    return 81 * pow(3.0, 5.0 / 6) * pow( 115 - 16 * x, 4.0/3) * (11.0 - x) / ( 812944 * pow(2.0, 1.0 / 6)) ;
}

long double tauIR ( long double z, long double t0, long double x)
{
    return t0 * exp( CI(x) * z) ;
}

#endif