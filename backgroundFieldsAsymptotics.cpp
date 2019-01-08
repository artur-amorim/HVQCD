#ifndef __BACKGROUND_FIELDS_ASYMPTOTICS_CPP
#define __BACKGROUND_FIELDS_ASYMPTOTICS_CPP

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

// // Definition of the UV asymptotics
long double V21 (long double x)
{
    return (88.0 - 16 * x) / ( 216.0 * M_PI * M_PI) ;
}

long double V22 (long double x)
{
    return (23.0 + 54 * (34 - 13 * x) / pow(11.0 - 2.0 * x,2.0) ) * pow(V21(x), 2.0) / 64.0 ;
}

long double rho (long double x)
{
    return 9.0 / (22 - 4 * x) ;
}

long double lambdaUV (long double z, long double x, long double Lambda)
{
    return ( - 8.0 / (9.0 * log( z * Lambda ) ) + log( - log( z * Lambda ))  * ( 46.0 / 81 - 128 * V22(x) / (81 * pow(V21(x),2.0)) )/ pow( log( z * Lambda), 2.0 ) ) / V21(x) ;
}

long double tauUV( long double z, long double x, long double Lambda, long double mq, long double sigma)
{
    return mq * z * pow( - log( z * Lambda) , -rho(x) ) + sigma * pow(z, 3.0) * pow( - log( z * Lambda), rho(x) ) ;
}


#endif