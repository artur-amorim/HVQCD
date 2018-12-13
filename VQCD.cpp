#include <iostream>
#include <vector>
#include <math.h>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <fstream>

using namespace std;
using namespace boost::numeric::odeint;

// Define the type of the state. We have X = {dz, lambda, tau, dlambda, dtau}
typedef boost::array< long double , 5 > state_type ;

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
        cout << A << '\t' << X[0] << '\t' <<  X[1] << '\t'<< X[2] << '\t' << X[3] << '\t' << X[4] << endl;
    }
};

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
long double Vf (long double l, long double t, long double x) 
{
    return exp( log( Vf0(l,x) ) - a0(x) * pow(t, 2.0) );
}
long double dVfdl(long double l, long double t, long double x)
{
    return exp( log( W0 * ( W1(x) + 2.0 * l * W2(x) )) - a0(x) * t * t);
}
long double dLogVfdt (long double t, long double x)
{
    return - 2.0 * t * a0(x) ;
}
long double dLogVfdl(long double l, long double x)
{
    return ( W1(x) + 2.0 * l * W2(x) ) / ( 1.0 + l * W1(x) + W2(x) * l * l) ;
}

// Definition of the IR assymptotics
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

// Solve VHQCD and get the background fields
vector<vector<long double> > solveVHQCD(long double xi, long double ti)
{
    vector<vector<long double> > bckFields ;
    vector <long double> As ;
    vector <long double> dzs ;
    vector <long double> ls ;
    vector <long double> taus ;
    vector <long double> dls;
    vector <long double> dtaus;
    // Boundary conditions in the IR
    VQCD model = VQCD(xi, ti);
    long double zIR = log(70.0 / model.tau0) / CI(model.x) ;
    long double aIR = AIR(zIR) ;
    long double lambdair = lambdaIR(zIR) ;
    long double tauir = 70.0 ;
    // dA/dz at zIR
    long double daIR =  173.0 / ( 1728.0 * pow(zIR, 3.0)) + 0.5 /zIR - 2.0 * zIR  ;
    // dtau/dz at zIR
    long double dtauir = CI(model.x) * tauir ;
    // dlambda/dz at zIR. Expression got from eom2 in Mathematica notebook.
    long double dlambdair = sqrt(1.5) * lambdair * sqrt( 6 * pow(daIR, 2.0) + exp( 2.0 * aIR) * model.x * Vf(lambdair,tauir,model.x) \
         / (2.0 * sqrt(1+ dtauir * dtauir * k(lambdair,model.x) / exp(2 * aIR ) ) ) - 0.5 * exp( 2.0 * aIR) * Vg(lambdair) );
    long double dzIR = 1.0 / daIR ;
    dlambdair = dlambdair / daIR;
    dtauir = dtauir / daIR;
    state_type X = {dzIR, lambdair, tauir, dlambdair, dtauir};
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
    cout << "Solving VHQCD for x = " << model.x << ", tau0 = " << model.tau0 << endl;
    long double AUV = 2.5;
    long double h = 0.00001 ;
    //runge_kutta4< state_type > stepper;
    //integrate_const( stepper , model , X , aIR , AUV , h, output_observer(bckFields) );
    integrate(model, X , aIR , AUV , h , output_observer( bckFields) );
    return bckFields;
} ;


int main()
{
    vector<vector<long double> > solution = solveVHQCD(1.0,1.0);
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