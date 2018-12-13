#include <iostream>
#include <vector>
#include <math.h>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <fstream>

using namespace std;
using namespace boost::numeric::odeint;

// Define the type of the state 
typedef boost::array< long double , 2 > state_type;

struct output_observer
{
    vector< vector< long double > > &data;
    output_observer( vector< vector < long double > > &fieldsData ) : data( fieldsData ) { }
    void operator()( const state_type &x , long double z ) const
    {
        // Add z to vector of zs
        data[0].push_back(z);
        // Add A to the vector of As
        data[1].push_back(x[0]);
        // Add lambda to the vector of lambdas
        data[2].push_back(x[1]);
        //cout << z << '\t' << x[0] << '\t' <<  x[1] << endl;
    }
};

// Definition of constants need to solve ODE
  const long double b0 = 4.2 ;
  const long double b1 = 51.0 * b0 * b0 /121.0;
  const long double alpha = 2.0;
  const long double aa = (3.0/8) * (alpha - 1) / alpha ;
  const long double LAdS = 1.0 ;

  // Definition of the system of ODEs necessary to solve IHQCD
  void odeSys( const state_type &x , state_type &dxdt , long double z )
  {
    // Given the pointer of the stat x, it updates the rate of change dxdt at z
    long double X = -b0 * x[1] / (3 + 2 * b0 * x[1]) - (2 * b0 * b0 + 3 * b1) * x[1] * x[1] / (9 * (1 + x[1] * x[1]) * (1 + (2 * b0 *b0 + 3 * b1) * log(1 + x[1] * x[1]) / (18 * aa))) ;
    long double W = (9.0/4) * pow(1 + (2.0/3) * b0 * x[1],2.0/3) * pow((1 + (2 * b0 * b0 + 3 * b1) * log(1 + x[1] * x[1]) / (18 * aa)),4.0 * aa / 3) ;
    dxdt[0] = -(4./9) * W * exp(x[0]) / LAdS;
    dxdt[1] = -(4./3) * X * W * x[1] * exp(x[0]) / LAdS;
  }

class IHQCD{
    // Initial value of the warp factor A
    long double A0;
    // Initial value of lambda
    long double lambda0;
    // The geometry is defined until zmax
    long double zmax;
    // Spacing between adjacent points
    long double h;

    public :
      // Initialize with default values
      IHQCD ()
      {
        A0 = 5.0;
        lambda0 = 0.0337462;
        zmax = 10.0;
        h = 0.002;
      } ;
      // Initialize with custom values except b0, b1, alpha, aa and LAdS
      IHQCD (long double a0, long double l0, long double zMax, long double spacing)
      {
        A0 = a0;
        lambda0 = l0;
        zmax = zMax;
        h = spacing;
      }
      
      // Solve IHQCD and get the background fields
      vector<vector<long double> > solveIHQCD()
      {
        vector<vector<long double> > bckFields ;
        vector <long double> zs ;
        vector <long double> As ;
        vector <long double> ls ;
        long double zmin = exp(-A0);
        zs.push_back(zmin);
        // Initial conditions
        state_type x = {A0, lambda0};
        As.push_back(x[0]);
        ls.push_back(x[1]);
        bckFields.push_back(zs);
        bckFields.push_back(As);
        bckFields.push_back(ls);
        cout << "Solving IHQCD for A0 = " << A0 << ", lambda0 = " << lambda0 << ", zmin = " << zmin << ", zmax = " << zmax << ", h = " << h << endl;
        runge_kutta4< state_type > stepper;
        integrate_const( stepper , odeSys , x , zmin , zmax , h, output_observer(bckFields) );
        //integrate(odeSys, x , zmin , zmax , h , output_observer( bckFields) );
        return bckFields;
      } ;
};

int main()
{
  IHQCD ihqcd = IHQCD();
  vector<vector<long double> > solution;
  solution = ihqcd.solveIHQCD();
  ofstream myfile;
  myfile.open ("data.txt");
  myfile << "z" << '\t' << "A" << '\t' << "lambda" << endl;
  int n = solution[0].size() ;
  for(int i = 0 ; i < n ; i++)
  {
    myfile << solution[0][i] << '\t' << solution[1][i] << '\t' << solution[2][i] << endl;
  }
  myfile.close();
  return 0;
}