#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

class IHQCD{
    // Initial value of the warp factor A
    float A0;
    // Initial value of lambda
    float lambda0;
    // The geometry is defined until zmax
    float zmax;
    // Spacing between adjacent points
    float h;

    public :
      // Initialize with default values
      IHQCD ()
      {
        A0 = 5;
        lambda0 = 0.0337462;
        zmax = 10;
        h = 0.002;
      } ;
      // Initialize with custom values
      IHQCD (float a0, float l0, float zMax, float spacing)
      {
        A0 = a0;
        lambda0 = l0;
        zmax = zMax;
        h = spacing;
      }
      // Solve IHQCD and get the background fields
      vector <vector<float> > solveIHQCD()
      {
        vector <vector <float> > bckFields ;
        // Definition of parameters needed to solve the ODEs
        float b0 = 4.2 ;
        float b1 = 51 * pow(b0,2) /121;
        float alpha = 2;
        float aa = (3/8) * (alpha - 1) / alpha ;
        float LAdS = 1;
        float zmin = exp(-A0);
        cout << "Solving IHQCD for A0 = " << A0 << ", lambda0 = " << lambda0 << ", zmin = " << zmin << ", zmax = " << zmax << ", h = " << h << endl;
        return bckFields;
      } ;
};

int main(){
    IHQCD ihqcd = IHQCD();
    ihqcd.solveIHQCD();
    return 0;
}