#ifndef ROOOT
#define ROOOT

#include <vector>
using namespace std; 

template <class T>
bool bracketRoot1(T &f, double &x1, double &x2)
{
    const int nTrials = 50;
    int i = 0 ;
    const float ScaleFactor = 1.6 ;
    if ( x1 == x2) throw("Bad guess of x1 and x2!") ;
    double f1 = f(x1);
    double f2 = f(x2);
    while (i <= nTrials)
    {
        // Check if there is a zero between x1 and x2
        if (f1 * f2 < 0) return true;
        // If there is no zero between x1 and x2 we change x1 and x2 until we find one
        if (abs(f1) < abs(f2))
        {
            x1 += ScaleFactor * ( x1 - x2 ) ;
            f1 = f(x1);
        }
        else
        {
            x2 += ScaleFactor * ( x2 - x1 );
            f2 = f(x2) ;
        }
        i += 1;
    }
    return false;
}

template <class T>
void bracketRoot2(T &f, double x1, double x2, const int n, int &nRoots, vector<double> &xlb, vector<double> &xub)
{
    // Given a function or functor f it checks if in each of the  n intervals between x1 and x2 there is a root of it.
    // If a root exists it adds the lower bound to the vector xlb and the upper bound to xup and updates nRoots.
    // It is assumed x1 < x2
    // nRoots tells us how many roots were found in that interval
    // f - function or functor to inspect if it has zeros
    // x1 - lowest lower bound where we want to search for zeros. It's a double
    // x2 - biggest upper bound where we want to search for zeros. It's a double
    // n  - integer number that represents in how many subintervals we want to search for a zero
    // nRoots - integer number that tells how many roots were found.
    // xlb - empty vector of doubles
    // xup - empty vector of doubles
    nRoots = 0;
    // spacing of each interval
    double dx = (x2 - x1) / n ;
    double f1, f2;
    // Iterate over the n different subsets between x1 and x2
    for(int i = 0; i < n; i++)
    {
        x2 = x1 + dx ;
        f1 = f(x1) ;
        f2 = f(x2) ;
        if (f1 * f2 < 0)
        {
            xlb.push_back(x1) ;
            xub.push_back(x2) ;
            nRoots += 1 ;
        }
        x1 = x2;
    }

}

template <class T>
double findRootBisection(T &f, const double x1, const double x2, const double xacc)
{
    // Computes the single root between x1 and x2 to accuracy xacc of the function f using the bisection method
    // x1 - lower limit for the root, double
    // x2 - upper limit for the root, double
    // xacc - desired accuracy for the root
    double f1 = f(x1);
    double f2 = f(x2);
    if ( f1 * f2 >= 0) throw("No root of the function in the given interval") ;
    // Number of maximum iterations
    int nItMax = 50 ;
    double fc, xc, xa, xb, tol ;
    xa = x1 ;
    xb = x2 ;
    xc = 0.5 * ( xa + xb );
    fc = f(xc) ;
    tol = 1000000000;
    for (int i = 0; i < nItMax ; i++)
    {
        if ( f1 * fc < 0)
        {
            f2 = fc ;
            xb = xc ;
            xc = 0.5 * ( xa + xb ) ;
            fc = f(xc) ;
            tol = 0.5 * abs( xb - xa );
        }

        else
        {
            f1 = fc ;
            xa = xc ;
            xc = 0.5 * ( xa + xb) ;
            fc = f(xc) ;
            tol = 0.5 * abs( xb - xa );
        }
        if (abs(tol) < xacc || fc == 0) return xc ;
    }
    throw("The algorithm didn't converge.") ;
}

// The following two methods are good if the function is sufficiently smooth
template <class T>
double findRootFalsePosition(T &f, const double x1, const double x2, const double xacc)
{
    // Computes the single root between x1 and x2 to accuracy xacc of the function f using the False Position method
    // x1 - lower limit for the root, double
    // x2 - upper limit for the root, double
    // xacc - desired accuracy for the root
    // Number of maximum iterations
    int nItMax = 100 ;
    double f1 = f(x1) ;
    double f2 = f(x2) ;
    if (f1 * f2 > 0.0) throw("No root in the given interval.") ;
    double xa = x1 ;
    double xb = x2 ;
    double tol = 1000000000;
    double root, froot;
    for (int n = 0; n < nItMax; n++)
    {
        root = (xb * f1 - xa * f2) / (f1 - f2) ;
        froot = f(root) ;
        if ( f1 * froot < 0 )
        {
            xb = root ;
            f2 = froot ;
        }
        else
        {
            xa = root;
            f1 = froot ;
        }
        tol = 0.5 * (xb - xa) ;
        if (abs(tol) < xacc || froot == 0) return root;
    }
    throw("Root not found!") ;
}

template <class T>
double findRootSecantMethod(T &f, const double x1, const double x2, const double xacc)
{
    // Computes the single root between x1 and x2 to accuracy xacc of the function f using the False Position method
    // x1 - lower limit for the root, double
    // x2 - upper limit for the root, double
    // xacc - desired accuracy for the root
    // Number of maximum iterations
    int nItMax = 100 ;
    double f1 = f(x1) ;
    double f2 = f(x2) ;
    if (f1 * f2 > 0.0) throw("No root in the given interval.") ;
    double xp = x1;
    double root = x2;
    double tol = 100000000000;
    double xn;
    for(int n = 0; n < nItMax; n++)
    {
        xn = (root * f1 - xp * f2) / (f1 - f2) ;
        xp = root;
        root = xn ;
        f2 = f(root) ;
        f1 = f(xp) ;
        tol = 0.5 * (root - xp);
        if (abs(tol) < xacc || f2 == 0) return root ;
    }
    throw("No root found.") ;
}

#endif