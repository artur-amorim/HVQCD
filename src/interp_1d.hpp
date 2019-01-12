#ifndef __INTERP_1D_HPP
#define __INTERP_1D_HPP

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>

using namespace std;

template<class T>
class Base_Interp {
	protected:
		int n, m, jsav, cor, dj;
		const vector<T> *XX, *YY; // pointers to the vectors containing data 
	public:
		Base_Interp(const vector<T> &X, const vector<T> &Y, const int mm)
		{
			n = X.size() ; 
			m = mm ;
			jsav = 0 ;
			cor = 0 ;
			XX = &X ;
			YY = &Y ;
			dj = min(1,(int)pow((T)n,0.25));
		}

		int locate(const T x)
		{
			int ju, jm, jl;
			if ( n < 2 || m < 2 || m > n) throw("locate size error");
			bool ascnd = ( (*XX)[n-1] >= (*XX)[0]) ;
			jl = 0;
			ju = n - 1 ;
			while (ju - jl > 1)
			{
				jm = (ju + jl) / 2 ;
				if ( x >= (*XX)[jm] == ascnd) jl = jm ;
				else ju = jm;
			}
			cor = abs(jl - jsav) > dj ? 0 : 1 ;
			jsav = jl;
			return max(0, min(n-m, jl - (m-2)/2)) ;

		}
		
		int hunt(const T x)
		{
			int jl = jsav, jm, ju, inc = 1;
			if ( n < 2 || m < 2 || m > n) throw("locate size error");
			bool ascnd = ((*XX)[n-1] >= (*XX)[0]) ;
			// Check if input guess is useful
			if ( jl < 0 || jl > n - 1)
			{
				jl = 0 ;
				ju = n - 1;
			}
			else
			{
				if ( x >= (*XX)[jl] == ascnd )
				{
					for(;;)
					{
						ju = jl + inc ;
						if( ju >= n-1 )
						{
							ju = n - 1 ;
							break ;
						}
						else if ( x < (*XX)[ju] == ascnd ) break ;
						else
						{
							jl = ju ;
							inc += inc ;
						}
					}
				}
				else
				{
					ju = jl;
					for(;;)
					{
						jl = jl - inc ;
						if (jl <= 0)
						{
							jl = 0 ;
							break ;
						}
						else if ( x >= (*XX)[jl] == ascnd) break;
						else
						{
							ju = jl;
							inc += inc;
						}
					}
				}
			}
			while ( ju - jl > 1)
			{
				jm = (ju + jl) / 2 ;
				if( x >= (*XX)[jm] == ascnd) jl = jm ;
				else ju = jm ;
			}
			cor = abs( jl -jsav ) > dj ? 0 : 1;
			jsav = jl ;
			return max(0, min(n-m, jl - (m-2)/2 ));
		}

		T interp( const T x )
		{
			// Given x, return interpolated value, using data pointed to xx and yy
			int jlo = cor ? hunt(x) : locate(x) ;
			return rawinterp(jlo, x) ;
		}

		virtual T rawinterp(int jlo, T x) = 0 ; //Method to be provided by derived classes

};

template<class T> 
class Spline_Interp : public Base_Interp<T> {
	private:
		vector<T> Y2 ; // Values of the 2nd derivatives
	public:
		Spline_Interp(const vector<T> &X, const vector<T> &Y, T yp1 = 1e99, T ypn = 1e99) : Base_Interp<T>(X, Y, 2)
		{
			Y2.resize(X.size()) ;
			sety2(X, Y, yp1, ypn); // Computes the values of y2
		}

		void sety2(const vector<T> &X, const vector<T> &Y, T yp1, T ypn)
		{
			int i, k ;
			T p, qn, sig, un ;
			int n = Y2.size() ;
			vector<T> U(n-1) ;
			if (yp1 > 0.99e99)
			{
				Y2[0] = U[0] = 0.0 ;
			}
			else
			{
				Y2[0] = - 0.5 ;
				U[0] = (3.0/(X[1]-X[0])) * ((Y[1]-Y[0])/(X[1]-X[0])-yp1);
			}
			for(i = 1; i < n - 1; i++)
			{
				sig = (X[i]-X[i-1]) / (X[i+1]-X[i-1]) ;
				p = sig * Y2[i-1] + 2.0 ;
				Y2[i] = (sig -1.0) / p ;
				U[i] = (Y[i+1]-Y[i]) / (X[i+1]-X[i]) - (Y[i]-Y[i-1]) / (X[i] - X[i-1]) ;
				U[i] = (6.0 * U[i] / (X[i+1]-X[i-1]) - sig * U[i-1]) / p ;
			}
			if (ypn > 0.99e99)
			{
				qn = un = 0.0 ;
			}
			else
			{
				qn = 0.5 ;
				un = (3.0 / (X[n-1] - X[n-2])) * (ypn - (Y[n-1] - Y[n-2]) / (X[n-1] - X[n-2])) ;
			}
			Y2[n-1] = (un - qn * U[n-2]) / (qn * Y2[n-2] + 1.0 );
			for (k = n - 2; k >= 0 ; k--)
			{
				Y2[k] = Y2[k] * Y2[k+1] + U[k];
			}
		}
		
		T rawinterp(int jl, T x)
		{
			int klo = jl, khi = jl + 1 ;
			T y, h, b, a ;
			h = (*Base_Interp<T>::XX)[khi] - (*Base_Interp<T>::XX)[klo] ;
			if ( h == 0.0) throw("Bad input to routine splint");
			a = ( (*Base_Interp<T>::XX)[khi] - x ) / h ;
			b = ( x - (*Base_Interp<T>::XX)[klo] ) / h ;
			y = a * (*Base_Interp<T>::YY)[klo] + b * (*Base_Interp<T>::YY)[khi] + ( (a*a*a - a) * Y2[klo] + (b*b*b - b) * Y2[khi] ) * ( h * h) / 6.0;
			return y ;
		}
};

#endif