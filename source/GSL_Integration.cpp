#include <iostream>
#include <string>
#include "GSL_Math.hpp"
#include "Precision.hpp"
#include "Error.hpp"

using namespace std;

namespace GSL_Math{

//  a VERY simple numerical integrator
    double Simple_Integrator( double (*f)(double, void *), double x1, double x2, void *param){
        double result=0;
        int N = 500;
        double ftemp[N];
        double dx = (x2-x1)/(N-1.);
        ftemp[0] = f(x1, param);

        for(int i=1; i<N; ++i){ //  to reduce evaluations of f()
            ftemp[i] = f(x1+i*dx, param);
            result += 0.5*(ftemp[i-1] + ftemp[i])*dx;
        }
		return result;
    }

//	GSL qags integrator
    double GSL_Integrator( double (*f)(double, void *), double x1, double x2, void *param ){
        double result, error;
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(GSL_QAGS_ITER);
        gsl_function F;
        F.function  = f;
        F.params    = param;
        gsl_integration_qags( &F, x1, x2, GSL_QAGS_ABSTOT, GSL_QAGS_RELTOT, GSL_QAGS_ITER, w, &result, &error );
        gsl_integration_workspace_free(w);
        return result;
    }

    double GSL_Integrator( double (*f)(double, void *), double x1, double x2, double abs_eps, double rel_eps, void *param ){
        double result, error;
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(GSL_QAGS_ITER);
        gsl_function F;
        F.function  = f;
        F.params    = param;
        gsl_integration_qags( &F, x1, x2, abs_eps, rel_eps, GSL_QAGS_ITER, w, &result, &error );
        gsl_integration_workspace_free(w);
        return result;
    }

//	added @ 4-16-2015
	double Romberg_Integrator( 	double (*f)(double, void *),
								double a, double b,
								void *param,
								double eps ){
		int m,n,i,j;
		double y[_ROMBERG_ITER_MAX_];
		double h,ep, p, x, s, q=-9999;

		h = b-a;
		y[0] = h*( (*f)(a,param) + (*f)(b,param) )/2.0;
		m = 1;
		n = 1;
		ep = eps + 1.0;

		while( (ep >= eps) && (m < _ROMBERG_ITER_MAX_) ){
			p = 0.0;
			for( i=0; i<=n-1; ++i ){
				x = a + (i+0.5)*h;
				p = p + (*f)(x,param);
			}

			p = (y[0] + h*p)/2.0;
			s = 1.0;

			for( j=1; j<=m; ++j ){
				s = 4.0*s;
				q = (s*p - y[j-1])/(s-1.0);
				y[j-1] = p;
				p = q;
			}

			ep = fabs(q-y[m-1]);
			m = m + 1;
			y[m-1] = q;
			n = n + n;
			h = h/2.0;
		}

		//	one possible problem is that even when the iteration number has exceeded _ROMBERG_ITER_MAX_,
		//	the precision requirement still has not been satisified
		if( (m >= _ROMBERG_ITER_MAX_) && (ep >= eps)){
			print_icosmo_warning;
			std::string err = "\n";
            err += "*** double Romberg_Integrator( double (*f)(double, void *), double a, double b, void *param, double eps ) ==>\n";
			err += "*** the iteration has exceeded the maximum number, but the presion still has not reached!";
		}

		return q;
	}

	double Romberg_Integrator( double xi[], double yi[], int size, double a, double b, double eps ){
	//	use cubic spline
		gsl_interp_accel 	*acc 	= gsl_interp_accel_alloc();
		gsl_spline			*spline	= gsl_spline_alloc(gsl_interp_cspline, size);
		gsl_spline_init(spline, xi, yi, size);

		int m,n,i,j;
		double y[_ROMBERG_ITER_MAX_];
		double h,ep, p, x, s, q=-9999;

		//	check integration range
		if( a < xi[0] ){
			print_icosmo_error;
			std::string err = "\n";
            err += "*** Romberg_Integration( double xi[], double yi[], int size, double a, double b, double eps )\n";
			err += "*** a should be >= xi[0]";
			throw std::runtime_error(err);
		}

		if( b > xi[size-1] ){
			print_icosmo_error;
			std::string err = "\n";
            err += "*** Romberg_Integration( double xi[], double yi[], int size, double a, double b, double eps )\n";
			err += "*** b should be <= xi[size-1]";
			throw std::runtime_error(err);
		}

		h = b - a;

		if( h < 0 ){
			std::string err = "\n";
            err += "double Romberg_Integrator( double xi[], double yi[], int size, double xx, double eps )\n";
			err += "==> xx must be greater than xi[0], so that h > 0!";
			throw std::runtime_error(err);
		}

		y[0] = h*( yi[0] + yi[size-1] )/2.0;
		m = n = 1;
		ep = eps + 1.0;

		while( (ep >= eps) && (m < _ROMBERG_ITER_MAX_) ){
			p = 0.0;
			for( i=0; i<n; ++i ){
				x = a + (i+0.5)*h;
				p = p + gsl_spline_eval(spline, x, acc);
			}

			p = (y[0] + h*p)/2.0;
			s = 1.0;

			for( j=1; j<=m; ++j ){
				s = 4.0*s;
				q = (s*p - y[j-1])/(s-1.0);
				y[j-1] = p;
				p = q;
			}

			ep = fabs(q-y[m-1]);
			++m;
			y[m-1] = q;
			n += n;
			h = h/2.0;
		}

		gsl_spline_free(spline);
		gsl_interp_accel_free(acc);

		//	one possible problem is that even when the iteration number has exceeded _ROMBERG_ITER_MAX_,
		//	the precision requirement still has not been satisified
		if( (m >= _ROMBERG_ITER_MAX_) && (ep >= eps)){
			print_icosmo_warning;
			std::string err = "\n";
            err += "*** double Romberg_Integrator( double (*f)(double, void *), double a, double b, void *param, double eps ) ==>\n";
			err += "*** the iteration has exceeded the maximum number, but the presion still has not reached!";
            throw std::runtime_error(err);
		}

		return q;
	}

}
