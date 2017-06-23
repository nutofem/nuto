#include <stdlib.h>
#include <iostream>
#include <vector>
#include <math.h>
#include "boost/math/special_functions/legendre.hpp"
#include <algorithm>
#include "math/Polynomial.h"

double NuTo::Math::Polynomial::legendre_p_deriv(int n, double x, int k) {
    if (n < k) {
        return(0.);
    } else if (n==k) {
        double fac_k = 1.;
        double fac_2k = 1.;
        double pow2_k = 1.;
        for (int i=1; i<=k;i++) {
            fac_k *= i;
            pow2_k *= 2;
        }
        for (int i=1; i<= (2*k);i++) {
            fac_2k *= i;
        }
        return( fac_2k/(pow2_k * fac_k) );
    } else {
        return( ( (2.*n-1.)*x* legendre_p_deriv(n-1,x, k) - (n-1.+k) * legendre_p_deriv(n-2,x, k) )/(n-k) );
    }
}

std::vector<double> NuTo::Math::Polynomial::legendre_roots(int n) {
    // first guess
    std::vector<double> x;
    for(int i=1; i<=n; i++) {
        x.push_back(cos(M_PI*(4.*i-1)/(4.*n + 2)));
    }
    // newton iterations
    double dx;
    double tol = 1e-16;
    for(int i=0; i<n; i++) 
    {
        double xnew;
        do
        {
            xnew = x[i] - boost::math::legendre_p(n,x[i])/(legendre_p_deriv(n,x[i],1));
            dx = abs(x[i] - xnew);
            x[i] = xnew;
        }  
        while (dx > tol);
    }
    std::sort(x.begin(),x.end());
    return(x);
}

std::vector<double> NuTo::Math::Polynomial::legendre_deriv_roots(int n) {
    // first guess
    std::vector<double> y;
    for(int i=1; i<=n; i++) {
        y.push_back(cos(M_PI*(4.*i-1)/(4.*n + 2)));
    }
    std::vector<double> x;
    for(int i=1; i<n; i++) {
        x.push_back(0.5 * (y[i-1] + y[i]) );
    }
    // newton iterations
    double dx;
    double tol = 1e-8;
    for(int i=0; i<(n-1); i++)
    {
        double xnew;
        do
        {
            xnew = x[i] - legendre_p_deriv(n,x[i],1)/(legendre_p_deriv(n,x[i],2));
            dx = abs(x[i] - xnew);
            x[i] = xnew;
        }
        while (dx > tol);
    }
    std::sort(x.begin(),x.end());
    return(x);
}
