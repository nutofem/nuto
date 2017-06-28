#include "math/Polynomial.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

double NuTo::Math::Polynomial::Legendre(int n, double x) {
  if (n == 0)
    return 1.;
  if (n == 1)
    return x;
  return ( ((2. * n - 1.) * x * Legendre(n - 1, x) -
           (n - 1.) * Legendre(n - 2, x)) / n );
}

double NuTo::Math::Polynomial::LegendreDeriv(int n, double x, int k) {
  if (n < k)
    return 0.;
  if (n == k) {
    int fac_k = 1;
    int fac_2k = 1.;
    int pow2_k = 1.;
    for (int i = 1; i <= k; i++) {
      fac_k *= i;
      pow2_k *= 2;
    }
    for (int i = 1; i <= (2 * k); i++) {
      fac_2k *= i;
    }
    return (fac_2k / (pow2_k * fac_k));
  }
  return (((2. * n - 1.) * x * LegendreDeriv(n - 1, x, k) -
           (n - 1. + k) * LegendreDeriv(n - 2, x, k)) /
          (n - k));
}

std::vector<double> NuTo::Math::Polynomial::LegendreRoots(int n) {
  // first guess
  std::vector<double> xs(n);
  for (int i = 1; i <= n; i++) {
    xs[i-1] = (cos(M_PI * (4. * i - 1) / (4. * n + 2)));
  }
  // newton iterations
  double dx;
  double tol = 1e-15;
  for (double& x : xs) {
    double xnew;
    do {
      xnew = x - Legendre(n, x) / (LegendreDeriv(n, x, 1));
      dx = std::abs(x - xnew);
      x = xnew;
    } while (dx > tol);
  }
  std::sort(xs.begin(), xs.end());
  return (xs);
}

std::vector<double> NuTo::Math::Polynomial::LegendreDerivRoots(int n) {
  // first guess
  std::vector<double> ys(n);
  for (int i = 1; i <= n; i++) {
    ys[i-1] = (cos(M_PI * (4. * i - 1) / (4. * n + 2)));
  }
  std::vector<double> xs(n - 1);
  for (int i = 1; i < n; i++) {
    xs[i-1] = (0.5 * (ys[i - 1] + ys[i]));
  }
  // newton iterations
  double dx;
  double tol = 1e-15;
  for (double& x : xs) {
    double xnew;
    do {
      xnew = x - LegendreDeriv(n, x, 1) / (LegendreDeriv(n, x, 2));
      dx = std::abs(x - xnew);
      x = xnew;
    } while (dx > tol);
  }
  std::sort(xs.begin(), xs.end());
  return (xs);
}
