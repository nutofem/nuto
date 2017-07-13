#include "math/Legendre.h"
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/algorithm/transform.hpp>
#include <cmath>
#include <vector>


double NuTo::Math::Polynomial::Legendre(int n, double x, int k)
{
    auto Factorial = [=](int i) -> int { return std::tgamma(i + 1); };

    if (n < k)
        return 0.;
    if (n == k)
        return Factorial(2 * k) / (std::pow(2, k) * Factorial(k));
    return (((2. * n - 1.) * x * Legendre(n - 1, x, k) - (n - 1. + k) * Legendre(n - 2, x, k)) / (n - k));
}

double FindLegendreRoot(double guess, int n, int derivative)
{
    using namespace NuTo::Math::Polynomial;
    const double tol = 1e-15;
    double x = guess;
    while (true)
    {
        double xnew = x - Legendre(n, x, derivative) / (Legendre(n, x, derivative + 1));
        if (std::abs(x - xnew) < tol)
            return xnew;
        x = xnew;
    }
}

std::vector<double> FindLegendreRoots(std::vector<double> guess, int n, int derivative)
{
    auto findRoot = [=](double guess) { return FindLegendreRoot(guess, n, derivative); };
    boost::range::transform(guess, guess.begin(), findRoot);
    return boost::range::sort(guess);
}

std::vector<double> FirstGuess(int n)
{
    std::vector<double> xs(n);
    for (int i = 1; i <= n; i++)
        xs[i - 1] = (cos(M_PI * (4. * i - 1) / (4. * n + 2)));
    return xs;
}

std::vector<double> NuTo::Math::Polynomial::LegendreRoots(int n)
{
    return FindLegendreRoots(FirstGuess(n), n, 0);
}

std::vector<double> NuTo::Math::Polynomial::LegendreDerivRoots(int n)
{
    std::vector<double> ys = FirstGuess(n);
    std::vector<double> xs(n - 1);
    for (int i = 1; i < n; i++)
    {
        xs[i - 1] = (0.5 * (ys[i - 1] + ys[i]));
    }
    return FindLegendreRoots(xs, n, 1);
}
