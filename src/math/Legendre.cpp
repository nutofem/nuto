#include "math/Legendre.h"
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/algorithm/transform.hpp>
#include <cmath>
#include <vector>


double NuTo::Math::Polynomial::Legendre(int n, double x, int k)
{
    auto Factorial = [=](int i) -> int { return std::tgamma(i + 1); };

    std::vector<double> vals(n + 1, 0.);
    if (n >= k)
        vals[k] = Factorial(2 * k) / (std::pow(2, k) * Factorial(k));
    for (int i = k + 1; i <= n; i++)
    {
        vals[i] = ((2. * i - 1.) * x * vals[i - 1] - (i - 1. + k) * vals[i - 2]) / (i - k);
    }
    return vals.back();
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

std::pair<std::vector<double>, std::vector<double>> NuTo::Math::Polynomial::ComputeWeightsAndPoints1DLobatto(int nIps)
{
    if (nIps <= 1)
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Ip number out of range.");

    std::vector<double> points = NuTo::Math::Polynomial::LegendreDerivRoots(nIps - 1);
    points.insert(points.begin(), -1.);
    points.push_back(1.);

    std::vector<double> weights;
    weights.push_back(2. / (nIps * (nIps - 1)));
    for (int i = 1; i < nIps - 1; i++)
    {
        double lp = NuTo::Math::Polynomial::Legendre(nIps - 1, points[i]);
        weights.push_back(2. / (nIps * (nIps - 1)) / (lp * lp));
    }
    weights.push_back(2. / (nIps * (nIps - 1)));
    return std::make_pair(weights, points);
}

std::pair<std::vector<double>, std::vector<double>> NuTo::Math::Polynomial::ComputeWeightsAndPoints1DGauss(int nIps)
{
    if (nIps <= 0)
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Ip number out of range.");

    std::vector<double> points = NuTo::Math::Polynomial::LegendreRoots(nIps);
    std::vector<double> weights;
    for (int i = 0; i < nIps; i++)
    {
        double dl = NuTo::Math::Polynomial::Legendre(nIps, points[i], 1);
        weights.push_back(2. / (1. - points[i] * points[i]) / (dl * dl));
    }
    return std::make_pair(weights, points);
}
