#include "math/Quadrature.h"
#include "math/Legendre.h"
#include <vector>
#include "base/Exception.h"

std::pair<std::vector<double>, std::vector<double>> NuTo::Math::ComputeWeightsAndPoints1DLobatto(int nIps)
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

std::pair<std::vector<double>, std::vector<double>> NuTo::Math::ComputeWeightsAndPoints1DGauss(int nIps)
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
