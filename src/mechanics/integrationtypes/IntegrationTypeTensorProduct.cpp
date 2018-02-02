#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "math/Legendre.h"

std::pair<std::vector<double>, std::vector<double>> ComputeWeightsAndPoints1D(size_t nIps,
                                                                              NuTo::eIntegrationMethod method)
{
    if (method == NuTo::eIntegrationMethod::GAUSS)
        return NuTo::Math::Polynomial::ComputeWeightsAndPoints1DGauss(nIps);
    if (method == NuTo::eIntegrationMethod::LOBATTO)
        return NuTo::Math::Polynomial::ComputeWeightsAndPoints1DLobatto(nIps);
    throw NuTo::Exception(__PRETTY_FUNCTION__, "Integration method not supported.");
}

template <int TDim>
NuTo::IntegrationTypeTensorProduct<TDim>::IntegrationTypeTensorProduct(size_t nIps, NuTo::eIntegrationMethod method)
{
    std::pair<std::vector<double>, std::vector<double>> weightsAndPoints1D = ComputeWeightsAndPoints1D(nIps, method);

    std::vector<double> weights1D = weightsAndPoints1D.first;
    mIPts1D = weightsAndPoints1D.second;

    int numIpsDim = std::pow(nIps, TDim);

    // calculate the tensor product integration point coordinates and weights
    mWeights.resize(numIpsDim, 1.); // initialize with 1
    mIPts.resize(numIpsDim);

    for (int i = 0; i < numIpsDim; i++)
    {
        int power = 1;
        for (int dim = 0; dim < TDim; dim++)
        {
            int index = (i / power) % nIps;
            mWeights[i] *= weights1D[index];
            mIPts[i][dim] = mIPts1D[index];
            power *= nIps;
        }
    }
}

template <int TDim>
Eigen::VectorXd NuTo::IntegrationTypeTensorProduct<TDim>::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    if (rIpNum >= 0 && rIpNum < (int)mIPts.size())
        return mIPts[rIpNum];
    throw Exception(__PRETTY_FUNCTION__, "Ip number out of range.");
}

template <int TDim>
int NuTo::IntegrationTypeTensorProduct<TDim>::GetNumIntegrationPoints() const
{
    return mIPts.size();
}

template <int TDim>
double NuTo::IntegrationTypeTensorProduct<TDim>::GetIntegrationPointWeight(int rIpNum) const
{
    if (rIpNum >= 0 && rIpNum < (int)mIPts.size())
        return mWeights[rIpNum];
    throw Exception(__PRETTY_FUNCTION__, "Ip number out of range.");
}

template class NuTo::IntegrationTypeTensorProduct<1>;
template class NuTo::IntegrationTypeTensorProduct<2>;
template class NuTo::IntegrationTypeTensorProduct<3>;
