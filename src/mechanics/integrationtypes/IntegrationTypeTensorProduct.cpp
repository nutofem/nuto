#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "math/Legendre.h"

//! @brief computes points and weights for Lobatto quadrature in 1D		
//! @param numIPs number of integration points		
//! @return pair of quadrature weights and points range [-1,1] including boundary points
std::pair<std::vector<double>, std::vector<double>> ComputeWeightsAndPoints1DLobatto(int nIps)
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

//! @brief computes points and weights for Gauss quadrature in 1D		
//! @param numIPs number of integration points		
//! @return pair of quadrature weights and points range (-1,1)
std::pair<std::vector<double>, std::vector<double>> ComputeWeightsAndPoints1DGauss(int nIps)
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

std::pair<std::vector<double>, std::vector<double>> ComputeWeightsAndPoints1D(size_t nIps, NuTo::eIntegrationMethod method)
{
    if (method == NuTo::eIntegrationMethod::GAUSS)
        return ComputeWeightsAndPoints1DGauss(nIps);
    if (method == NuTo::eIntegrationMethod::LOBATTO)
        return ComputeWeightsAndPoints1DLobatto(nIps);
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

#ifdef ENABLE_VISUALIZE
template <int TDim>
void NuTo::IntegrationTypeTensorProduct<TDim>::GetVisualizationPoints(
        unsigned int& NumVisualizationPoints, std::vector<double>& VisualizationPointLocalCoordinates) const
{
    std::vector<double> VisualizationPoints1D;
    size_t NumVisualizationPoints1D = mIPts1D.size() + 1;

    VisualizationPoints1D.push_back(-1.);
    for (size_t i = 1; i < NumVisualizationPoints1D - 1; i++)
    {
        VisualizationPoints1D.push_back(0.5 * (mIPts1D[i - 1] + mIPts1D[i]));
    }
    VisualizationPoints1D.push_back(1.);


    NumVisualizationPoints = 1;
    for (int dim = 0; dim < TDim; dim++)
        NumVisualizationPoints *= NumVisualizationPoints1D;

    VisualizationPointLocalCoordinates.reserve(NumVisualizationPoints);

    for (unsigned int i = 0; i < NumVisualizationPoints; i++)
    {
        int power = 1;
        for (int dim = 0; dim < TDim; dim++)
        {
            size_t index = (i / power) % NumVisualizationPoints1D;
            VisualizationPointLocalCoordinates.push_back(VisualizationPoints1D[index]);
            power *= NumVisualizationPoints1D;
        }
    }
}

namespace NuTo
{
template <>
void NuTo::IntegrationTypeTensorProduct<1>::GetVisualizationCells(
        unsigned int& NumVisualizationPoints, std::vector<double>& VisualizationPointLocalCoordinates,
        unsigned int& NumVisualizationCells, std::vector<NuTo::eCellTypes>& VisualizationCellType,
        std::vector<unsigned int>& VisualizationCellsIncidence, std::vector<unsigned int>& VisualizationCellsIP) const
{
    GetVisualizationPoints(NumVisualizationPoints, VisualizationPointLocalCoordinates);

    NumVisualizationCells = mIPts.size();
    for (size_t i = 0; i < NumVisualizationCells; i++)
    {
        VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
    }

    for (size_t i = 0; i < NumVisualizationCells; i++)
    {
        VisualizationCellsIncidence.push_back(i);
        VisualizationCellsIncidence.push_back(i + 1);
    }

    for (size_t i = 0; i < NumVisualizationCells; i++)
    {
        VisualizationCellsIP.push_back(i);
    }
}

template <>
void NuTo::IntegrationTypeTensorProduct<2>::GetVisualizationCells(
        unsigned int& NumVisualizationPoints, std::vector<double>& VisualizationPointLocalCoordinates,
        unsigned int& NumVisualizationCells, std::vector<NuTo::eCellTypes>& VisualizationCellType,
        std::vector<unsigned int>& VisualizationCellsIncidence, std::vector<unsigned int>& VisualizationCellsIP) const
{
    GetVisualizationPoints(NumVisualizationPoints, VisualizationPointLocalCoordinates);

    NumVisualizationCells = 0;

    size_t numIPs1D = mIPts1D.size();
    for (size_t row = 0; row < numIPs1D; row++)
    {
        for (size_t col = 0; col < numIPs1D; col++)
        {
            size_t start = row * (numIPs1D + 1) + col;

            VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);

            VisualizationCellsIncidence.push_back(start);
            VisualizationCellsIncidence.push_back(start + 1);
            VisualizationCellsIncidence.push_back(start + numIPs1D + 2); // the number of 1D cell points is numIPs1D+1
            VisualizationCellsIncidence.push_back(start + numIPs1D + 1);

            VisualizationCellsIP.push_back(NumVisualizationCells);
            NumVisualizationCells++;
        }
    }
}

template <>
void NuTo::IntegrationTypeTensorProduct<3>::GetVisualizationCells(
        unsigned int& NumVisualizationPoints, std::vector<double>& VisualizationPointLocalCoordinates,
        unsigned int& NumVisualizationCells, std::vector<NuTo::eCellTypes>& VisualizationCellType,
        std::vector<unsigned int>& VisualizationCellsIncidence, std::vector<unsigned int>& VisualizationCellsIP) const
{
    GetVisualizationPoints(NumVisualizationPoints, VisualizationPointLocalCoordinates);

    NumVisualizationCells = 0;

    size_t numIPs1D = mIPts1D.size();
    for (size_t height = 0; height < numIPs1D; height++)
    {
        for (size_t row = 0; row < numIPs1D; row++)
        {
            for (size_t col = 0; col < numIPs1D; col++)
            {
                size_t start = row * (numIPs1D + 1) + col + height * ((numIPs1D + 1) * (numIPs1D + 1));

                VisualizationCellsIncidence.push_back(start);
                VisualizationCellsIncidence.push_back(start + 1);
                VisualizationCellsIncidence.push_back(start + numIPs1D + 2); 
                VisualizationCellsIncidence.push_back(start + numIPs1D + 1);

                start = row * (numIPs1D + 1) + col + (height + 1) * ((numIPs1D + 1) * (numIPs1D + 1));

                VisualizationCellsIncidence.push_back(start);
                VisualizationCellsIncidence.push_back(start + 1);
                VisualizationCellsIncidence.push_back(start + numIPs1D + 2);
                VisualizationCellsIncidence.push_back(start + numIPs1D + 1);

                VisualizationCellType.push_back(NuTo::eCellTypes::HEXAHEDRON);
                VisualizationCellsIP.push_back(NumVisualizationCells);
                NumVisualizationCells++;
            }
        }
    }
}

} // namespace NuTo

#endif // ENABLE_VISUALIZE

template class NuTo::IntegrationTypeTensorProduct<1>;
template class NuTo::IntegrationTypeTensorProduct<2>;
template class NuTo::IntegrationTypeTensorProduct<3>;
