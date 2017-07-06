#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

#include "mechanics/integrationtypes/IntegrationTypeTensorProductGauss.h"
#include "mechanics/elements/ElementShapeFunctions.h"
#include "math/Legendre.h"
#include <iostream>

// constructor
template<int TDim>
NuTo::IntegrationTypeTensorProductGauss<TDim>::IntegrationTypeTensorProductGauss(size_t nIps)
{
    mIPts1D = NuTo::Math::Polynomial::LegendreRoots(nIps);
    std::vector<double> weights(nIps);

    if (nIps > 0)
    {
        for (size_t i = 0; i < nIps; i++)
        {
            double dl = NuTo::Math::Polynomial::Legendre(nIps, mIPts1D[i], 1);
            weights.push_back(2. / (1. - mIPts1D[i] * mIPts1D[i]) / (dl * dl));
        }
    }
    else throw MechanicsException(__PRETTY_FUNCTION__, "Ip number out of range.");

    int numIpsDim = 1;
    for(int dim = 0; dim < TDim; dim++) numIpsDim *= nIps;

    // calculate the tensor product integratration point coordinates and weights
    mWeights.resize(numIpsDim);
    mIPts.resize(numIpsDim);

    for(int i = 0; i < numIpsDim; i++)
    {
        mWeights[i] = 0.;
        int power = 1;
        for(int dim = 0; dim < TDim; dim++)
        {
            int index = (i/power) % nIps;
            mWeights[i] *= weights[index];
            mIPts[i][dim] = mIPts1D[index];
            power *= nIps;
        }
    }
}

template<int TDim>
Eigen::VectorXd NuTo::IntegrationTypeTensorProductGauss<TDim>::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    if (rIpNum >= 0 && rIpNum < (int)mIPts.size())
        return mIPts[rIpNum];
    else
        throw MechanicsException(
                "[NuTo::IntegrationTypeTensorProductGauss::GetLocalIntegrationPointCoordinates] Ip number out of range.");
}

template<int TDim>
int NuTo::IntegrationTypeTensorProductGauss<TDim>::GetNumIntegrationPoints() const
{
    return mIPts.size();
}

template<int TDim>
double NuTo::IntegrationTypeTensorProductGauss<TDim>::GetIntegrationPointWeight(int rIpNum) const
{
    if (rIpNum >= 0 && rIpNum < (int)mIPts.size())
        return mWeights[rIpNum];
    throw MechanicsException("[NuTo::IntegrationTypeTensorProductGauss::GetIntegrationPointWeight] Ip number out of range.");
}

#ifdef ENABLE_VISUALIZE
template<int TDim>
void NuTo::IntegrationTypeTensorProductGauss<TDim>::GetVisualizationPoints(unsigned int& NumVisualizationPoints, std::vector<double>& VisualizationPointLocalCoordinates) const
{
    std::vector<double> VisualizationPoints1D;
    size_t NumVisualizationPoints1D = mIPts1D.size()+1;

    VisualizationPoints1D.push_back(-1.);
    for (size_t i = 1; i < NumVisualizationPoints1D - 1 ; i++)
    {
        VisualizationPoints1D.push_back(0.5 * (mIPts1D[i - 1] + mIPts1D[i]));
    }
    VisualizationPoints1D.push_back(1.);


    NumVisualizationPoints = 1;
    for(int dim = 0; dim < TDim; dim++) NumVisualizationPoints *= NumVisualizationPoints1D;

    VisualizationPointLocalCoordinates.resize(NumVisualizationPoints);

    for(unsigned int i = 0; i < NumVisualizationPoints; i++)
    {
        int power = 1;
        for(int dim = 0; dim < TDim; dim++)
        {
            size_t index = (i/power) % NumVisualizationPoints1D;
            VisualizationPointLocalCoordinates.push_back(VisualizationPoints1D[index]);
            power *= NumVisualizationPoints1D;
        }
        std::cout << std::endl;
    }
}

namespace NuTo
{
template <>
void NuTo::IntegrationTypeTensorProductGauss<1>::GetVisualizationCells(unsigned int& NumVisualizationPoints,
                                                                       std::vector<double>& VisualizationPointLocalCoordinates,
                                                                       unsigned int& NumVisualizationCells,
                                                                       std::vector<NuTo::eCellTypes>& VisualizationCellType,
                                                                       std::vector<unsigned int>& VisualizationCellsIncidence,
                                                                       std::vector<unsigned int>& VisualizationCellsIP) const
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
void NuTo::IntegrationTypeTensorProductGauss<2>::GetVisualizationCells(unsigned int& NumVisualizationPoints,
                                                                       std::vector<double>& VisualizationPointLocalCoordinates,
                                                                       unsigned int& NumVisualizationCells,
                                                                       std::vector<NuTo::eCellTypes>& VisualizationCellType,
                                                                       std::vector<unsigned int>& VisualizationCellsIncidence,
                                                                       std::vector<unsigned int>& VisualizationCellsIP) const
{
    GetVisualizationPoints(NumVisualizationPoints, VisualizationPointLocalCoordinates);

    NumVisualizationCells = 0;

    size_t numIPs1D = mIPts1D.size();
    for(size_t row = 0; row < numIPs1D; row++)
    {
        for(size_t col = 0; col < numIPs1D; col++)
        {
            size_t start = row*(numIPs1D+1) + col;

            VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);

            VisualizationCellsIncidence.push_back(start);
            VisualizationCellsIncidence.push_back(start + 1);
            VisualizationCellsIncidence.push_back(start + numIPs1D + 2); // the number of 1D cell points is numIPs1D+1
            VisualizationCellsIncidence.push_back(start + numIPs1D + 1 );

            VisualizationCellsIP.push_back(NumVisualizationCells);
            NumVisualizationCells++;
        }
    }
}

template <>
void NuTo::IntegrationTypeTensorProductGauss<3>::GetVisualizationCells(unsigned int& NumVisualizationPoints,
                                                                       std::vector<double>& VisualizationPointLocalCoordinates,
                                                                       unsigned int& NumVisualizationCells,
                                                                       std::vector<NuTo::eCellTypes>& VisualizationCellType,
                                                                       std::vector<unsigned int>& VisualizationCellsIncidence,
                                                                       std::vector<unsigned int>& VisualizationCellsIP) const
{
    GetVisualizationPoints(NumVisualizationPoints, VisualizationPointLocalCoordinates);

    NumVisualizationCells = 0;

    size_t numIPs1D = mIPts1D.size();
    for(size_t height = 0; height < numIPs1D; height++)
    {
        for(size_t row = 0; row < numIPs1D; row++)
        {
            for(size_t col = 0; col < numIPs1D; col++)
            {
                size_t start = row*(numIPs1D + 1) + col + height*((numIPs1D + 1)*(numIPs1D + 1));

                VisualizationCellsIncidence.push_back(start);
                VisualizationCellsIncidence.push_back(start + 1);
                VisualizationCellsIncidence.push_back(start + numIPs1D + 2); // the number of 1D cell points is numIPs1D+1
                VisualizationCellsIncidence.push_back(start + numIPs1D + 1);

                start = row*(numIPs1D + 1) + col + (height+1)*((numIPs1D + 1)*(numIPs1D + 1));

                VisualizationCellsIncidence.push_back(start);
                VisualizationCellsIncidence.push_back(start + 1);
                VisualizationCellsIncidence.push_back(start + numIPs1D + 2);
                VisualizationCellsIncidence.push_back(start + numIPs1D + 1);

                VisualizationCellType.push_back(NuTo::eCellTypes::HEXAHEDRON);
                VisualizationCellsIP.push_back(NumVisualizationCells);
                NumVisualizationCells++;

                std::cout << std::endl;
            }
        }
    }
}

} // namespace NuTo

#endif // ENABLE_VISUALIZE

template class NuTo::IntegrationTypeTensorProductGauss<1>;
template class NuTo::IntegrationTypeTensorProductGauss<2>;
template class NuTo::IntegrationTypeTensorProductGauss<3>;
