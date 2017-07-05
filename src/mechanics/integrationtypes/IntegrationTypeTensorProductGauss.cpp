#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

#include "mechanics/integrationtypes/IntegrationTypeTensorProductGauss.h"
#include "mechanics/elements/ElementShapeFunctions.h"
#include "math/Legendre.h"

// constructor
template<int TDim>
NuTo::IntegrationTypeTensorProductGauss<TDim>::IntegrationTypeTensorProductGauss(int nIps)
{
    mNumIPs1D = nIps;
    std::vector<double> IPts = NuTo::Math::Polynomial::LegendreRoots(nIps);
    std::vector<double> weights(nIps);

    if (nIps > 0)
    {
        for (int i = 0; i < nIps; i++)
        {
            double dl = NuTo::Math::Polynomial::Legendre(nIps, IPts[i], 1);
            weights.push_back(2. / (1. - IPts[i] * IPts[i]) / (dl * dl));
        }
    }
    else throw MechanicsException(__PRETTY_FUNCTION__, "Ip number out of range.");

    int numIpsDim = std::pow(nIps, TDim);
    // calculate the tensor product integratration point coordinates and weights
    mWeights.resize(numIpsDim);
    mIPts.resize(numIpsDim);

    for(int i = 0; i < numIpsDim; i++)
    {
        mWeights[i] = 0.;
        int power = 1;
        for(int dim = 0; dim < TDim; dim++)
        {
            int index = i/power;
            mWeights[i] *= weights[index];
            mIPts[i][dim] = IPts[index];
            power *= TDim;
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
    NumVisualizationPoints = mIPts.size() + 1;
    VisualizationPointLocalCoordinates.push_back(-1.);
    for (unsigned int i = 1; i < NumVisualizationPoints - 1; i++)
    {
        VisualizationPointLocalCoordinates.push_back(0.5 * (mIPts[i - 1](0,0) + mIPts[i](0,0)));
    }
    VisualizationPointLocalCoordinates.push_back(1.);

    NumVisualizationCells = mIPts.size();
    for (unsigned int i = 0; i < NumVisualizationCells; i++)
    {
        VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
    }

    for (unsigned int i = 0; i < NumVisualizationCells; i++)
    {
        VisualizationCellsIncidence.push_back(i);
        VisualizationCellsIncidence.push_back(i + 1);
    }

    for (unsigned int i = 0; i < NumVisualizationCells; i++)
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
//    std::vector<double> VisualizationPoints1D;
//    VisualizationPoints1D.resize(mNumIPs1D+1);

//    VisualizationPoints1D.push_back(-1.);
//    for (unsigned int i = 1; i < NumVisualizationPoints - 1; i++)
//    {
//        VisualizationPoints1D.push_back(0.5 * (mIPts[i - 1](0,0) + mIPts[i](0,0)));
//    }
//    VisualizationPoints1D.push_back(1.);

//    int numCoordinates = std::pow(mNumIPs1D+1, TDim);

//    VisualizationPointLocalCoordinates.resize(numCoordinates);

//    for(int i = 0; i < numCoordinates; i++)
//    {
//        int power = 1;
//        for(int dim = 0; dim < TDim; dim++)
//        {
//            int index = i/power;
//            VisualizationPointLocalCoordinates.push_back(VisualizationPoints1D[index]);
//            power *= TDim;
//        }
//    }

//    NumVisualizationCells =

//    for()
//    {
//        VisualizationCellType.push_back(NuTo::eCellTypes::QUAD); NuTo::eCellTypes::HEXAHEDRON

//    }
}

template <>
void NuTo::IntegrationTypeTensorProductGauss<3>::GetVisualizationCells(unsigned int& NumVisualizationPoints,
                                                                       std::vector<double>& VisualizationPointLocalCoordinates,
                                                                       unsigned int& NumVisualizationCells,
                                                                       std::vector<NuTo::eCellTypes>& VisualizationCellType,
                                                                       std::vector<unsigned int>& VisualizationCellsIncidence,
                                                                       std::vector<unsigned int>& VisualizationCellsIP) const
{

}

} // namespace NuTo

#endif // ENABLE_VISUALIZE

template class NuTo::IntegrationTypeTensorProductGauss<1>;
template class NuTo::IntegrationTypeTensorProductGauss<2>;
template class NuTo::IntegrationTypeTensorProductGauss<3>;
