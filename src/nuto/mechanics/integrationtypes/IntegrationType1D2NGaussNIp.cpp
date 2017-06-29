// $Id$

#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGaussNIp.h"
#include <assert.h>
#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

// constructor
template<int TnumIPs>
NuTo::IntegrationType1D2NGaussNIp<TnumIPs>::IntegrationType1D2NGaussNIp()
{
    mCoordinates.setZero(TnumIPs);
    mWeights.setZero(TnumIPs);
    gauleg(mCoordinates, mWeights);
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)#
template<int TnumIPs>
void NuTo::IntegrationType1D2NGaussNIp<TnumIPs>::GetLocalIntegrationPointCoordinates1D(int rIpNum, double& rCoordinates)const
{
    assert(rIpNum < TnumIPs);

    rCoordinates = mCoordinates(rIpNum);
}

//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
template<int TnumIPs>
int NuTo::IntegrationType1D2NGaussNIp<TnumIPs>::GetNumIntegrationPoints()const
{
    return TnumIPs;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
template<int TnumIPs>
double NuTo::IntegrationType1D2NGaussNIp<TnumIPs>::GetIntegrationPointWeight(int rIpNum)const
{
    assert(rIpNum < TnumIPs);
    return mWeights(rIpNum);
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
template<int TnumIPs>
std::string NuTo::IntegrationType1D2NGaussNIp<TnumIPs>::GetStrIdentifier()const
{
    return GetStrIdentifierStatic();
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
template<int TnumIPs>
std::string NuTo::IntegrationType1D2NGaussNIp<TnumIPs>::GetStrIdentifierStatic()
{
    std::string str = "1D2NGAUSS";
    str.append(std::to_string(TnumIPs));
    str.append("IP");
    return str;
}

#ifdef ENABLE_VISUALIZE
template<int TnumIPs>
void NuTo::IntegrationType1D2NGaussNIp<TnumIPs>::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    double coordinates[2];
    NumVisualizationPoints = GetNumIntegrationPoints()+1;
    VisualizationPointLocalCoordinates.push_back(-1.);

    for(int i = 1; i < GetNumIntegrationPoints(); i++)
    {
        GetLocalIntegrationPointCoordinates1D(i-1, coordinates[0]);
        GetLocalIntegrationPointCoordinates1D(i, coordinates[1]);
        VisualizationPointLocalCoordinates.push_back((coordinates[0] + coordinates[1])/2);
    }
    VisualizationPointLocalCoordinates.push_back(1.);

    NumVisualizationCells = GetNumIntegrationPoints();
    for(int i = 0; i < GetNumIntegrationPoints(); i++)  VisualizationCellType.push_back(NuTo::eCellTypes::LINE);

    VisualizationCellsIncidence.push_back(0);
    for(int i = 1; i < GetNumIntegrationPoints(); i++)
        for(int j = 0; j < 2; j++) VisualizationCellsIncidence.push_back(i);
    VisualizationCellsIncidence.push_back(TnumIPs);

    for(int i = 0; i < GetNumIntegrationPoints(); i++) VisualizationCellsIP.push_back(i);
}
#endif // ENABLE_VISUALIZE

template class NuTo::IntegrationType1D2NGaussNIp<6>;
template class NuTo::IntegrationType1D2NGaussNIp<7>;
template class NuTo::IntegrationType1D2NGaussNIp<8>;
template class NuTo::IntegrationType1D2NGaussNIp<9>;
template class NuTo::IntegrationType1D2NGaussNIp<10>;
template class NuTo::IntegrationType1D2NGaussNIp<11>;
template class NuTo::IntegrationType1D2NGaussNIp<12>;
template class NuTo::IntegrationType1D2NGaussNIp<13>;
template class NuTo::IntegrationType1D2NGaussNIp<14>;
template class NuTo::IntegrationType1D2NGaussNIp<15>;
template class NuTo::IntegrationType1D2NGaussNIp<16>;
template class NuTo::IntegrationType1D2NGaussNIp<17>;
template class NuTo::IntegrationType1D2NGaussNIp<18>;
template class NuTo::IntegrationType1D2NGaussNIp<19>;
template class NuTo::IntegrationType1D2NGaussNIp<20>;


#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType1D2NGaussNIp)
#endif
