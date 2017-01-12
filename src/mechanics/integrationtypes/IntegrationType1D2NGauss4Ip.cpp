// $Id$

#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif //ENABLE_VISUALIZE

#include "mechanics/integrationtypes/IntegrationType1D2NGauss4Ip.h"
#include <assert.h>


// constructor
NuTo::IntegrationType1D2NGauss4Ip::IntegrationType1D2NGauss4Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
Eigen::VectorXd NuTo::IntegrationType1D2NGauss4Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    switch (rIpNum)
    {
    case 0 : return Eigen::Matrix<double, 1, 1>::Constant(-0.861136311594052575224);
    case 1 : return Eigen::Matrix<double, 1, 1>::Constant(-0.339981043584856264803);
    case 2 : return Eigen::Matrix<double, 1, 1>::Constant( 0.339981043584856264803);
    case 3 : return Eigen::Matrix<double, 1, 1>::Constant( 0.861136311594052575224);
    default:
        throw MechanicsException("[NuTo::IntegrationType1D2NGauss4Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}

//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType1D2NGauss4Ip::GetNumIntegrationPoints()const
{
    return 4;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType1D2NGauss4Ip::GetIntegrationPointWeight(int rIpNum)const
{
    switch (rIpNum)
    {
    case 0 :
        return 0.34785484513745385737;
    case 1 :
    	return 0.65214515486254614262;
    case 2 :
    	return 0.65214515486254614262;
    case 3 :
        return 0.34785484513745385737;
    default:
        throw MechanicsException("[NuTo::IntegrationType1D2NGauss4Ip::GetIntegrationPointWeight] Ip number out of range.");
    }
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType1D2NGauss4Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 4;
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-0.3043);
    VisualizationPointLocalCoordinates.push_back(0.3043);
    VisualizationPointLocalCoordinates.push_back(1);
    NumVisualizationCells = 3;
    VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
    VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
    VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIP.push_back(0);
    VisualizationCellsIP.push_back(1);
    VisualizationCellsIP.push_back(2);
}
#endif // ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType1D2NGauss4Ip)
#endif
