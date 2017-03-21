#include "mechanics/integrationtypes/IntegrationType1D2NGauss2Ip.h"
#include <cassert>

#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

//! @brief constructor
NuTo::IntegrationType1D2NGauss2Ip::IntegrationType1D2NGauss2Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
Eigen::VectorXd NuTo::IntegrationType1D2NGauss2Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    assert(rIpNum==0 || rIpNum==1);
    switch (rIpNum)
    {
    case 0 : return Eigen::Matrix<double, 1, 1>::Constant(-0.577350269189626);
    case 1 : return Eigen::Matrix<double, 1, 1>::Constant( 0.577350269189626);
    default:
        throw MechanicsException("[NuTo::IntegrationType1D2NGauss2Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType1D2NGauss2Ip::GetNumIntegrationPoints()const
{
    return 2;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType1D2NGauss2Ip::GetIntegrationPointWeight(int rIpNum)const
{
    return 1;
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType1D2NGauss2Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 3;
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(1);
    NumVisualizationCells = 2;
    VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
    VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIP.push_back(0);
    VisualizationCellsIP.push_back(1);
}
#endif // ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType1D2NGauss2Ip)
#endif
