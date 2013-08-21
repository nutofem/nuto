// $Id$

#include "nuto/mechanics/integrationtypes/IntegrationType3D4NGauss1Ip.h"
#include <assert.h>


//! @brief constructor
NuTo::IntegrationType3D4NGauss1Ip::IntegrationType3D4NGauss1Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
void NuTo::IntegrationType3D4NGauss1Ip::GetLocalIntegrationPointCoordinates3D(int rIpNum, double rCoordinates[3])const
{
    assert(rIpNum==0);
    rCoordinates[0] = 0.25;
    rCoordinates[1] = 0.25;
    rCoordinates[2] = 0.25;
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType3D4NGauss1Ip::GetNumIntegrationPoints()const
{
    return 1;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType3D4NGauss1Ip::GetIntegrationPointWeight(int rIpNum)const
{
    return 1/6.;
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType3D4NGauss1Ip::GetStrIdentifier()const
{
    return GetStrIdentifierStatic();
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType3D4NGauss1Ip::GetStrIdentifierStatic()
{
    return std::string("3D4NGAUSS1IP");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType3D4NGauss1Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::CellBase::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 4;

    // Point 0
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 1
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 2
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 3
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(1);

    NumVisualizationCells = 1;

    // cell 0
    VisualizationCellType.push_back(NuTo::CellBase::TETRAEDER);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIP.push_back(0);
}
#endif // ENABLE_VISUALIZE
