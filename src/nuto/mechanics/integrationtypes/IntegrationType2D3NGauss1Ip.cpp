// $Id: IntegrationType2D3NGauss1Ip.cpp 276 2010-06-30 13:04:32Z unger3 $

#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss1Ip.h"
#include <assert.h>


//! @brief constructor
NuTo::IntegrationType2D3NGauss1Ip::IntegrationType2D3NGauss1Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
void NuTo::IntegrationType2D3NGauss1Ip::GetLocalIntegrationPointCoordinates2D(int rIpNum, double rCoordinates[2])const
{
    if (rIpNum==0)
    {
	    rCoordinates[0] = 1./3.;
        rCoordinates[1] = 1./3.;
    }
    else
    {
        throw MechanicsException("[NuTo::IntegrationType2D3NGauss1Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType2D3NGauss1Ip::GetNumIntegrationPoints()const
{
    return 1;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType2D3NGauss1Ip::GetIntegrationPointWeight(int rIpNum)const
{
    return 0.5;
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType2D3NGauss1Ip::GetStrIdentifier()const
{
    return GetStrIdentifierStatic();
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType2D3NGauss1Ip::GetStrIdentifierStatic()
{
    return std::string("2D3NGAUSS1IP");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType2D3NGauss1Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::CellBase::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 3;

    // Point 0
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 1
    VisualizationPointLocalCoordinates.push_back(1.0);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 2
    VisualizationPointLocalCoordinates.push_back(0.);
    VisualizationPointLocalCoordinates.push_back(1.0);

    NumVisualizationCells = 1;

    // cell 0
    VisualizationCellType.push_back(NuTo::CellBase::TRIANGLE);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIP.push_back(0);
}
#endif // ENABLE_VISUALIZE
