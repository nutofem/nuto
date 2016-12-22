// $Id: IntegrationType3D4NGauss4Ip.cpp 344 2010-10-19 07:48:41Z arnold2 $

#include "mechanics/integrationtypes/IntegrationType3D4NGauss4Ip.h"
#include <assert.h>

#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

//! @brief constructor
NuTo::IntegrationType3D4NGauss4Ip::IntegrationType3D4NGauss4Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
void NuTo::IntegrationType3D4NGauss4Ip::GetLocalIntegrationPointCoordinates3D(int rIpNum, double rCoordinates[3])const
{
    assert(rIpNum>=0 && rIpNum<4);
    switch (rIpNum)
    {
    case 0:
        rCoordinates[0] = 0.13819660;
        rCoordinates[1] = 0.13819660;
        rCoordinates[2] = 0.13819660;
        break;
    case 1:
        rCoordinates[0] = 0.58541020;
        rCoordinates[1] = 0.13819660;
        rCoordinates[2] = 0.13819660;
        break;
    case 2:
        rCoordinates[0] = 0.13819660;
        rCoordinates[1] = 0.58541020;
        rCoordinates[2] = 0.13819660;
        break;
    case 3:
        rCoordinates[0] = 0.13819660;
        rCoordinates[1] = 0.13819660;
        rCoordinates[2] = 0.58541020;
        break;
    default:
    	throw MechanicsException("[NuTo::IntegrationType3D4NGauss4Ip::GetLocalIntegrationPointCoordinates3D] number of ip out of range.");
    }
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType3D4NGauss4Ip::GetNumIntegrationPoints()const
{
    return 4;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType3D4NGauss4Ip::GetIntegrationPointWeight(int rIpNum)const
{
    return 1./24.;
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType3D4NGauss4Ip::GetStrIdentifier()const
{
    return GetStrIdentifierStatic();
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType3D4NGauss4Ip::GetStrIdentifierStatic()
{
    return std::string("3D4NGAUSS4IP");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType3D4NGauss4Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
	NumVisualizationPoints = 15; //Identical with the nodes + centroid of faces + centroid

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

    // Point 4
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 5
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 6
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 7
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0.5);

    // Point 8
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(0.5);

    // Point 9
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(0.0);
    VisualizationPointLocalCoordinates.push_back(0.5);

    // Point 10 (centroid of  0 1 3)
    VisualizationPointLocalCoordinates.push_back(0.3333333);
    VisualizationPointLocalCoordinates.push_back(0.0);
    VisualizationPointLocalCoordinates.push_back(0.3333333);

    // Point 11 (centroid of 1 2 3)
    VisualizationPointLocalCoordinates.push_back(0.3333333);
    VisualizationPointLocalCoordinates.push_back(0.3333333);
    VisualizationPointLocalCoordinates.push_back(0.3333333);

    // Point 12 (centoid 0 2 3)
    VisualizationPointLocalCoordinates.push_back(0.);
    VisualizationPointLocalCoordinates.push_back(0.3333333);
    VisualizationPointLocalCoordinates.push_back(0.3333333);

    // Point 13 (centroid 0 1 2)
    VisualizationPointLocalCoordinates.push_back(0.3333333);
    VisualizationPointLocalCoordinates.push_back(0.3333333);
    VisualizationPointLocalCoordinates.push_back(0.0);

    // Point 14 (centroid)
    VisualizationPointLocalCoordinates.push_back(0.25);
    VisualizationPointLocalCoordinates.push_back(0.25);
    VisualizationPointLocalCoordinates.push_back(0.25);

    NumVisualizationCells = 4;

    // cell 0
    VisualizationCellType.push_back(NuTo::eCellTypes::HEXAHEDRON);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(13);
    VisualizationCellsIncidence.push_back(6);
    VisualizationCellsIncidence.push_back(7);
    VisualizationCellsIncidence.push_back(10);
    VisualizationCellsIncidence.push_back(14);
    VisualizationCellsIncidence.push_back(12);
    VisualizationCellsIP.push_back(0);

    // cell 1
    VisualizationCellType.push_back(NuTo::eCellTypes::HEXAHEDRON);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(13);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(9);
    VisualizationCellsIncidence.push_back(11);
    VisualizationCellsIncidence.push_back(14);
    VisualizationCellsIncidence.push_back(10);
    VisualizationCellsIP.push_back(1);

    // cell 2
    VisualizationCellType.push_back(NuTo::eCellTypes::HEXAHEDRON);
    VisualizationCellsIncidence.push_back(6);
    VisualizationCellsIncidence.push_back(13);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(12);
    VisualizationCellsIncidence.push_back(14);
    VisualizationCellsIncidence.push_back(11);
    VisualizationCellsIncidence.push_back(8);
    VisualizationCellsIP.push_back(2);

    // cell 3
    VisualizationCellType.push_back(NuTo::eCellTypes::HEXAHEDRON);
    VisualizationCellsIncidence.push_back(7);
    VisualizationCellsIncidence.push_back(12);
    VisualizationCellsIncidence.push_back(8);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIncidence.push_back(10);
    VisualizationCellsIncidence.push_back(14);
    VisualizationCellsIncidence.push_back(11);
    VisualizationCellsIncidence.push_back(9);
    VisualizationCellsIP.push_back(3);
}
#endif // ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType3D4NGauss4Ip)
#endif
