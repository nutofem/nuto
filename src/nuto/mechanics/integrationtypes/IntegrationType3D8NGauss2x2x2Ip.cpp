// $Id$

#include "nuto/mechanics/integrationtypes/IntegrationType3D8NGauss2x2x2Ip.h"
#include <assert.h>


//! @brief constructor
NuTo::IntegrationType3D8NGauss2x2x2Ip::IntegrationType3D8NGauss2x2x2Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
void NuTo::IntegrationType3D8NGauss2x2x2Ip::GetLocalIntegrationPointCoordinates3D(int rIpNum, double rCoordinates[3])const
{
    assert(rIpNum>=0 && rIpNum<=7);
    switch (rIpNum)
    {
    case 0 :
        rCoordinates[0] = -0.577350269189626;
        rCoordinates[1] = -0.577350269189626;
        rCoordinates[2] = -0.577350269189626;
        break;
    case 1 :
        rCoordinates[0] = +0.577350269189626;
        rCoordinates[1] = -0.577350269189626;
        rCoordinates[2] = -0.577350269189626;
        break;
    case 2 :
        rCoordinates[0] = -0.577350269189626;
        rCoordinates[1] = +0.577350269189626;
        rCoordinates[2] = -0.577350269189626;
        break;
    case 3 :
        rCoordinates[0] = +0.577350269189626;
        rCoordinates[1] = +0.577350269189626;
        rCoordinates[2] = -0.577350269189626;
        break;
    case 4 :
        rCoordinates[0] = -0.577350269189626;
        rCoordinates[1] = -0.577350269189626;
        rCoordinates[2] = +0.577350269189626;
        break;
    case 5 :
        rCoordinates[0] = +0.577350269189626;
        rCoordinates[1] = -0.577350269189626;
        rCoordinates[2] = +0.577350269189626;
        break;
    case 6 :
        rCoordinates[0] = -0.577350269189626;
        rCoordinates[1] = +0.577350269189626;
        rCoordinates[2] = +0.577350269189626;
        break;
    case 7 :
        rCoordinates[0] = +0.577350269189626;
        rCoordinates[1] = +0.577350269189626;
        rCoordinates[2] = +0.577350269189626;
        break;
    default:
        throw MechanicsException("[NuTo::IntegrationType1D2NGauss2Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType3D8NGauss2x2x2Ip::GetNumIntegrationPoints()const
{
    return 8;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType3D8NGauss2x2x2Ip::GetIntegrationPointWeight(int rIpNum)const
{
    return 1;
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType3D8NGauss2x2x2Ip::GetStrIdentifier()const
{
    return GetStrIdentifierStatic();
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType3D8NGauss2x2x2Ip::GetStrIdentifierStatic()
{
    return std::string("3D8NGAUSS2x2x2IP");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType3D8NGauss2x2x2Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::CellBase::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 27;

    // Point 0
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 1
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 2
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 3
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 4
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 5
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 6
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 7
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 8
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 9
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 10
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 11
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 12
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 13
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 14
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 15
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 16
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 17
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 18
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(1);

    // Point 19
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(1);

    // Point 20
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(1);

    // Point 21
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(1);

    // Point 22
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(1);

    // Point 23
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(1);

    // Point 24
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(1);

    // Point 25
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(1);

    // Point 26
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(1);


    NumVisualizationCells = 8;

    // cell 1
    VisualizationCellType.push_back(NuTo::CellBase::HEXAHEDRON);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIncidence.push_back(9);
    VisualizationCellsIncidence.push_back(10);
    VisualizationCellsIncidence.push_back(13);
    VisualizationCellsIncidence.push_back(12);
    VisualizationCellsIP.push_back(0);


    // cell 2
    VisualizationCellType.push_back(NuTo::CellBase::HEXAHEDRON);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(10);
    VisualizationCellsIncidence.push_back(11);
    VisualizationCellsIncidence.push_back(14);
    VisualizationCellsIncidence.push_back(13);
    VisualizationCellsIP.push_back(1);

    // cell 3
    VisualizationCellType.push_back(NuTo::CellBase::HEXAHEDRON);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(7);
    VisualizationCellsIncidence.push_back(6);
    VisualizationCellsIncidence.push_back(12);
    VisualizationCellsIncidence.push_back(13);
    VisualizationCellsIncidence.push_back(16);
    VisualizationCellsIncidence.push_back(15);
    VisualizationCellsIP.push_back(2);

    // cell 4
    VisualizationCellType.push_back(NuTo::CellBase::HEXAHEDRON);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(8);
    VisualizationCellsIncidence.push_back(7);
    VisualizationCellsIncidence.push_back(13);
    VisualizationCellsIncidence.push_back(14);
    VisualizationCellsIncidence.push_back(17);
    VisualizationCellsIncidence.push_back(16);
    VisualizationCellsIP.push_back(3);

    // cell 5
    VisualizationCellType.push_back(NuTo::CellBase::HEXAHEDRON);
    VisualizationCellsIncidence.push_back(9);
    VisualizationCellsIncidence.push_back(10);
    VisualizationCellsIncidence.push_back(13);
    VisualizationCellsIncidence.push_back(12);
    VisualizationCellsIncidence.push_back(18);
    VisualizationCellsIncidence.push_back(19);
    VisualizationCellsIncidence.push_back(22);
    VisualizationCellsIncidence.push_back(21);
    VisualizationCellsIP.push_back(4);

    // cell 6
    VisualizationCellType.push_back(NuTo::CellBase::HEXAHEDRON);
    VisualizationCellsIncidence.push_back(10);
    VisualizationCellsIncidence.push_back(11);
    VisualizationCellsIncidence.push_back(14);
    VisualizationCellsIncidence.push_back(13);
    VisualizationCellsIncidence.push_back(19);
    VisualizationCellsIncidence.push_back(20);
    VisualizationCellsIncidence.push_back(23);
    VisualizationCellsIncidence.push_back(22);
    VisualizationCellsIP.push_back(5);

    // cell 7
    VisualizationCellType.push_back(NuTo::CellBase::HEXAHEDRON);
    VisualizationCellsIncidence.push_back(12);
    VisualizationCellsIncidence.push_back(13);
    VisualizationCellsIncidence.push_back(16);
    VisualizationCellsIncidence.push_back(15);
    VisualizationCellsIncidence.push_back(21);
    VisualizationCellsIncidence.push_back(22);
    VisualizationCellsIncidence.push_back(25);
    VisualizationCellsIncidence.push_back(24);
    VisualizationCellsIP.push_back(6);

    // cell 8
    VisualizationCellType.push_back(NuTo::CellBase::HEXAHEDRON);
    VisualizationCellsIncidence.push_back(13);
    VisualizationCellsIncidence.push_back(14);
    VisualizationCellsIncidence.push_back(17);
    VisualizationCellsIncidence.push_back(16);
    VisualizationCellsIncidence.push_back(22);
    VisualizationCellsIncidence.push_back(23);
    VisualizationCellsIncidence.push_back(26);
    VisualizationCellsIncidence.push_back(25);
    VisualizationCellsIP.push_back(7);
}
#endif // ENABLE_VISUALIZE
