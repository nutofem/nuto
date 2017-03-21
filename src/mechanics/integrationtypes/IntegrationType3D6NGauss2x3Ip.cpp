#include "mechanics/integrationtypes/IntegrationType3D6NGauss2x3Ip.h"
#include <cassert>

#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

//! @brief constructor
NuTo::IntegrationType3D6NGauss2x3Ip::IntegrationType3D6NGauss2x3Ip()
{
}

Eigen::VectorXd NuTo::IntegrationType3D6NGauss2x3Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    /**
     *  x and y [index 0 and 1] taken from 2D3N3IP triangle
     *  z taken from 1D2N2IP truss
     */

    assert(rIpNum>=0 && rIpNum<6);
    switch (rIpNum)
    {
    case 0 : return Eigen::Vector3d({1./6., 1./6., -0.577350269189626});
    case 1 : return Eigen::Vector3d({4./6., 1./6., -0.577350269189626});
    case 2 : return Eigen::Vector3d({1./6., 4./6., -0.577350269189626});
    case 3 : return Eigen::Vector3d({1./6., 1./6.,  0.577350269189626});
    case 4 : return Eigen::Vector3d({4./6., 1./6.,  0.577350269189626});
    case 5 : return Eigen::Vector3d({1./6., 4./6.,  0.577350269189626});
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Ip number out of range.");
    }
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType3D6NGauss2x3Ip::GetNumIntegrationPoints()const
{
    return 6;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType3D6NGauss2x3Ip::GetIntegrationPointWeight(int rIpNum)const
{
    return 1/6.;
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType3D6NGauss2x3Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 21;

    // Point 0
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 1
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 2
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 3
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 4
    VisualizationPointLocalCoordinates.push_back(1./3.);
    VisualizationPointLocalCoordinates.push_back(1./3.);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 5
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 6
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 7
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 8
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 9
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 10
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 11
    VisualizationPointLocalCoordinates.push_back(1./3.);
    VisualizationPointLocalCoordinates.push_back(1./3.);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 12
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 13
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 14
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(1);

    // Point 15
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(1);

    // Point 16
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(1);

    // Point 17
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(1);

    // Point 18
    VisualizationPointLocalCoordinates.push_back(1./3.);
    VisualizationPointLocalCoordinates.push_back(1./3.);
    VisualizationPointLocalCoordinates.push_back(1);

    // Point 19
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(1);

    // Point 20
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(1);


    NumVisualizationCells = 6;

    // cell 1
    VisualizationCellType.push_back(NuTo::eCellTypes::HEXAHEDRON);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIncidence.push_back(7);
    VisualizationCellsIncidence.push_back(8);
    VisualizationCellsIncidence.push_back(11);
    VisualizationCellsIncidence.push_back(10);
    VisualizationCellsIP.push_back(0);


    // cell 2
    VisualizationCellType.push_back(NuTo::eCellTypes::HEXAHEDRON);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(8);
    VisualizationCellsIncidence.push_back(9);
    VisualizationCellsIncidence.push_back(12);
    VisualizationCellsIncidence.push_back(11);
    VisualizationCellsIP.push_back(1);

    // cell 3
    VisualizationCellType.push_back(NuTo::eCellTypes::HEXAHEDRON);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(6);
    VisualizationCellsIncidence.push_back(10);
    VisualizationCellsIncidence.push_back(11);
    VisualizationCellsIncidence.push_back(12);
    VisualizationCellsIncidence.push_back(13);
    VisualizationCellsIP.push_back(2);

    // cell 4
    VisualizationCellType.push_back(NuTo::eCellTypes::HEXAHEDRON);
    VisualizationCellsIncidence.push_back(7);
    VisualizationCellsIncidence.push_back(8);
    VisualizationCellsIncidence.push_back(11);
    VisualizationCellsIncidence.push_back(10);
    VisualizationCellsIncidence.push_back(14);
    VisualizationCellsIncidence.push_back(15);
    VisualizationCellsIncidence.push_back(18);
    VisualizationCellsIncidence.push_back(17);
    VisualizationCellsIP.push_back(3);

    // cell 5
    VisualizationCellType.push_back(NuTo::eCellTypes::HEXAHEDRON);
    VisualizationCellsIncidence.push_back(8);
    VisualizationCellsIncidence.push_back(9);
    VisualizationCellsIncidence.push_back(12);
    VisualizationCellsIncidence.push_back(11);
    VisualizationCellsIncidence.push_back(15);
    VisualizationCellsIncidence.push_back(16);
    VisualizationCellsIncidence.push_back(19);
    VisualizationCellsIncidence.push_back(18);
    VisualizationCellsIP.push_back(4);

    // cell 6
    VisualizationCellType.push_back(NuTo::eCellTypes::HEXAHEDRON);
    VisualizationCellsIncidence.push_back(10);
    VisualizationCellsIncidence.push_back(12);
    VisualizationCellsIncidence.push_back(12);
    VisualizationCellsIncidence.push_back(13);
    VisualizationCellsIncidence.push_back(17);
    VisualizationCellsIncidence.push_back(18);
    VisualizationCellsIncidence.push_back(19);
    VisualizationCellsIncidence.push_back(20);
    VisualizationCellsIP.push_back(5);

}
#endif // ENABLE_VISUALIZE
