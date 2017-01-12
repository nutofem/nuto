#include "mechanics/integrationtypes/IntegrationType3D6NGauss2x3Ip.h"
#include <assert.h>

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
    NumVisualizationPoints = 0;

//    // Point 0
//    VisualizationPointLocalCoordinates.push_back(-1);
//    VisualizationPointLocalCoordinates.push_back(-1);
//    VisualizationPointLocalCoordinates.push_back(-1);
//
//    // Point 1
//    VisualizationPointLocalCoordinates.push_back(+1);
//    VisualizationPointLocalCoordinates.push_back(-1);
//    VisualizationPointLocalCoordinates.push_back(-1);
//
//    // Point 2
//    VisualizationPointLocalCoordinates.push_back(+1);
//    VisualizationPointLocalCoordinates.push_back(+1);
//    VisualizationPointLocalCoordinates.push_back(-1);
//
//    // Point 3
//    VisualizationPointLocalCoordinates.push_back(-1);
//    VisualizationPointLocalCoordinates.push_back(+1);
//    VisualizationPointLocalCoordinates.push_back(-1);
//
//    // Point 4
//    VisualizationPointLocalCoordinates.push_back(-1);
//    VisualizationPointLocalCoordinates.push_back(-1);
//    VisualizationPointLocalCoordinates.push_back(+1);
//
//    // Point 5
//    VisualizationPointLocalCoordinates.push_back(+1);
//    VisualizationPointLocalCoordinates.push_back(-1);
//    VisualizationPointLocalCoordinates.push_back(+1);
//
//    // Point 6
//    VisualizationPointLocalCoordinates.push_back(+1);
//    VisualizationPointLocalCoordinates.push_back(+1);
//    VisualizationPointLocalCoordinates.push_back(+1);
//
//    // Point 7
//    VisualizationPointLocalCoordinates.push_back(-1);
//    VisualizationPointLocalCoordinates.push_back(+1);
//    VisualizationPointLocalCoordinates.push_back(+1);

    NumVisualizationCells = 0;

//    // cell 0
//    VisualizationCellType.push_back(NuTo::eCellTypes::HEXAHEDRON);
//    VisualizationCellsIncidence.push_back(0);
//    VisualizationCellsIncidence.push_back(1);
//    VisualizationCellsIncidence.push_back(2);
//    VisualizationCellsIncidence.push_back(3);
//    VisualizationCellsIncidence.push_back(4);
//    VisualizationCellsIncidence.push_back(5);
//    VisualizationCellsIncidence.push_back(6);
//    VisualizationCellsIncidence.push_back(7);
//    VisualizationCellsIP.push_back(0);
}
#endif // ENABLE_VISUALIZE
