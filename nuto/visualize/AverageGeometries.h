#pragma once

#include <vector>
#include "nuto/visualize/VisualizeEnum.h"
#include "Eigen/Core"

namespace NuTo
{
namespace Visualize
{

//! Geometry description of one cell for the AverageHandler
struct AverageGeometry
{
    std::vector<Eigen::VectorXd> cornerCoordinates;
    eCellTypes cellType;
};

//! @return quad geometry from (-1,-1) to (1,1)
inline AverageGeometry AverageGeometryLine()
{
    std::vector<Eigen::VectorXd> cornerCoordinates = {Eigen::Vector3d(-1, 0, 0), Eigen::Vector3d(1, 0, 0)};
    return {cornerCoordinates, eCellTypes::LINE};
}

//! @return quad geometry from (-1,-1) to (1,1)
inline AverageGeometry AverageGeometryQuad()
{
    std::vector<Eigen::VectorXd> cornerCoordinates = {Eigen::Vector3d(-1, -1, 0), Eigen::Vector3d(1, -1, 0),
                                                      Eigen::Vector3d(1, 1, 0), Eigen::Vector3d(-1, 1, 0)};
    return {cornerCoordinates, eCellTypes::QUAD};
}

//! @return triangle geometry (0,0 -- 1,0 -- 0,1)
inline AverageGeometry AverageGeometryTriangle()
{
    std::vector<Eigen::VectorXd> cornerCoordinates = {Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 0, 0),
                                                      Eigen::Vector3d(0, 1, 0)};
    return {cornerCoordinates, eCellTypes::TRIANGLE};
}

//! @return tetrahedron geometry (0,0,0 -- 1,0,0 -- 0,1,0 -- 0,0,1)
inline AverageGeometry AverageGeometryTetrahedron()
{
    std::vector<Eigen::VectorXd> cornerCoordinates = {Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 0, 0),
                                                      Eigen::Vector3d(0, 1, 0), Eigen::Vector3d(0, 0, 1)};
    return {cornerCoordinates, eCellTypes::TETRAEDER};
}

//! @return prism geometry (product of Triangle() in x,y, Line() in z)
inline AverageGeometry AverageGeometryPrism()
{
    std::vector<Eigen::VectorXd> cornerCoordinates = {Eigen::Vector3d(0, 0, -1), Eigen::Vector3d(1, 0, -1),
                                                      Eigen::Vector3d(0, 1, -1), Eigen::Vector3d(0, 0, 1),
                                                      Eigen::Vector3d(1, 0, 1),  Eigen::Vector3d(0, 1, 1)};
    return {cornerCoordinates, eCellTypes::WEDGE};
}
} // namespace Visualize
} // namespace NuTo
