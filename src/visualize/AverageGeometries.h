#pragma once

#include "visualize/AverageHandler.h"

namespace NuTo
{
namespace Visualize
{
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
    std::vector<Eigen::VectorXd> cornerCoordinates = {Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 1, 0),
                                                      Eigen::Vector3d(1, 0, 0)};
    return {cornerCoordinates, eCellTypes::TRIANGLE};
}

//! @return tetrahedron geometry (0,0,0 -- 1,0,0 -- 0,1,0 -- 0,0,1)
inline AverageGeometry AverageGeometryTetrahedron()
{
    std::vector<Eigen::VectorXd> cornerCoordinates = {Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 0, 0),
                                                      Eigen::Vector3d(0, 1, 0), Eigen::Vector3d(0, 0, 1)};
    return {cornerCoordinates, eCellTypes::TETRAEDER};
}
} // namespace Visualize
} // namespace NuTo
