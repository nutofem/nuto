#pragma once


#include <functional>
#include <eigen3/Eigen/Core>
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include <vector>

namespace NuTo
{
class Structure;

class MeshGenerator
{
public:
    //! @brief creates a grid from origin to rEnd
    //! @param rS structure
    //! @param rEnd coordinates of the end point
    //! @param rNumDivisions number of grid cells in each direction
    //! @param rElementShape element shape (truss, quad, triangle, ...)
    //! @return pair<groupId, interpolationTypeId>
    static std::pair<int, int> Grid(Structure& rS, std::vector<double> rEnd, std::vector<int> rNumDivisions,
                                    Interpolation::eShapeType rElementShape);

    //! @brief creates a grid from origin to rEnd with a default element shape
    //!        1D = Truss | 2D = Quad | 3D = Brick
    //! @param rS structure
    //! @param rEnd coordinates of the end point
    //! @param rNumDivisions number of grid cells in each direction
    //! @return pair<groupId, interpolationTypeId>
    static std::pair<int, int> Grid(Structure& rS, std::vector<double> rEnd, std::vector<int> rNumDivisions);

    //! @brief creates a grid from rStart to rEnd
    //! @param rS structure
    //! @param rEnd coordinates of the start point
    //! @param rEnd coordinates of the end point
    //! @param rNumDivisions number of grid cells in each direction
    //! @param rElementShape element shape (truss, quad, triangle, ...)
    //! @return pair<groupId, interpolationTypeId>
    static std::pair<int, int> Grid(Structure& rS, std::vector<double> rStart, std::vector<double> rEnd,
                             std::vector<int> rNumDivisions);

    //! @brief creates a grid from rStart to rEnd with a default element shape
    //!        1D = Truss | 2D = Quad | 3D = Brick
    //! @param rS structure
    //! @param rEnd coordinates of the start point
    //! @param rEnd coordinates of the end point
    //! @param rNumDivisions number of grid cells in each direction
    //! @param rElementShape element shape (truss, quad, triangle, ...)
    //! @return pair<groupId, interpolationTypeId>
    static std::pair<int, int> Grid(Structure& rS, std::vector<double> rStart, std::vector<double> rEnd,
                             std::vector<int> rNumDivisions, Interpolation::eShapeType rElementShape);

    //! @param rRadius ... radius of the cylinder
    //! @param rHeight ... height of the cylinder
    static std::function<Eigen::Vector3d(Eigen::Vector3d)> GetCylinderMapping(double rRadius, double rHeight);

};
}
