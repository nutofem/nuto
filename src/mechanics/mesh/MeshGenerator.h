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
    //! @brief creates a grid from origin to end
    //! @param s structure
    //! @param end coordinates of the end point
    //! @param numDivisions number of grid cells in each direction
    //! @param elementShape element shape (truss, quad, triangle, ...)
    //! @return pair<groupId, interpolationTypeId>
    static std::pair<int, int> Grid(Structure& s, std::vector<double> end, std::vector<int> numDivisions,
                                    Interpolation::eShapeType elementShape);

    //! @brief creates a grid from origin to end with a default element shape
    //!        1D = Truss | 2D = Quad | 3D = Brick
    //! @param s structure
    //! @param end coordinates of the end point
    //! @param numDivisions number of grid cells in each direction
    //! @return pair<groupId, interpolationTypeId>
    static std::pair<int, int> Grid(Structure& s, std::vector<double> end, std::vector<int> numDivisions);

    //! @brief creates a grid from start to end
    //! @param s structure
    //! @param end coordinates of the start point
    //! @param end coordinates of the end point
    //! @param numDivisions number of grid cells in each direction
    //! @param elementShape element shape (truss, quad, triangle, ...)
    //! @return pair<groupId, interpolationTypeId>
    static std::pair<int, int> Grid(Structure& s, std::vector<double> start, std::vector<double> end,
                                    std::vector<int> numDivisions);

    //! @brief creates a grid from start to end with a default element shape
    //!        1D = Truss | 2D = Quad | 3D = Brick
    //! @param s structure
    //! @param end coordinates of the start point
    //! @param end coordinates of the end point
    //! @param numDivisions number of grid cells in each direction
    //! @param elementShape element shape (truss, quad, triangle, ...)
    //! @return pair<groupId, interpolationTypeId>
    static std::pair<int, int> Grid(Structure& s, std::vector<double> start, std::vector<double> end,
                                    std::vector<int> numDivisions, Interpolation::eShapeType elementShape);

    static std::pair<int, int> Grid(Structure& s, Eigen::VectorXd start, Eigen::VectorXd end,
                                    Eigen::VectorXi numDivisions, Interpolation::eShapeType elementShape);

    //! @param radius ... radius of the cylinder
    //! @param height ... height of the cylinder
    static std::function<Eigen::Vector3d(Eigen::Vector3d)> GetCylinderMapping(double radius, double height);
};
}
