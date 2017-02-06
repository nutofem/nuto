#pragma once


#include <array>
#include <functional>
#include <vector>
#include <eigen3/Eigen/Dense>
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"

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
    template <int TDim>
    static std::pair<int, int> Grid(Structure& rS,
                                    std::array<double, TDim> rEnd,
                                    std::array<int, TDim> rNumDivisions,
                                    Interpolation::eShapeType rElementShape)
    {
        Eigen::Matrix<double, TDim, 1> end;
        Eigen::Matrix<int, TDim, 1> divs;

        for (int i = 0; i < TDim; ++i)
        {
            end[i] = rEnd[i];
            divs[i] = rNumDivisions[i];
        }
        return CreateGrid(rS, end, divs, rElementShape);
    }


    //! @brief creates a grid from origin to rEnd with a default element shape
    //!        1D = Truss | 2D = Quad | 3D = Brick
    //! @param rS structure
    //! @param rEnd coordinates of the end point
    //! @param rNumDivisions number of grid cells in each direction
    //! @return pair<groupId, interpolationTypeId>
    template <int TDim>
    static std::pair<int, int> Grid(Structure& rS,
                                    std::array<double, TDim> rEnd,
                                    std::array<int, TDim> rNumDivisions);

    //! @param rRadius ... radius of the cylinder
    //! @param rHeight ... height of the cylinder
    static std::function<Eigen::Vector3d(Eigen::Vector3d)> GetCylinderMapping(double rRadius, double rHeight);


private:

    //! @brief creates a grid from origin to rEnd
    //! @param rS structure
    //! @param rEnd coordinates of the end point
    //! @param rNumDivisions number of grid cells in each direction
    //! @param rElementShape element shape (truss, quad, triangle, ...)
    //! @return pair<groupId, interpolationTypeId>
    template <int TDim>
    static std::pair<int, int> CreateGrid(Structure& rS,
                                          Eigen::Matrix<double, TDim, 1> rEnd,
                                          Eigen::Matrix<int, TDim, 1> rNumDivisions,
                                          Interpolation::eShapeType rElementShape);
};
}
