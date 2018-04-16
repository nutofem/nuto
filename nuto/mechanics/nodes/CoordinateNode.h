#pragma once

#include <Eigen/Core>
#include <vector>
#include <cassert>


namespace NuTo
{
//! @brief Store node values and its dof
//! @todo fix sized nodes?
class CoordinateNode
{
public:
    //! @brief initizalizes the nodes coordintes with the passed values
    //! @param coordinates : inititial coordinates
    CoordinateNode(Eigen::VectorXd coordinates)
        : mCoordinates({coordinates})
    {
    }

    //! @brief initializes a 1D coordinate node with the passed value
    //! @param coordinate : initial coordinate value
    CoordinateNode(double coordinate)
        : mCoordinates(Eigen::VectorXd::Constant(1, coordinate))
    {
    }


    const Eigen::VectorXd& GetCoordinates() const
    {
        return mCoordinates;
    }

    void SetCoordinates(Eigen::VectorXd coordinates)
    {
        assert(coordinates.rows() == mCoordinates.rows());
        mCoordinates = coordinates;
    }

    void SetCoordinate(int component, double value)
    {
        assert(component < mCoordinates.rows());
        mCoordinates[component] = value;
    }

    int GetNumValues() const
    {
        return mCoordinates.rows();
    }

private:
    Eigen::VectorXd mCoordinates;
};
} /* NuTo */

