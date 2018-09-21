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
    //! @brief initizalizes the nodes coordintes with the passed values and id
    //! @param coordinates : inititial coordinates
    //! @param id
    CoordinateNode(Eigen::VectorXd coordinates, int id)
        : mCoordinates({coordinates})
        , mId(id)
    {
    }

    //! @brief initializes a 1D coordinate node with the passed value and id
    //! @param coordinate : initial coordinate value
    //! @param id
    CoordinateNode(double coordinate, int id)
        : mCoordinates(Eigen::VectorXd::Constant(1, coordinate))
        , mId(id)
    {
    }

    //! @brief initizalizes the nodes coordintes with the passed values
    //! and invalid id (-1)
    //! @param coordinates : inititial coordinates
    CoordinateNode(Eigen::VectorXd coordinates)
        : mCoordinates({coordinates})
        , mId(-1)
    {
    }

    //! @brief initializes a 1D coordinate node with the passed value
    //! and invalid id (-1)
    //! @param coordinate : initial coordinate value
    CoordinateNode(double coordinate)
        : mCoordinates(Eigen::VectorXd::Constant(1, coordinate))
        , mId(-1)
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

    int Id() const
    {
        return mId;
    }

private:
    Eigen::VectorXd mCoordinates;
    int mId;
};
} /* NuTo */
