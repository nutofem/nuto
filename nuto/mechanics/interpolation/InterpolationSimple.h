#pragma once
#include <Eigen/Core>
#include "nuto/mechanics/interpolation/TypeDefs.h"
#include <memory>
#include "nuto/math/shapes/Shape.h"
#include "nuto/base/Exception.h"
#include <vector>

namespace NuTo
{


//! @brief Base class for the interpolation. The derived classes provide information about the actual interpolation.
class InterpolationSimple
{
public:
    virtual ~InterpolationSimple() = default; // virtual destructor needed, rule of 5 below...

    virtual std::unique_ptr<InterpolationSimple> Clone() const = 0;

    //! @brief calculates the shape functions
    //! @param naturalIpCoords integration point coordinates in the natural coordinate system
    //! @return vector of shape functions, dimension: [GetNumNodes() x 1]
    virtual Eigen::VectorXd GetShapeFunctions(const NaturalCoords& naturalIpCoords) const = 0;

    //! @brief calculates the derivative shape functions
    //! @param naturalIpCoords integration point coordinates in the natural coordinate system
    //! @return matrix of derivate shape functions, dimension: [GetNumNodes() x local dimension]
    //! @remark 'local dimension' above is the dimension of rLocalIPCoords
    virtual Eigen::MatrixXd GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const = 0;

    //! @brief returns the local node coordinates
    //! @param nodeId local node number
    //! @return node coordinates in the natural coordinate system
    virtual NaturalCoords GetLocalCoords(int nodeId) const = 0;

    //! @brief returns the number of nodes
    //! @return number of nodes
    virtual int GetNumNodes() const = 0;

    virtual const Shape& GetShape() const = 0;

    virtual std::vector<int> EdgeNodeIds(int /* edgeIndex */) const
    {
        throw Exception(__PRETTY_FUNCTION__, "Not implemented");
    }

    virtual std::unique_ptr<InterpolationSimple> EdgeInterpolation(int /* edgeIndex*/) const
    {
        throw Exception(__PRETTY_FUNCTION__, "Not implemented");
    }

    virtual int NumEdges() const
    {
        throw Exception(__PRETTY_FUNCTION__, "Not implemented");
    }
};

//! @brief clone methods that enables a boost::ptr_container<this> to copy itself
//! @param interpolation reference to the InterpolationSimple
//! @return cloned owning raw pointer of interpolation
inline NuTo::InterpolationSimple* new_clone(const NuTo::InterpolationSimple& interpolation)
{
    return interpolation.Clone().release();
}

} /* NuTo */
