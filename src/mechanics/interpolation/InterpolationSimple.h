#pragma once
#include <Eigen/Core>
#include "mechanics/interpolation/TypeDefs.h"
#include <memory>

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
    virtual ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const = 0;

    //! @brief calculates the derivative shape functions
    //! @param rNaturalPCoords integration point coordinates in the natural coordinate system
    //! @return matrix of derivate shape functions, dimension: [GetNumNodes() x local dimension]
    //! @remark 'local dimension' above is the dimension of rLocalIPCoords
    virtual DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const = 0;

    //! @brief returns the local node coordinates
    //! @param nodeId local node number
    //! @return node coordinates in the natural coordinate system
    virtual NaturalCoords GetLocalCoords(int nodeId) const = 0;

    //! @brief returns the number of nodes
    //! @return number of nodes
    virtual int GetNumNodes() const = 0;
};

//! @brief clone methods that enables a boost::ptr_container<this> to copy itself
//! @param interpolation reference to the InterpolationSimple
//! @return cloned owning raw pointer of interpolation
inline NuTo::InterpolationSimple* new_clone(const NuTo::InterpolationSimple& interpolation)
{
    return interpolation.Clone().release();
}

} /* NuTo */
