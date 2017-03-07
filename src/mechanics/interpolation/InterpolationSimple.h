#pragma once
#include <eigen3/Eigen/Core>
#include "mechanics/interpolation/TypeDefs.h"

namespace NuTo
{


//! @brief Base class for the interpolation. The derived classes provide information about the actual interpolation.
//!        [TODO] Implement 'caching' of the results. Requests for the same local node coordinates
//!        must be calculated only once.
class InterpolationSimple
{
public:
    InterpolationSimple()                     = default;
    virtual ~InterpolationSimple()            = default; // virtual destructor needed, rule of 5 below...
    InterpolationSimple(const InterpolationSimple&) = default;
    InterpolationSimple(InterpolationSimple&&)      = default;
    InterpolationSimple& operator=(const InterpolationSimple&) = default;
    InterpolationSimple& operator=(InterpolationSimple&&) = default;


    NMatrix GetN(const NaturalCoords& rNaturalIPCoords) const
    {
        int dim = GetDofDimension();
        Eigen::MatrixXd N(dim, dim * GetNumNodes());

        auto shapeFunctions = GetShapeFunctions(rNaturalIPCoords);

        for (size_t i = 0; i < GetNumNodes(); ++i)
            N.block(0, i * dim, dim, dim) = Eigen::MatrixXd::Identity(dim, dim) * shapeFunctions[i];
        return N;
    }

    //! @brief calculates the shape functions
    //! @param rNaturalIPCoords integration point coordinates in the natural coordinate system
    //! @return vector of shape functions, dimension: [GetNumNodes() x 1]
    virtual ShapeFunctions GetShapeFunctions(const NaturalCoords& rNaturalIPCoords) const = 0;

    //! @brief calculates the derivative shape functions
    //! @param rNaturalPCoords integration point coordinates in the natural coordinate system
    //! @return matrix of derivate shape functions, dimension: [GetNumNodes() x local dimension]
    //! @remark 'local dimension' above is the dimension of rLocalIPCoords
    virtual DerivativeShapeFunctionsNatural
    GetDerivativeShapeFunctions(const NaturalCoords& rNaturalIPCoords) const = 0;

    //! @brief returns the local node coordinates
    //! @param rNodeId local node number
    //! @return node coordinates in the natural coordinate system
    virtual NaturalCoords GetLocalCoords(int rNodeId) const = 0;

    //! @brief returns the number of nodes
    //! @return number of nodes
    virtual int GetNumNodes() const = 0;

    virtual int GetDofDimension() const = 0;
};
} /* NuTo */
