/*
 *
 * InterpolationBase.cpp
 *
 *  Created on: 17 Mar 2015
 *      Author: ttitsche
 */

#include <Eigen/Dense> // for ::determinant()
#include "mechanics/interpolationtypes/InterpolationBase.h"

NuTo::InterpolationBase::InterpolationBase(Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder,
                                           int rDimension)
    : mDofType(rDofType)
    , mIsConstitutiveInput(true)
    , mIsActive(true)
    , mNumDofs(-1)
    , mNumNodes(-1)
    , mTypeOrder(rTypeOrder)
    , mDimension(rDimension)
{
}

NuTo::Interpolation::eTypeOrder NuTo::InterpolationBase::GetTypeOrder() const
{
    return mTypeOrder;
}

bool NuTo::InterpolationBase::IsActive() const
{
    return mIsActive;
}

bool NuTo::InterpolationBase::IsConstitutiveInput() const
{
    return mIsConstitutiveInput;
}

int NuTo::InterpolationBase::GetNumDofs() const
{
    assert(mNumNodes != -1);
    return mNumDofs;
}

int NuTo::InterpolationBase::GetNumNodes() const
{
    assert(mNumNodes != -1);
    return mNumNodes;
}

int NuTo::InterpolationBase::GetNodeIndex(int rNodeDofIndex) const
{
    assert(rNodeDofIndex < mNumNodes);
    assert((unsigned int)rNodeDofIndex < mNodeIndices.size());

    return mNodeIndices[rNodeDofIndex];
}

int NuTo::InterpolationBase::GetNumSurfaceNodes(int rSurface) const
{
    assert((unsigned int)rSurface < mSurfaceNodeIndices.size() && "Surface node indices not build.");
    return mSurfaceNodeIndices[rSurface].size();
}

int NuTo::InterpolationBase::GetSurfaceNodeIndex(int rSurface, int rNodeDofIndex) const
{
    assert((unsigned int)rSurface < mSurfaceNodeIndices.size() && "Surface node indices not build.");
    assert((unsigned int)rNodeDofIndex < mSurfaceNodeIndices[rSurface].size() && "Surface node indices not build.");

    return mSurfaceNodeIndices[rSurface][rNodeDofIndex];
}

bool NuTo::InterpolationBase::NodeIsOnSurface(int rSurface, const Eigen::VectorXd& rNaturalNodeCoordinate) const
{
    auto surfaceEdgeCoordinates = GetSurfaceEdgesCoordinates(rSurface);

    int unsigned dim = rNaturalNodeCoordinate.rows(); // bit of a hack to get the global dimension

    // the following calculations require at least
    // 1 surface egde for 1D  (distance point to point)
    // 2 surface egdes for 2D (distance point to line)
    // 3 surface edges for 3D (distance point to plane)
    assert(surfaceEdgeCoordinates.size() >= dim);

    const auto& pP = rNaturalNodeCoordinate;

    // using the determinant equation
    Eigen::MatrixXd matrix(dim, dim);

    matrix.col(0) = pP - surfaceEdgeCoordinates[0]; // dimension 0 outside the loop, since it involves pP
    for (unsigned int iDim = 1; iDim < dim; ++iDim) // note: loop starting at 1
        matrix.col(iDim) = surfaceEdgeCoordinates[iDim] - surfaceEdgeCoordinates[0];

    double det = matrix.determinant(); // determinant is the distance of the point to point/line/plane

    return std::abs(det) < 1.e-10;
}
