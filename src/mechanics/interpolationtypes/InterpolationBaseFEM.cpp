/*
 *
 * InterpolationBaseFEM.cpp
 *
 *  Created on: 17 Mar 2015
 *      Author: ttitsche
 */

#include <eigen3/Eigen/Dense> // for determinant
#include "mechanics/interpolationtypes/InterpolationBaseFEM.h"
#include "mechanics/nodes/NodeEnum.h"

using namespace NuTo;

InterpolationBaseFEM::InterpolationBaseFEM(Node::eDof rDofType, Interpolation::eTypeOrder rTypeOrder, int rDimension)
    : InterpolationBase::InterpolationBase(rDofType, rTypeOrder, rDimension)
    , mShapeFunctions([=](const Eigen::VectorXd& v) { return CalculateShapeFunctions(v); })
    , mMatrixN([=](const Eigen::VectorXd& v) { return CalculateMatrixN(v); })
    , mDerivativeShapeFunctionsNatural(
              [=](const Eigen::VectorXd& v) { return CalculateDerivativeShapeFunctionsNatural(v); })
{
}

void InterpolationBaseFEM::ClearCache() const
{
    mShapeFunctions.ClearCache();
    mMatrixN.ClearCache();
    mDerivativeShapeFunctionsNatural.ClearCache();
}

Eigen::MatrixXd InterpolationBaseFEM::CalculateMatrixN(const Eigen::VectorXd& rCoordinates) const
{

    auto shapeFunctions = CalculateShapeFunctions(rCoordinates);
    if (GetLocalDimension() == 1)
        return shapeFunctions.transpose();

    int numNodes = GetNumNodes();
    int dimBlock = Node::GetNumComponents(mDofType, mDimension);

    assert(shapeFunctions.rows() == numNodes);

    Eigen::MatrixXd matrixN(dimBlock, numNodes * dimBlock);
    for (int iNode = 0, iBlock = 0; iNode < numNodes; ++iNode, iBlock += dimBlock)
    {
        matrixN.block(0, iBlock, dimBlock, dimBlock) =
                Eigen::MatrixXd::Identity(dimBlock, dimBlock) * shapeFunctions(iNode);
    }

    return matrixN;
}

int InterpolationBaseFEM::GetSplineDegree(int dir) const
{
    throw Exception(__PRETTY_FUNCTION__, "Use 'GetTypeOrder' instead!");
}

const Eigen::VectorXd& InterpolationBaseFEM::GetNaturalNodeCoordinates(int rNodeIndex) const
{
    assert(rNodeIndex < mNumNodes);
    assert((unsigned int)rNodeIndex < mNodeCoordinates.size());
    return mNodeCoordinates[rNodeIndex];
}

const Eigen::VectorXd& InterpolationBaseFEM::ShapeFunctions(const Eigen::VectorXd& naturalCoordinates) const
{
    return mShapeFunctions.Get(naturalCoordinates);
}

const Eigen::MatrixXd& InterpolationBaseFEM::MatrixN(const Eigen::VectorXd& naturalCoordinates) const
{
    return mMatrixN.Get(naturalCoordinates);
}

const Eigen::MatrixXd&
InterpolationBaseFEM::DerivativeShapeFunctionsNatural(const Eigen::VectorXd& naturalCoordinates) const
{
    return mDerivativeShapeFunctionsNatural.Get(naturalCoordinates);
}

void InterpolationBaseFEM::CalculateSurfaceNodeIds()
{
    mSurfaceNodeIndices.clear();
    mSurfaceNodeIndices.resize(GetNumSurfaces());
    for (int iSurface = 0; iSurface < GetNumSurfaces(); ++iSurface)
    {
        auto& surfaceNodeIndices = mSurfaceNodeIndices[iSurface];
        surfaceNodeIndices.reserve(GetNumSurfaceNodes(iSurface));

        for (int iNode = 0; iNode < GetNumNodes(); ++iNode)
        {
            if (NodeIsOnSurface(iSurface, GetNaturalNodeCoordinates(iNode)))
            {
                surfaceNodeIndices.push_back(iNode);
            }
        }
    }
}

bool InterpolationBaseFEM::NodeIsOnSurface(int rSurface, const Eigen::VectorXd& rNaturalNodeCoordinate) const
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

void InterpolationBaseFEM::Initialize()
{
    mNumNodes = CalculateNumNodes();
    mNumDofs = mNumNodes * Node::GetNumComponents(mDofType, mDimension);

    mNodeCoordinates.resize(mNumNodes);
    for (int i = 0; i < mNumNodes; ++i)
        mNodeCoordinates[i] = CalculateNaturalNodeCoordinates(i);
}
