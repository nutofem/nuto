/*
 *
 * InterpolationBaseFEM.cpp
 *
 *  Created on: 17 Mar 2015
 *      Author: ttitsche
 */

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "math/CustomBoostSerializationExtensions.h"
#endif // ENABLE_SERIALIZATION

#include <eigen3/Eigen/Dense> // for determinant
#include "mechanics/interpolationtypes/InterpolationBaseFEM.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"
#include "mechanics/nodes/NodeEnum.h"

using namespace NuTo;

InterpolationBaseFEM::InterpolationBaseFEM(Node::eDof rDofType, Interpolation::eTypeOrder rTypeOrder, int rDimension)
    : InterpolationBase::InterpolationBase(rDofType, rTypeOrder, rDimension)
{
}

//! @brief calculate and store the shape functions and their derivatives
//! @param rIntegrationType ... integration type
void InterpolationBaseFEM::UpdateIntegrationType(const IntegrationTypeBase& rIntegrationType)
{
    int numIPs = rIntegrationType.GetNumIntegrationPoints();

    mShapeFunctions.clear();
    mMatrixN.clear();
    mDerivativeShapeFunctionsNatural.clear();

    mShapeFunctions.resize(numIPs);
    mMatrixN.resize(numIPs);
    mDerivativeShapeFunctionsNatural.resize(numIPs);

    for (int iIP = 0; iIP < numIPs; ++iIP)
    {
        Eigen::VectorXd IPcoordinates = rIntegrationType.GetLocalIntegrationPointCoordinates(iIP);
        mShapeFunctions[iIP] = CalculateShapeFunctions(IPcoordinates);
        mMatrixN[iIP] = CalculateMatrixN(IPcoordinates);
        mDerivativeShapeFunctionsNatural[iIP] = CalculateDerivativeShapeFunctionsNatural(IPcoordinates);
    }
    mUpdateRequired = false;
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
    throw MechanicsException(__PRETTY_FUNCTION__, "Use 'GetTypeOrder' instead!");
}

const Eigen::VectorXd& InterpolationBaseFEM::GetNaturalNodeCoordinates(int rNodeIndex) const
{
    assert(rNodeIndex < mNumNodes);
    assert((unsigned int)rNodeIndex < mNodeCoordinates.size());
    return mNodeCoordinates[rNodeIndex];
}

const Eigen::VectorXd& InterpolationBaseFEM::GetShapeFunctions(int rIP) const
{
    assert(rIP < (int)mShapeFunctions.size());
    assert(not mUpdateRequired);
    return mShapeFunctions.at(rIP);
}

const Eigen::MatrixXd& InterpolationBaseFEM::GetMatrixN(int rIP) const
{
    assert(rIP < (int)mMatrixN.size());
    assert(not mUpdateRequired);
    return mMatrixN.at(rIP);
}

const Eigen::MatrixXd& InterpolationBaseFEM::GetDerivativeShapeFunctionsNatural(int rIP) const
{
    assert(rIP < (int)mDerivativeShapeFunctionsNatural.size());
    assert(not mUpdateRequired);
    return mDerivativeShapeFunctionsNatural.at(rIP);
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

#ifdef ENABLE_SERIALIZATION
InterpolationBaseFEM::InterpolationBaseFEM()
{
}

template void InterpolationBaseFEM::serialize(boost::archive::binary_oarchive& ar, const unsigned int version);
template void InterpolationBaseFEM::serialize(boost::archive::xml_oarchive& ar, const unsigned int version);
template void InterpolationBaseFEM::serialize(boost::archive::text_oarchive& ar, const unsigned int version);
template void InterpolationBaseFEM::serialize(boost::archive::binary_iarchive& ar, const unsigned int version);
template void InterpolationBaseFEM::serialize(boost::archive::xml_iarchive& ar, const unsigned int version);
template void InterpolationBaseFEM::serialize(boost::archive::text_iarchive& ar, const unsigned int version);
template <class Archive>
void InterpolationBaseFEM::serialize(Archive& ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize InterpolationBaseFEM" << std::endl;
#endif
    ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(InterpolationBase);
    ar& boost::serialization::make_nvp("mNodeCoordinates", mNodeCoordinates);
    ar& boost::serialization::make_nvp("mShapeFunctions", mShapeFunctions);
    ar& boost::serialization::make_nvp("mMatrixN", mMatrixN);
    ar& boost::serialization::make_nvp("mDerivativeShapeFunctionsNatural", mDerivativeShapeFunctionsNatural);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize InterpolationBaseFEM" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(InterpolationBaseFEM)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(InterpolationBaseFEM)
#endif // ENABLE_SERIALIZATION
