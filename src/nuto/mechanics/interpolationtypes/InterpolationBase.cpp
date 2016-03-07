/*
 *
 * InterpolationBase.cpp
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
#include "nuto/math/EigenBoostSerialization.h"
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/interpolationtypes/InterpolationBase.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"

NuTo::InterpolationBase::InterpolationBase(Node::eAttributes rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension) :
    mTypeOrder(rTypeOrder),
    mDofType(rDofType),
    mIsConstitutiveInput(true),
    mIsActive(true),
    mNumDofs(-1),
    mNumNodes(-1),
    mLocalStartIndex(0),
    mUpdateRequired(true),
    mDimension(rDimension)
{
}


NuTo::InterpolationBase::InterpolationBase():
    mTypeOrder(NuTo::Interpolation::eTypeOrder::EQUIDISTANT1),
    mDofType(NuTo::Node::eAttributes::COORDINATES),
    mDimension(0)
{
}

const NuTo::Interpolation::eTypeOrder NuTo::InterpolationBase::GetTypeOrder() const
{
    return mTypeOrder;
}

//! @brief calculate and store the shape functions and their derivatives
//! @param rIntegrationType ... integration type
void NuTo::InterpolationBase::UpdateIntegrationType(const IntegrationTypeBase& rIntegrationType)
{
    int numIPs = rIntegrationType.GetNumIntegrationPoints();

    mShapeFunctions.clear();
    mDerivativeShapeFunctionsNatural.clear();

    mShapeFunctions = std::vector<Eigen::VectorXd>(numIPs);
    mDerivativeShapeFunctionsNatural = std::vector<Eigen::MatrixXd>(numIPs);

    int dim = rIntegrationType.GetCoordinateDimension();
    for (int iIP = 0; iIP < numIPs; ++iIP)
    {
        Eigen::VectorXd IPcoordinates;
        switch (dim)
        {
        case 1:
        {
            double coordinate1D;
            rIntegrationType.GetLocalIntegrationPointCoordinates1D(iIP, coordinate1D);
            IPcoordinates.resize(1);
            IPcoordinates(0) = coordinate1D;
        }
            break;
        case 2:
        {
            double coordinate2D[2];
            rIntegrationType.GetLocalIntegrationPointCoordinates2D(iIP, coordinate2D);
            IPcoordinates = Eigen::Matrix<double, 2, 1>(coordinate2D);
        }
            break;
        case 3:
        {
            double coordinate3D[3];
            rIntegrationType.GetLocalIntegrationPointCoordinates3D(iIP, coordinate3D);
            IPcoordinates = Eigen::Matrix<double, 3, 1>(coordinate3D);
        }
            break;
        default:
            throw MechanicsException("[NuTo::InterpolationBase::UpdateIntegrationType] only implemented for dimension 1,2 and 3");
        }
        mShapeFunctions[iIP] = CalculateShapeFunctions(IPcoordinates);
        mDerivativeShapeFunctionsNatural[iIP] = CalculateDerivativeShapeFunctionsNatural(IPcoordinates);
    }
    mUpdateRequired = false;
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

int NuTo::InterpolationBase::GetLocalStartIndex() const
{
    return mLocalStartIndex;
}

int NuTo::InterpolationBase::GetNumNodes() const
{
    assert(mNumNodes != -1);
    return mNumNodes;
}

int NuTo::InterpolationBase::GetNodeIndex(int rNodeDofIndex) const
{
    assert(rNodeDofIndex < mNumNodes);
    assert((unsigned int) rNodeDofIndex < mNodeIndices.size());
    assert(not mUpdateRequired);

    return mNodeIndices[rNodeDofIndex];
}

int NuTo::InterpolationBase::GetNumSurfaceNodes(int rSurface) const
{
    assert((unsigned int) rSurface < mSurfaceNodeIndices.size() && "Surface node indices not build.");
    return mSurfaceNodeIndices[rSurface].size();
}

int NuTo::InterpolationBase::GetSurfaceNodeIndex(int rSurface, int rNodeDofIndex) const
{
    assert((unsigned int) rSurface < mSurfaceNodeIndices.size() && "Surface node indices not build.");
    assert((unsigned int) rNodeDofIndex < mSurfaceNodeIndices[rSurface].size() &&  "Surface node indices not build.");
    assert(not mUpdateRequired);

    return mSurfaceNodeIndices[rSurface][rNodeDofIndex];
}

const Eigen::VectorXd& NuTo::InterpolationBase::GetNaturalNodeCoordinates(int rNodeIndex) const
{
    assert(rNodeIndex < mNumNodes);
    assert((unsigned int) rNodeIndex < mNodeCoordinates.size());
    return mNodeCoordinates[rNodeIndex];
}

const Eigen::VectorXd& NuTo::InterpolationBase::GetShapeFunctions(int rIP) const
{
    assert(rIP < (int )mShapeFunctions.size());
    assert(not mUpdateRequired);
    return mShapeFunctions.at(rIP);
}

const Eigen::MatrixXd& NuTo::InterpolationBase::GetDerivativeShapeFunctionsNatural(int rIP) const
{
    assert(rIP < (int )mDerivativeShapeFunctionsNatural.size());
    assert(not mUpdateRequired);
    return mDerivativeShapeFunctionsNatural.at(rIP);
}

void NuTo::InterpolationBase::CalculateSurfaceNodeIds()
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

bool NuTo::InterpolationBase::NodeIsOnSurface(int rSurface, const Eigen::VectorXd& rNaturalNodeCoordinate) const
{
    auto surfaceEdgeCoordinates = GetSurfaceEdgesCoordinates(rSurface);

    int unsigned dim = rNaturalNodeCoordinate.rows(); // bit of a hack to get the global dimension

    // the following calculations require at least
    // 1 surface egde for 1D  (distance point to point)
    // 2 surface egdes for 2D (distance point to line)
    // 3 surface edges for 3D (distance point to plane)
    assert (surfaceEdgeCoordinates.size() >= dim);

    const auto& pP = rNaturalNodeCoordinate;

    // using the determinant equation
    Eigen::MatrixXd matrix(dim, dim);

    matrix.col(0) = pP - surfaceEdgeCoordinates[0];         // dimension 0 outside the loop, since it involves pP
    for (unsigned int iDim = 1; iDim < dim; ++iDim)                  // note: loop starting at 1
        matrix.col(iDim) = surfaceEdgeCoordinates[iDim] - surfaceEdgeCoordinates[0];

    double det = matrix.determinant(); // determinant is the distance of the point to point/line/plane

    return std::abs(det) < 1.e-10;
}

void NuTo::InterpolationBase::Initialize()
{
    mNumNodes = CalculateNumNodes();
    mNumDofs = mNumNodes*GetNumDofsPerNode();


    mNodeCoordinates.resize(mNumNodes);
    for (int i = 0; i < mNumNodes; ++i)
        mNodeCoordinates[i] = CalculateNaturalNodeCoordinates(i);

}

#ifdef ENABLE_SERIALIZATION
template void NuTo::InterpolationBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::InterpolationBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::InterpolationBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::InterpolationBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::InterpolationBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::InterpolationBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::InterpolationBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize InterpolationBase" << std::endl;
#endif
    ar & boost::serialization::make_nvp("mTypeOrder", const_cast<NuTo::Interpolation::eTypeOrder&>(mTypeOrder));
    ar & boost::serialization::make_nvp("mDofType", const_cast<NuTo::Node::eAttributes&>(mDofType));
    ar & BOOST_SERIALIZATION_NVP(mIsConstitutiveInput);
    ar & BOOST_SERIALIZATION_NVP(mIsActive);
    ar & BOOST_SERIALIZATION_NVP(mNumDofs);
    ar & BOOST_SERIALIZATION_NVP(mNumNodes);

    ar & boost::serialization::make_nvp("mNodeIndices", mNodeIndices);
    ar & BOOST_SERIALIZATION_NVP(mLocalStartIndex);
    ar & boost::serialization::make_nvp("mNodeCoordinates", mNodeCoordinates);
    ar & boost::serialization::make_nvp("mShapeFunctions", mShapeFunctions);
    ar & boost::serialization::make_nvp("mNodeCoordinates", mDerivativeShapeFunctionsNatural);
    ar & BOOST_SERIALIZATION_NVP(mUpdateRequired);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize InterpolationBase" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::InterpolationBase)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::InterpolationBase)
#endif  // ENABLE_SERIALIZATION
