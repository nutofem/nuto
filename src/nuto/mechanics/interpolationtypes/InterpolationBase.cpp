/*
 *
 * InterpolationBase.cpp
 *
 *  Created on: 17 Mar 2015
 *      Author: ttitsche
 */

#include "nuto/mechanics/interpolationtypes/InterpolationBase.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"

NuTo::InterpolationBase::InterpolationBase(Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension) :
    mDofType(rDofType),
    mIsConstitutiveInput(true),
    mIsActive(true),
    mNumDofs(-1),
    mNumNodes(-1),
    mTypeOrder(rTypeOrder),
    mUpdateRequired(true),
    mDimension(rDimension)
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

#ifdef ENABLE_SERIALIZATION
NuTo::InterpolationBase::InterpolationBase():
    mTypeOrder(NuTo::Interpolation::eTypeOrder::EQUIDISTANT1),
    mDofType(NuTo::Node::eDof::COORDINATES),
    mDimension(0)
{
}

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
    ar & boost::serialization::make_nvp("mDofType", const_cast<NuTo::Node::eDof&>(mDofType));
    ar & BOOST_SERIALIZATION_NVP(mIsConstitutiveInput);
    ar & BOOST_SERIALIZATION_NVP(mIsActive);
    ar & BOOST_SERIALIZATION_NVP(mNumDofs);
    ar & BOOST_SERIALIZATION_NVP(mNumNodes);

    ar & boost::serialization::make_nvp("mTypeOrder", const_cast<NuTo::Interpolation::eTypeOrder&>(mTypeOrder));

    ar & BOOST_SERIALIZATION_NVP(mNodeIndices);

    ar & BOOST_SERIALIZATION_NVP(mSurfaceNodeIndices);
    ar & BOOST_SERIALIZATION_NVP(mUpdateRequired);
    ar & boost::serialization::make_nvp("mDimension", const_cast<int&>(mDimension));

#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize InterpolationBase" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::InterpolationBase)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::InterpolationBase)
#endif  // ENABLE_SERIALIZATION
