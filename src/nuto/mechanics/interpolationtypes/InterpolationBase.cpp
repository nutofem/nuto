/*
 * InterpolationBase.cpp
 *
 *  Created on: 17 Mar 2015
 *      Author: ttitsche
 */

#include "nuto/mechanics/interpolationtypes/InterpolationBase.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"

NuTo::InterpolationBase::InterpolationBase(const StructureBase* rStructure, Node::eAttributes rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder) :
    mTypeOrder(rTypeOrder),
    mDofType(rDofType),
    mIsConstitutiveInput(true),
    mIsActive(true),
    mNumDofs(-1),
    mNumNodes(-1),
    mLocalStartIndex(0),
    mUpdateRequired(true),
    mStructure(rStructure)
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

void NuTo::InterpolationBase::Initialize()
{
    mNumNodes = CalculateNumNodes();
    mNumDofs = mNumNodes*GetNumDofsPerNode();


    mNodeCoordinates.resize(mNumNodes);
    for (int i = 0; i < mNumNodes; ++i)
        mNodeCoordinates[i] = CalculateNaturalNodeCoordinates(i);

}
