/*
 * InterpolationType.cpp
 *
 *  Created on: 31 Mar 2015
 *      Author: ttitsche
 */

#include <iostream>

#include "mechanics/MechanicsException.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"

#include "mechanics/interpolationtypes/Interpolation2DTriangle.h"
#include "mechanics/interpolationtypes/Interpolation2DQuad.h"
#include "mechanics/interpolationtypes/Interpolation3DTetrahedron.h"
#include "mechanics/interpolationtypes/Interpolation3DBrick.h"
#include "mechanics/interpolationtypes/Interpolation3DPrism.h"
#include "mechanics/interpolationtypes/Interpolation1DTruss.h"
#include "mechanics/interpolationtypes/Interpolation1DInterface.h"
#include "mechanics/interpolationtypes/Interpolation1DIGA.h"
#include "mechanics/interpolationtypes/Interpolation2DIGA.h"

#include "mechanics/nodes/NodeEnum.h"

#include <iomanip>

NuTo::InterpolationType::InterpolationType(NuTo::Interpolation::eShapeType rShapeType, int rDimension)
    : mShapeType(rShapeType)
    , mNumDofs(0)
    , mNumActiveDofs(0)
    , mDimension(rDimension)
{
}

NuTo::InterpolationType::~InterpolationType()
{
}

const NuTo::InterpolationBase& NuTo::InterpolationType::Get(const Node::eDof& rDofType) const
{
    try
    {
        return mInterpolations.at(rDofType);
    }
    catch (boost::bad_ptr_container_operation& e)
    {
        std::cout << e.what() << std::endl;
        throw NuTo::MechanicsException("[NuTo::InterpolationType::Get] Dof " + Node::DofToString(rDofType) +
                                       " is not a member of this interpolation type. Add it first.");
    }
}

NuTo::InterpolationBase& NuTo::InterpolationType::GetNonConst(Node::eDof rDofType)
{
    auto interpolationTypeIterator = mInterpolations.find(rDofType);

    if (interpolationTypeIterator == mInterpolations.end())
        throw NuTo::MechanicsException("[NuTo::InterpolationType::Get] Dof " + Node::DofToString(rDofType) +
                                       " is not a member of this interpolation type. Add it first.");

    return *(interpolationTypeIterator->second);
}

void NuTo::InterpolationType::AddDofInterpolation(Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder,
                                                  const Eigen::VectorXi& rDegree,
                                                  const std::vector<Eigen::VectorXd>& rKnots,
                                                  const Eigen::MatrixXd& rWeights)
{
    if (IsDof(rDofType))
        throw NuTo::MechanicsException("[NuTo::InterpolationTypeBase::AddDofInterpolation] Dof " +
                                       NuTo::Node::DofToString(rDofType) + " exists.");

    InterpolationBase* newType;
    switch (mShapeType)
    {
    case Interpolation::eShapeType::SPRING:
    case Interpolation::eShapeType::TRUSS1D:
    case Interpolation::eShapeType::TRUSSXD:
    case Interpolation::eShapeType::TRIANGLE2D:
    case Interpolation::eShapeType::QUAD2D:
    case Interpolation::eShapeType::TETRAHEDRON3D:
    case Interpolation::eShapeType::BRICK3D:
    case Interpolation::eShapeType::INTERFACE:
        throw NuTo::MechanicsException("[NuTo::InterpolationTypeBase::AddDofInterpolation] This method is for IG "
                                       "interpolation, please use the 'AddDofInterpolation(Node::eDof rDofType, "
                                       "NuTo::Interpolation::eTypeOrder rTypeOrder)'.");
        break;
    case Interpolation::eShapeType::IGA1D:
        newType = new Interpolation1DIGA(rDofType, rTypeOrder, mDimension, rDegree(0), rKnots[0], rWeights);
        break;
    case Interpolation::eShapeType::IGA2D:
        newType = new Interpolation2DIGA(rDofType, rTypeOrder, mDimension, rDegree, rKnots[0], rKnots[1], rWeights);
        break;
    default:
        throw NuTo::MechanicsException("[NuTo::InterpolationType::AddDofInterpolation] ShapeType " +
                                       NuTo::Interpolation::ShapeTypeToString(mShapeType) + " not implemented.");
    }

    mInterpolations.insert(rDofType, newType);
    mDofs.insert(rDofType);

    mNumDofs += newType->GetNumDofs();

    if (rDofType == Node::eDof::COORDINATES)
    {
        SetIsActive(false, rDofType);
        SetIsConstitutiveInput(false, rDofType);
    }
    else
    {
        SetIsActive(true, rDofType);
        SetIsConstitutiveInput(true, rDofType);
    }

    // IGA: same number of dofs per element => just insert the new dof into the set
    // calculate mNodeDofs, mNodeCoordinates and mNodeIndices

    size_t size = newType->GetNumNodes();

    if (mNodeDofs.size() == 0)
        mNodeDofs.resize(size);
    else if (mNodeDofs.size() != size)
        throw NuTo::MechanicsException(
                "[NuTo::InterpolationType::AddDofInterpolation] The number of dofs per IGA element isn't equal!");

    mNodeCoordinates.resize(size);

    newType->mNodeIndices.resize(size);

    for (size_t i = 0; i < size; i++)
    {
        newType->mNodeIndices[i] = i;
        mNodeDofs[i].insert(rDofType);
    }
}

void NuTo::InterpolationType::AddDofInterpolation(Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder)
{
    if (IsDof(rDofType))
        throw NuTo::MechanicsException("[NuTo::InterpolationTypeBase::AddDofInterpolation] Dof " +
                                       NuTo::Node::DofToString(rDofType) + " exists.");

    InterpolationBase* newType;
    switch (mShapeType)
    {
    case Interpolation::eShapeType::SPRING:
    case Interpolation::eShapeType::TRUSS1D:
    case Interpolation::eShapeType::TRUSSXD:
        newType = new Interpolation1DTruss(rDofType, rTypeOrder, mDimension);
        break;
    case Interpolation::eShapeType::TRIANGLE2D:
        newType = new Interpolation2DTriangle(rDofType, rTypeOrder, mDimension);
        break;
    case Interpolation::eShapeType::QUAD2D:
        newType = new Interpolation2DQuad(rDofType, rTypeOrder, mDimension);
        break;
    case Interpolation::eShapeType::TETRAHEDRON3D:
        newType = new Interpolation3DTetrahedron(rDofType, rTypeOrder, mDimension);
        break;
    case Interpolation::eShapeType::BRICK3D:
        newType = new Interpolation3DBrick(rDofType, rTypeOrder, mDimension);
        break;
    case Interpolation::eShapeType::PRISM3D:
        newType = new Interpolation3DPrism(rDofType, rTypeOrder, mDimension);
        break;
    case Interpolation::eShapeType::INTERFACE:
        newType = new Interpolation1DInterface(rDofType, rTypeOrder, mDimension);
        break;
    default:
        throw NuTo::MechanicsException("[NuTo::InterpolationType::AddDofInterpolation] ShapeType " +
                                       NuTo::Interpolation::ShapeTypeToString(mShapeType) + " not implemented.");
    }

    mInterpolations.insert(rDofType, newType);
    mDofs.insert(rDofType);

    mNumDofs += newType->GetNumDofs();

    if (rDofType == Node::eDof::COORDINATES)
    {
        SetIsActive(false, rDofType);
        SetIsConstitutiveInput(false, rDofType);
    }
    else
    {
        SetIsActive(true, rDofType);
        SetIsConstitutiveInput(true, rDofType);
    }

    // calculate mNodeDofs, mNodeCoordinates and mNodeIndices
    mNodeDofs.resize(mNodeDofs.size() + newType->GetNumNodes());

    newType->mNodeIndices.resize(newType->GetNumNodes());

    for (int iNode = 0; iNode < newType->GetNumNodes(); ++iNode)
    {
        Eigen::VectorXd newCoordinates = newType->GetNaturalNodeCoordinates(iNode);

        // check if unique:
        bool isUnique = true;
        for (unsigned int iExistingNode = 0; iExistingNode < mNodeCoordinates.size(); ++iExistingNode)
        {
            const Eigen::VectorXd& existingCoordinates = mNodeCoordinates[iExistingNode];

            if (CoordinatesAreEqual(existingCoordinates, newCoordinates))
            {
                isUnique = false;
                // link the node to the global node indexing
                newType->mNodeIndices[iNode] = iExistingNode;
                // add the corresponding dof
                mNodeDofs[iExistingNode].insert(rDofType);
            }
        }

        if (isUnique)
        {
            mNodeCoordinates.push_back(newCoordinates);
            int newIndex = mNodeCoordinates.size() - 1;
            // link the node to the global node indexing
            newType->mNodeIndices[iNode] = newIndex;
            // add the corresponding dof
            mNodeDofs[newIndex].insert(rDofType);
        }
    }

    mNumActiveDofs = 0;
    for (auto dof : mActiveDofs)
        mNumActiveDofs += Get(dof).GetNumDofs();

    if (rDofType == Node::eDof::COORDINATES)
        UpdateNodeRenumberingIndices();

    if (mShapeType != Interpolation::eShapeType::INTERFACE)
    {
        // update the surface node ids for each dof
        newType->CalculateSurfaceNodeIds();

        // update the surface node ids for all nodes
        mSurfaceNodeIndices.clear();
        mSurfaceNodeIndices.resize(GetNumSurfaces());
        for (int iSurface = 0; iSurface < GetNumSurfaces(); ++iSurface)
        {
            auto& surfaceNodeIndices = mSurfaceNodeIndices[iSurface];

            for (int iNode = 0; iNode < GetNumNodes(); ++iNode)
                if (newType->NodeIsOnSurface(iSurface, GetNaturalNodeCoordinates(iNode)))
                    surfaceNodeIndices.push_back(iNode);
        }
    }
}


const NuTo::Interpolation::eShapeType NuTo::InterpolationType::GetShapeType() const
{
    return mShapeType;
}

void NuTo::InterpolationType::ClearCache() const
{
    for (const auto& interpolation : mInterpolations)
    {
        interpolation.second->ClearCache();
    }
}

NuTo::eIntegrationType NuTo::InterpolationType::GetStandardIntegrationType() const
{
    return Get(GetDofWithHighestStandardIntegrationOrder()).GetStandardIntegrationType();
}

NuTo::Node::eDof NuTo::InterpolationType::GetDofWithHighestStandardIntegrationOrder() const
{
    int maxOrder = 0;
    Node::eDof maxDof = *GetDofs().begin();
    for (Node::eDof dof : GetDofs())
    {
        Interpolation::eTypeOrder order = Get(dof).GetTypeOrder();
        int currentOrder = 0;
        switch (order)
        {
        case Interpolation::eTypeOrder::EQUIDISTANT1:
            currentOrder = 1;
            break;
        case Interpolation::eTypeOrder::LOBATTO2:
        case Interpolation::eTypeOrder::EQUIDISTANT2:
            currentOrder = 2;
            break;
        case Interpolation::eTypeOrder::LOBATTO3:
        case Interpolation::eTypeOrder::EQUIDISTANT3:
            currentOrder = 3;
            break;
        case Interpolation::eTypeOrder::LOBATTO4:
        case Interpolation::eTypeOrder::EQUIDISTANT4:
            currentOrder = 4;
            break;
        case Interpolation::eTypeOrder::SPLINE:
        {
            int dim = Get(dof).GetLocalDimension();
            if (dim == 1)
                currentOrder = Get(dof).GetSplineDegree(0) + 1;
            else if (dim == 2)
                currentOrder = Get(dof).GetSplineDegree(0) + Get(dof).GetSplineDegree(1) + 2;
            else
                throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "3D IGA is not implemented yet.");
        }
        break;
        default:
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Standard integration type for " +
                                                                        Interpolation::TypeOrderToString(order) +
                                                                        " is not implemented.");
        }

        if (currentOrder > maxOrder)
        {
            maxOrder = currentOrder;
            maxDof = dof;
        }
    }
    return maxDof;
}

Eigen::VectorXi NuTo::InterpolationType::GetSurfaceNodeIndices(int rSurface) const
{
    assert(IsDof(Node::eDof::COORDINATES));
    const InterpolationBase& interpolation = Get(Node::eDof::COORDINATES);

    Interpolation::eTypeOrder order = interpolation.GetTypeOrder();

    if (order != Interpolation::eTypeOrder::SPLINE)
    {
        const auto& surfaceEdgesCoordinates = interpolation.GetSurfaceEdgesCoordinates(rSurface);

        Eigen::VectorXi surfaceNodeIndices(surfaceEdgesCoordinates.size());

        for (unsigned int i = 0; i < surfaceEdgesCoordinates.size(); ++i)
        {
            const Eigen::VectorXd& naturalSurfaceEdgeCoordinate = surfaceEdgesCoordinates[i];
            Eigen::VectorXd shapeFunctions = interpolation.CalculateShapeFunctions(naturalSurfaceEdgeCoordinate);
            assert(std::abs(shapeFunctions.norm() - 1) < 1.e-8);
            int indexNodeCoordinate; // index in the coordinate interpolation

#ifndef NDEBUG
            double value =
                    shapeFunctions.maxCoeff(&indexNodeCoordinate); // find the index where the shape function is 1
            assert(std::abs(value - 1) < 1.e-8);
#else
            shapeFunctions.maxCoeff(&indexNodeCoordinate); // find the index where the shape function is 1
#endif
            surfaceNodeIndices(i) =
                    interpolation.GetNodeIndex(indexNodeCoordinate); // index in the Element::mNodes vector
        }
        return surfaceNodeIndices;
    }
    else
    {
        return interpolation.GetSurfaceNodeIndices(rSurface);
    }
}

int NuTo::InterpolationType::GetNumSurfaces() const
{
    assert(mInterpolations.size() != 0);
    return mInterpolations.begin()->second->GetNumSurfaces();
}

void NuTo::InterpolationType::SetIsActive(bool rIsActiveDof, Node::eDof rDofType)
{
    if (not IsDof(rDofType))
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__,
                                       "Dof " + Node::DofToString(rDofType) +
                                               " is not a member of this interpolation type. Add it first.");

    if (rIsActiveDof)
        mActiveDofs.insert(rDofType);
    else
        mActiveDofs.erase(rDofType);

    GetNonConst(rDofType).mIsActive = rIsActiveDof;

    mNumActiveDofs = 0;
    for (auto dof : mActiveDofs)
        mNumActiveDofs += Get(dof).GetNumDofs();
}

bool NuTo::InterpolationType::IsActive(const Node::eDof& rDofType) const
{
    return Get(rDofType).IsActive();
}

void NuTo::InterpolationType::SetIsConstitutiveInput(bool rIsConstitutiveInput, Node::eDof rDofType)
{
    GetNonConst(rDofType).mIsConstitutiveInput = rIsConstitutiveInput;
}

bool NuTo::InterpolationType::IsConstitutiveInput(const Node::eDof& rDofType) const
{
    if (not IsDof(rDofType))
        return false;
    return Get(rDofType).IsConstitutiveInput();
}

bool NuTo::InterpolationType::IsDof(const Node::eDof& rDofType) const
{
    return mDofs.find(rDofType) != mDofs.end();
}

int NuTo::InterpolationType::GetNumDofs() const
{
    return mNumDofs;
}

//! @brief returns the number of active
int NuTo::InterpolationType::GetNumActiveDofs() const
{
    return mNumActiveDofs;
}

const std::set<NuTo::Node::eDof>& NuTo::InterpolationType::GetActiveDofs() const
{
    return mActiveDofs;
}

const std::set<NuTo::Node::eDof>& NuTo::InterpolationType::GetDofs() const
{
    return mDofs;
}

std::set<NuTo::Node::eDof> NuTo::InterpolationType::GetNodeDofs(int rNodeIndex) const
{
    assert((unsigned int)rNodeIndex < mNodeDofs.size());
    return mNodeDofs[rNodeIndex];
}

const Eigen::VectorXd& NuTo::InterpolationType::GetNaturalNodeCoordinates(int rNodeIndex) const
{
    assert((unsigned int)rNodeIndex < mNodeCoordinates.size());
    return mNodeCoordinates[rNodeIndex];
}

int NuTo::InterpolationType::GetNumNodes() const
{
    return mNodeCoordinates.size();
}

int NuTo::InterpolationType::GetNumSurfaceNodes(int rSurface) const
{
    assert((unsigned int)rSurface < mSurfaceNodeIndices.size() && "Surface node indices not build.");
    return mSurfaceNodeIndices[rSurface].size();
}

int NuTo::InterpolationType::GetSurfaceNodeIndex(int rSurface, int rNodeIndex) const
{
    assert((unsigned int)rSurface < mSurfaceNodeIndices.size() && "Surface node indices not build.");
    assert((unsigned int)rNodeIndex < mSurfaceNodeIndices[rSurface].size() && "Surface node indices not build.");

    return mSurfaceNodeIndices[rSurface][rNodeIndex];
}

std::string NuTo::InterpolationType::Info() const
{
    std::stringstream out;

    for (auto dofType : mDofs)
    {
        out << Node::DofToString(dofType) << ": " << Get(dofType).GetNumDofs() << "|\t|"
            << "Type and Order: " << Interpolation::TypeOrderToString(Get(dofType).GetTypeOrder()) << "|\n";
    }

    return out.str();
}

void NuTo::InterpolationType::PrintNodeIndices() const
{

    std::cout << " ============ NODE INDICES ============ " << std::endl;
    for (auto dofType : mDofs)
    {
        std::cout << "dof type: " << Node::DofToString(dofType) << " with " << Get(dofType).GetNumNodes()
                  << " nodes:" << std::endl;
        for (int iNode = 0; iNode < Get(dofType).GetNumNodes(); ++iNode)
        {
            std::cout << "Node " << iNode << " = global node " << Get(dofType).GetNodeIndex(iNode) << std::endl;
        }
        std::cout << " ====================================== " << std::endl;
    }
}

void NuTo::InterpolationType::PrintNodeCoordinates() const
{

    std::cout << " ========== NODE COORDINATES ========== " << std::endl;
    for (int iNode = 0; iNode < GetNumNodes(); ++iNode)
    {
        std::cout << "Node " << std::setw(2) << iNode << " at ( ";
        auto nodeCoordinate = GetNaturalNodeCoordinates(iNode);
        for (unsigned int iDim = 0; iDim < nodeCoordinate.size(); ++iDim)
            std::cout << std::setprecision(3) << std::setw(5) << nodeCoordinate[iDim] << " ";
        std::cout << ") with ";
        for (auto dofType : GetNodeDofs(iNode))
            std::cout << Node::DofToString(dofType) << " ";
        std::cout << std::endl;
    }
}

bool NuTo::InterpolationType::CoordinatesAreEqual(const Eigen::VectorXd& rC1, const Eigen::VectorXd& rC2) const
{
    assert(rC1.rows() == rC2.rows());
    if ((rC1 - rC2).norm() > 1.e-10)
        return false;
    return true;
}

const Eigen::MatrixX2i& NuTo::InterpolationType::GetNodeRenumberingIndices() const
{
    return mNodeRenumberingIndices;
}

void NuTo::InterpolationType::UpdateNodeRenumberingIndices()
{
    int numSwaps = 0;
    mNodeRenumberingIndices.resize(0, 2);

    // loop over all points i (with coordinates)
    const InterpolationBase& it = Get(Node::eDof::COORDINATES);

    // why -1? the last check is done for i=num-2 vs j=num-1
    for (int i = 0; i < it.GetNumNodes() - 1; ++i)
    {
        // calculate swapped coordinates x_i --> x_i'
        const Eigen::VectorXd& x_i = it.GetNaturalNodeCoordinates(i);
        Eigen::VectorXd x_i_prime = x_i;

        switch (mShapeType)
        {
        case Interpolation::eShapeType::SPRING:
        case Interpolation::eShapeType::TRUSS1D:
        case Interpolation::eShapeType::TRUSSXD:
        case Interpolation::eShapeType::INTERFACE:
            // reflect at (0,0,0) n = (1,0,0)
            x_i_prime = -x_i;
            break;
        case Interpolation::eShapeType::TRIANGLE2D:
        case Interpolation::eShapeType::QUAD2D:
        case Interpolation::eShapeType::TETRAHEDRON3D:
        case Interpolation::eShapeType::BRICK3D:
            // reflect at (0,0,0) n = (1,-1,0)
            x_i_prime[0] = x_i[1];
            x_i_prime[1] = x_i[0];
            break;
        case Interpolation::eShapeType::PRISM3D:
            // reflect at (0,0,0) n = (1,-1,0) as well...
            x_i_prime[0] = x_i[1];
            x_i_prime[1] = x_i[0];
            break;

        default:
            throw NuTo::MechanicsException(
                    "[NuTo::InterpolationType::UpdateNodeRenumberingIndices] not implemented for " +
                    Interpolation::ShapeTypeToString(mShapeType));
        }

        // find a point j != i with x_j == x_i'

        // why j = i+1? to avoid the case i==j and to avoid duplicated entries [i,j] and [j,i]
        for (int j = i + 1; j < it.GetNumNodes(); ++j)
        {
            const auto& x_j = it.GetNaturalNodeCoordinates(j);
            if (CoordinatesAreEqual(x_i_prime, x_j))
            {
                // store i and j as swap pairs
                mNodeRenumberingIndices.conservativeResize(numSwaps + 1, 2);
                mNodeRenumberingIndices(numSwaps, 0) = i;
                mNodeRenumberingIndices(numSwaps, 1) = j;
                numSwaps++;
            }
        }
    }
}
