/*
 * InterpolationType.cpp
 *
 *  Created on: 31 Mar 2015
 *      Author: ttitsche
 */

#include "nuto/mechanics/interpolationtypes/InterpolationType.h"

#include "nuto/mechanics/interpolationtypes/Interpolation2DTriangle.h"
#include "nuto/mechanics/interpolationtypes/Interpolation2DQuad.h"
#include "nuto/mechanics/interpolationtypes/Interpolation3DTetrahedron.h"
#include "nuto/mechanics/interpolationtypes/Interpolation3DBrick.h"
#include "nuto/mechanics/interpolationtypes/Interpolation1DTruss.h"

#include <boost/foreach.hpp>

#include <iomanip>


NuTo::InterpolationType::InterpolationType(const StructureBase* rStructure, NuTo::Interpolation::eShapeType rShapeType) :
        mShapeType(rShapeType), mNumDofs(0), mNumActiveDofs(0), mIntegrationType(nullptr), mStructure(rStructure)
{
}

NuTo::InterpolationType::~InterpolationType()
{
}

const NuTo::InterpolationBase& NuTo::InterpolationType::Get(const Node::eAttributes& rDofType) const
{
    try
    {
        return mInterpolations.at(rDofType);
    } catch (boost::bad_ptr_container_operation& e)
    {
        std::cout << e.what() << std::endl;
        throw NuTo::MechanicsException("[NuTo::InterpolationType::Get] Dof " + Node::AttributeToString(rDofType) + " is not a member of this interpolation type. Add it first.");
    }
}

NuTo::InterpolationBase& NuTo::InterpolationType::GetNonConst(Node::eAttributes rDofType)
{
    auto interpolationTypeIterator = mInterpolations.find(rDofType);

    if (interpolationTypeIterator == mInterpolations.end())
        throw NuTo::MechanicsException("[NuTo::InterpolationType::Get] Dof " + Node::AttributeToString(rDofType) + " is not a member of this interpolation type. Add it first.");

    return *(interpolationTypeIterator->second);
}

void NuTo::InterpolationType::AddDofInterpolation(Node::eAttributes rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder)
{
    if (IsDof(rDofType))
        throw NuTo::MechanicsException("[NuTo::InterpolationTypeBase::AddDofInterpolation] Dof " + NuTo::Node::AttributeToString(rDofType) + " exists.");

    InterpolationBase* newType;
    switch (mShapeType)
    {
    case Interpolation::eShapeType::TRUSS1D:
        newType = new Interpolation1DTruss(mStructure, rDofType, rTypeOrder);
        break;
    case Interpolation::eShapeType::TRIANGLE2D:
        newType = new Interpolation2DTriangle(mStructure, rDofType, rTypeOrder);
        break;
    case Interpolation::eShapeType::QUAD2D:
        newType = new Interpolation2DQuad(mStructure, rDofType, rTypeOrder);
        break;
    case Interpolation::eShapeType::TETRAHEDRON3D:
        newType = new Interpolation3DTetrahedron(mStructure, rDofType, rTypeOrder);
        break;
    case Interpolation::eShapeType::BRICK3D:
        newType = new Interpolation3DBrick(mStructure, rDofType, rTypeOrder);
        break;
    default:
        throw NuTo::MechanicsException("[NuTo::InterpolationType::AddDofInterpolation] ShapeType " + NuTo::Interpolation::ShapeTypeToString(mShapeType) + " not implemented.");
    }

    mInterpolations.insert(rDofType, newType);
    mDofs.insert(rDofType);

    mNumDofs += newType->GetNumDofs();

    if (rDofType == Node::COORDINATES)
    {
        SetIsActive(false, rDofType);
        SetIsConstitutiveInput(false, rDofType);
    } else
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

    UpdateLocalStartIndices();

    if (rDofType == Node::COORDINATES)
        UpdateNodeRenumberingIndices();

    if (mIntegrationType != nullptr)
        newType->UpdateIntegrationType(*mIntegrationType);
}

void NuTo::InterpolationType::UpdateIntegrationType(const IntegrationTypeBase& rIntegrationType)
{
    if (mIntegrationType == &rIntegrationType)
        return;

    for (auto dofType : mDofs)
        GetNonConst(dofType).UpdateIntegrationType(rIntegrationType);

    mIntegrationType = &rIntegrationType;
}

//! @brief returns the pointer to the integration type that is currently used
const NuTo::IntegrationTypeBase* NuTo::InterpolationType::GetCurrentIntegrationType() const
{
    return mIntegrationType;
}

const NuTo::Interpolation::eShapeType NuTo::InterpolationType::GetShapeType() const
{
    return mShapeType;
}

NuTo::IntegrationType::eIntegrationType NuTo::InterpolationType::GetStandardIntegrationType() const
{
    NuTo::IntegrationType::eIntegrationType integrationType = IntegrationType::eIntegrationType::NotSet;
    int maxOrder = 0;

    for (Node::eAttributes dof : GetDofs())
    {
        Interpolation::eTypeOrder order = Get(dof).GetTypeOrder();
        int currentOrder = 0;
        switch (order)
        {
        case Interpolation::EQUIDISTANT1:
            currentOrder = 1;
            break;
        case Interpolation::LOBATTO2:
        case Interpolation::EQUIDISTANT2:
            currentOrder = 2;
            break;
        case Interpolation::LOBATTO3:
        case Interpolation::EQUIDISTANT3:
            currentOrder = 3;
            break;
        case Interpolation::LOBATTO4:
        case Interpolation::EQUIDISTANT4:
            currentOrder = 4;
            break;
        default:
            throw NuTo::MechanicsException("[NuTo::InterpolationType::GetStandardIntegrationType] Standard integration type for " + Interpolation::TypeOrderToString(order) + " is not implemented.");
        }

        if (currentOrder > maxOrder)
        {
            maxOrder = currentOrder;
            integrationType = Get(dof).GetStandardIntegrationType();
        }
    }

    return integrationType;
}

const Eigen::VectorXi NuTo::InterpolationType::GetSurfaceNodeIndices(int rSurface) const
{
    assert(IsDof(Node::COORDINATES));
    const InterpolationBase& interpolationType = Get(Node::COORDINATES);

    const auto& surfaceEdgesCoordinates = interpolationType.GetSurfaceEdgesCoordinates(rSurface);

    Eigen::VectorXi surfaceNodeIndices(surfaceEdgesCoordinates.size());

    for (unsigned int i = 0; i < surfaceEdgesCoordinates.size(); ++i)
    {
        const Eigen::VectorXd& naturalSurfaceEdgeCoordinate = surfaceEdgesCoordinates[i];
//        std::cout << naturalSurfaceEdgeCoordinate << std::endl;
        Eigen::VectorXd shapeFunctions = interpolationType.CalculateShapeFunctions(naturalSurfaceEdgeCoordinate);
//        std::cout << shapeFunctions << std::endl;
        assert(std::abs(shapeFunctions.norm() - 1) < 1.e-8);
        int indexNodeCoordinate;                                        // index in the coordinate interpolation

#ifdef DEBUG
        double value = shapeFunctions.maxCoeff(&indexNodeCoordinate);   // find the index where the shape function is 1
        assert(std::abs(value - 1) < 1.e-8);
#else
        shapeFunctions.maxCoeff(&indexNodeCoordinate);   // find the index where the shape function is 1
#endif
        surfaceNodeIndices(i) = interpolationType.GetNodeIndex(indexNodeCoordinate); // index in the Element::mNodes vector
    }

    return surfaceNodeIndices;
}

int NuTo::InterpolationType::GetNumSurfaces() const
{
    assert(mInterpolations.size() != 0);
    return mInterpolations.begin()->second->GetNumSurfaces();
}

void NuTo::InterpolationType::UpdateLocalStartIndices()
{
    // calculate local start indices
    // prescribe a specific order

    std::vector<Node::eAttributes> orderedDofs(
    { Node::COORDINATES, Node::DISPLACEMENTS, Node::TEMPERATURES, Node::NONLOCALEQPLASTICSTRAIN, Node::NONLOCALEQSTRAIN, Node::RELATIVEHUMIDITY, Node::WATERVOLUMEFRACTION });

    int currentStartIndex = 0;
    for (unsigned int i = 0; i < orderedDofs.size(); ++i)
    {
        auto dof = orderedDofs[i];
        if (IsDof(dof))
        {
            GetNonConst(dof).mLocalStartIndex = currentStartIndex;
            if (mActiveDofs.find(dof) != mActiveDofs.end())
                currentStartIndex += GetNonConst(dof).GetNumDofs();
        }
    }
}

void NuTo::InterpolationType::SetIsActive(bool rIsActiveDof, Node::eAttributes rDofType)
{
    if (not IsDof(rDofType))
        throw NuTo::MechanicsException("[NuTo::InterpolationType::SetIsActive] Dof " + Node::AttributeToString(rDofType) + " is not a member of this interpolation type. Add it first.");

    if (rIsActiveDof)
        mActiveDofs.insert(rDofType);
    else
        mActiveDofs.erase(rDofType);

    GetNonConst(rDofType).mIsActive = rIsActiveDof;

    UpdateLocalStartIndices();

    mNumActiveDofs = 0;
    for (auto dof : mActiveDofs)
    {
        mNumActiveDofs += Get(dof).GetNumDofs();
    }

}

bool NuTo::InterpolationType::IsActive(const Node::eAttributes& rDofType) const
{
    return Get(rDofType).IsActive();
}

void NuTo::InterpolationType::SetIsConstitutiveInput(bool rIsConstitutiveInput, Node::eAttributes rDofType)
{
    GetNonConst(rDofType).mIsConstitutiveInput = rIsConstitutiveInput;
}

bool NuTo::InterpolationType::IsConstitutiveInput(const Node::eAttributes& rDofType) const
{
    if (not IsDof(rDofType))
        return false;
    return Get(rDofType).IsConstitutiveInput();
}

bool NuTo::InterpolationType::IsDof(const Node::eAttributes& rDofType) const
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

const std::set<NuTo::Node::eAttributes>& NuTo::InterpolationType::GetActiveDofs() const
{
    return mActiveDofs;
}

const std::set<NuTo::Node::eAttributes>& NuTo::InterpolationType::GetDofs() const
{
    return mDofs;
}

std::set<NuTo::Node::eAttributes> NuTo::InterpolationType::GetNodeDofs(int rNodeIndex) const
{
    assert((unsigned int )rNodeIndex < mNodeDofs.size());
    return mNodeDofs[rNodeIndex];
}

const Eigen::VectorXd& NuTo::InterpolationType::GetNaturalNodeCoordinates(int rNodeIndex) const
{
    assert((unsigned int )rNodeIndex < mNodeCoordinates.size());
    return mNodeCoordinates[rNodeIndex];
}
int NuTo::InterpolationType::GetNumNodes() const
{
    return mNodeCoordinates.size();
}

std::string NuTo::InterpolationType::Info() const
{
    std::stringstream out;

    for (auto dofType : mDofs)
    {
        out << Node::AttributeToString(dofType) << ": " << Get(dofType).GetNumDofs() << "|\t|" << "Type and Order: " << Interpolation::TypeOrderToString(Get(dofType).GetTypeOrder()) << "|\n";
    }

    return out.str();
}

void NuTo::InterpolationType::PrintNodeIndices() const
{

    std::cout << " ============ NODE INDICES ============ " << std::endl;
    for (auto dofType : mDofs)
    {
        std::cout << "dof type: " << Node::AttributeToString(dofType) << " with " << Get(dofType).GetNumNodes() << " nodes:" << std::endl;
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
            std::cout << Node::AttributeToString(dofType) << " ";
        std::cout << std::endl;
    }

}

bool NuTo::InterpolationType::CoordinatesAreEqual(const Eigen::VectorXd& rC1, const Eigen::VectorXd& rC2) const
{
    assert(rC1.rows() == rC2.rows());
    for (unsigned int iDim = 0; iDim < rC1.rows(); ++iDim)
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
    const InterpolationBase& it = Get(Node::COORDINATES);

    // why -1? the last check is done for i=num-2 vs j=num-1
    for (int i = 0; i < it.GetNumNodes() - 1; ++i)
    {
        // calculate swapped coordinates x_i --> x_i'
        const Eigen::VectorXd& x_i = it.GetNaturalNodeCoordinates(i);
        Eigen::VectorXd x_i_prime = x_i;

        switch (mShapeType)
        {
        case Interpolation::TRUSS1D:
        case Interpolation::TRUSSXD:
            // reflect at (0,0,0) n = (1,0,0)
            x_i_prime = -x_i;
            break;
        case Interpolation::TRIANGLE2D:
        case Interpolation::QUAD2D:
        case Interpolation::TETRAHEDRON3D:
        case Interpolation::BRICK3D:
            // reflect at (0,0,0) n = (1,-1,0)
            x_i_prime[0] = x_i[1];
            x_i_prime[1] = x_i[0];
            break;
        default:
            throw NuTo::MechanicsException("[NuTo::InterpolationType::UpdateNodeRenumberingIndices] not implemented for " + Interpolation::ShapeTypeToString(mShapeType));
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
