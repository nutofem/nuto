#include <cassert>
#include <typeinfo>

#include "base/Timer.h"

#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/elements/ContinuumElement.h"
#include "mechanics/elements/ContinuumElementIGA.h"
#include "mechanics/elements/ContinuumBoundaryElement.h"
#include "mechanics/elements/ContinuumBoundaryElementConstrainedControlNode.h"
#include "mechanics/elements/Element1DInXD.h"
#include "mechanics/elements/Element2DInterface.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "mechanics/interpolationtypes/InterpolationBase.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"

#include "mechanics/groups/Group.h"
#include "mechanics/groups/GroupEnum.h"

#include "mechanics/mesh/MeshCompanion.h"

#include <eigen3/Eigen/Core>

//! @brief returns the number of nodes
//! @return number of nodes
int NuTo::Structure::GetNumElements() const
{
    return mElementMap.size();
}

//! @brief a reference to an element
//! @param identifier
//! @return reference to an element
NuTo::ElementBase* NuTo::Structure::ElementGetElementPtr(int rIdent)
{
    boost::ptr_map<int, ElementBase>::iterator it = mElementMap.find(rIdent);
    if (it != mElementMap.end())
        return it->second;
    else
    {
        throw Exception(__PRETTY_FUNCTION__,
                                 "Element with identifier " + std::to_string(rIdent) + " does not exist.");
    }
}

//! @brief a reference to an element
//! @param identifier
//! @return reference to an element
const NuTo::ElementBase* NuTo::Structure::ElementGetElementPtr(int rIdent) const
{
    boost::ptr_map<int, ElementBase>::const_iterator it = mElementMap.find(rIdent);
    if (it != mElementMap.end())
        return it->second;
    else
    {
        throw Exception(__PRETTY_FUNCTION__,
                                 "Element with identifier " + std::to_string(rIdent) + " does not exist.");
    }
}

//! @brief gives the identifier of an element
//! @param reference to an element
//! @return element number
int NuTo::Structure::ElementGetId(const ElementBase* rElement) const
{
    for (boost::ptr_map<int, ElementBase>::const_iterator it = mElementMap.begin(); it != mElementMap.end(); it++)
    {
        if (it->second == rElement)
            return it->first;
    }
    throw Exception(__PRETTY_FUNCTION__, "Element does not exist.");
}

//! @brief returns a vector with the node ids of an element
//! @param identifier
//! @return vector with node ids
std::vector<int> NuTo::Structure::ElementGetNodes(int rId)
{
    NuTo::ElementBase* elementPtr = ElementGetElementPtr(rId);
    std::vector<int> nodeVector(elementPtr->GetNumNodes());
    for (int count = 0; count < elementPtr->GetNumNodes(); count++)
        nodeVector[count] = this->NodeGetId(elementPtr->GetNode(count));
    return nodeVector;
}

//! @brief info about one single element
//! @param rElement (Input) ... pointer to the element
//! @param rVerboseLevel (Input) ... level of verbosity
void NuTo::Structure::ElementInfo(const ElementBase* rElement, int rVerboseLevel) const
{
    if (rVerboseLevel > 2)
    {
        if (rVerboseLevel > 3)
        {
            std::cout << "\tNodes:";
            for (unsigned short iNode = 0; iNode < rElement->GetNumNodes(); ++iNode)
                std::cout << "\t" << this->NodeGetId(rElement->GetNode(iNode));
            std::cout << std::endl;
            if (rVerboseLevel > 4)
            {
                std::cout << "\tintegration points :" << std::endl;
                for (int iIp = 0; iIp < rElement->GetNumIntegrationPoints(); ++iIp)
                {
                    Eigen::Vector3d coor = rElement->GetGlobalIntegrationPointCoordinates(iIp);
                    std::cout << "\t\t" << iIp << ": [" << coor(0, 0) << ";" << coor(1, 0) << ";" << coor(2, 0) << "]"
                              << std::endl;
                }
            }
        }
    }
}

//! @brief info about the elements in the Structure
void NuTo::Structure::ElementInfo(int rVerboseLevel) const
{
    mLogger << "number of elements: " << mElementMap.size() << "\n";
    if (rVerboseLevel > 3)
    {
        mLogger << "\t\telements :"
                << "\n";
        for (boost::ptr_map<int, ElementBase>::const_iterator it = mElementMap.begin(); it != mElementMap.end(); it++)
        {
            mLogger << "\t\t" << it->first;
            if (rVerboseLevel > 4)
            {
                mLogger << "\t:";
                for (unsigned short iNode = 0; iNode < it->second->GetNumNodes(); ++iNode)
                    mLogger << "\t" << this->NodeGetId(it->second->GetNode(iNode));
            }
            mLogger << "\n";
        }
    }
}


void NuTo::Structure::ElementTotalSetInterpolationType(const int rInterpolationTypeId)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, InterpolationType>::iterator itInterpolationType =
            mInterpolationTypeMap.find(rInterpolationTypeId);
    if (itInterpolationType == mInterpolationTypeMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Interpolation type with the given identifier does not exist.");

    for (const auto& elementPair : mElementMap)
        ElementSetInterpolationType(elementPair.second, itInterpolationType->second);
}

void NuTo::Structure::ElementTotalConvertToInterpolationType()
{
    MeshCompanion::ElementTotalConvertToInterpolationType(*this);
}

void NuTo::Structure::ElementConvertToInterpolationType(int rGroupNumberElements)
{
    MeshCompanion::ElementConvertToInterpolationType(*this, rGroupNumberElements);
}

void NuTo::Structure::ElementTotalConvertToInterpolationType(double rNodeDistanceMerge, double rMeshSize)
{
    (void)rMeshSize; // unsused
    MeshCompanion::ElementTotalConvertToInterpolationType(*this, rNodeDistanceMerge);
}

void NuTo::Structure::ElementConvertToInterpolationType(int rGroupNumberElements, double rNodeDistanceMerge,
                                                        double rMeshSize)
{
    (void)rMeshSize; // unsused
    MeshCompanion::ElementConvertToInterpolationType(*this, rGroupNumberElements, rNodeDistanceMerge);
}


void NuTo::Structure::ElementDelete(int rElementNumber)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());
    ElementDeleteInternal(rElementNumber);
}


int NuTo::Structure::ElementCreate(int rInterpolationTypeId, const Eigen::VectorXi& rNodeNumbers,
                                   const Eigen::MatrixXd& rKnots, const Eigen::VectorXi& rKnotIDs)
{
    int elementNumber = GetUnusedId(mElementMap);
    this->ElementCreate(elementNumber, rInterpolationTypeId, rNodeNumbers, rKnots, rKnotIDs);
    return elementNumber;
}

int NuTo::Structure::ElementCreate(int rInterpolationTypeId, const std::vector<int>& rNodeNumbers)
{
    int elementNumber = GetUnusedId(mElementMap);
    this->ElementCreate(elementNumber, rInterpolationTypeId, rNodeNumbers);
    return elementNumber;
}

int NuTo::Structure::ElementCreate(int rInterpolationTypeId, std::vector<NodeBase*> rNodes)
{
    int elementNumber = GetUnusedId(mElementMap);
    this->ElementCreate(elementNumber, rInterpolationTypeId, rNodes);
    return elementNumber;
}


void NuTo::Structure::ElementCreate(int rElementNumber, int rInterpolationTypeId, const Eigen::VectorXi& rNodeNumbers,
                                    const Eigen::MatrixXd& rKnots, const Eigen::VectorXi& rKnotIDs)
{
    // convert node numbers to pointer
    std::vector<NodeBase*> nodeVector;
    for (int iNode = 0; iNode < rNodeNumbers.rows(); iNode++)
        nodeVector.push_back(NodeGetNodePtr(rNodeNumbers(iNode)));

    boost::ptr_map<int, InterpolationType>::iterator itIterator = mInterpolationTypeMap.find(rInterpolationTypeId);

    if (itIterator == mInterpolationTypeMap.end())
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Interpolation type does not exist.");

    InterpolationType& interpolationType = *itIterator->second;

    if (not interpolationType.IsDof(Node::eDof::COORDINATES))
        throw NuTo::Exception(__PRETTY_FUNCTION__, "COORDINATE interpolation required.");

    unsigned int numNodesCoordinates = interpolationType.Get(Node::eDof::COORDINATES).GetNumNodes();
    if (numNodesCoordinates != nodeVector.size())
        throw NuTo::Exception(__PRETTY_FUNCTION__,
                                       "COORDINATE interpolation requires " + std::to_string(numNodesCoordinates) +
                                               " nodes. " + std::to_string(nodeVector.size()) + " are provided.");

    interpolationType.ClearCache();
    const auto& integrationType = *GetPtrIntegrationType(interpolationType.GetStandardIntegrationType());

    ElementBase* ptrElement = nullptr;
    switch (interpolationType.GetShapeType())
    {
    case NuTo::Interpolation::eShapeType::SPRING:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Element1DSpring currently not implemented.");
        break;
    case NuTo::Interpolation::eShapeType::TRUSS1D:
    case NuTo::Interpolation::eShapeType::TRUSSXD:
    case NuTo::Interpolation::eShapeType::TRIANGLE2D:
    case NuTo::Interpolation::eShapeType::QUAD2D:
    case NuTo::Interpolation::eShapeType::TETRAHEDRON3D:
    case NuTo::Interpolation::eShapeType::BRICK3D:
    case NuTo::Interpolation::eShapeType::INTERFACE:
        throw NuTo::Exception(
                __PRETTY_FUNCTION__,
                "Please use approriate functions for element creation, this is IGA implementation.");
        break;
    case NuTo::Interpolation::eShapeType::IGA1D:
        ptrElement = new ContinuumElementIGA<1>(nodeVector, rKnots, rKnotIDs, interpolationType, integrationType,
                                                GetDofStatus());
        ptrElement->CheckElement();
        break;
    case NuTo::Interpolation::eShapeType::IGA2D:
        ptrElement = new ContinuumElementIGA<2>(nodeVector, rKnots, rKnotIDs, interpolationType, integrationType,
                                                GetDofStatus());
        ptrElement->CheckElement();
        break;
    default:
        throw Exception(__PRETTY_FUNCTION__, "invalid dimension.");
    }

    mElementMap.insert(rElementNumber, ptrElement);
}

void NuTo::Structure::ElementCreate(int elementNumber, int interpolationTypeId, const std::vector<int>& rNodeNumbers)
{
    // convert node numbers to pointer
    std::vector<NodeBase*> nodeVector;
    for (auto nodeNumber : rNodeNumbers)
        nodeVector.push_back(NodeGetNodePtr(nodeNumber));

    ElementCreate(elementNumber, interpolationTypeId, nodeVector);
}

void NuTo::Structure::ElementCreate(int rElementNumber, int rInterpolationTypeId, std::vector<NodeBase*> rNodes)
{
    boost::ptr_map<int, InterpolationType>::iterator itIterator = mInterpolationTypeMap.find(rInterpolationTypeId);

    if (itIterator == mInterpolationTypeMap.end())
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Interpolation type does not exist.");

    InterpolationType& interpolationType = *itIterator->second;

    if (not interpolationType.IsDof(Node::eDof::COORDINATES))
        throw NuTo::Exception(__PRETTY_FUNCTION__, "COORDINATE interpolation required.");

    unsigned int numNodesCoordinates = interpolationType.Get(Node::eDof::COORDINATES).GetNumNodes();
    if (numNodesCoordinates != rNodes.size())
        throw NuTo::Exception(__PRETTY_FUNCTION__, "COORDINATE interpolation requires " +
                                                                    std::to_string(numNodesCoordinates) + " nodes. " +
                                                                    std::to_string(rNodes.size()) + " are provided.");

    interpolationType.ClearCache();
    const auto& integrationType = *GetPtrIntegrationType(interpolationType.GetStandardIntegrationType());

    ElementBase* ptrElement = nullptr;
    switch (interpolationType.GetShapeType())
    {
    case NuTo::Interpolation::eShapeType::SPRING:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Element1DSpring currently not implemented.");
        break;
    case NuTo::Interpolation::eShapeType::TRUSS1D:
        ptrElement = new ContinuumElement<1>(rNodes, interpolationType, integrationType, GetDofStatus());
        ptrElement->CheckElement();
        break;
    case NuTo::Interpolation::eShapeType::TRUSSXD:
        ptrElement = new Element1DInXD(rNodes, interpolationType, integrationType, GetDofStatus(), mDimension);
        ptrElement->CheckElement();
        break;
    case NuTo::Interpolation::eShapeType::TRIANGLE2D:
    case NuTo::Interpolation::eShapeType::QUAD2D:
        ptrElement = new ContinuumElement<2>(rNodes, interpolationType, integrationType, GetDofStatus());
        ptrElement->CheckElement();
        break;
    case NuTo::Interpolation::eShapeType::TETRAHEDRON3D:
    case NuTo::Interpolation::eShapeType::BRICK3D:
    case NuTo::Interpolation::eShapeType::PRISM3D:
        ptrElement = new ContinuumElement<3>(rNodes, interpolationType, integrationType, GetDofStatus());
        ptrElement->CheckElement();
        break;
    case NuTo::Interpolation::eShapeType::INTERFACE:
        ptrElement = new Element2DInterface(rNodes, interpolationType, integrationType, mDimension);
        ptrElement->CheckElement();
        break;
    case NuTo::Interpolation::eShapeType::IGA1D:
        throw NuTo::Exception(
                __PRETTY_FUNCTION__,
                "Please use the ElementCreate function for IGA elements, where the knot parameters are provided.");
        break;
    default:
        throw Exception(__PRETTY_FUNCTION__, "invalid dimension.");
    }

    mElementMap.insert(rElementNumber, ptrElement);
}


int NuTo::Structure::ElementsCreate(int rInterpolationTypeId, const Eigen::MatrixXi& rNodeNumbers)
{
    std::vector<int> newElementIds;
    // go through the elements
    for (int iNode = 0; iNode < rNodeNumbers.cols(); ++iNode)
    {
        auto column = rNodeNumbers.col(iNode);
        std::vector<int> incidence(column.data(), column.data() + column.size());
        int newElementId = ElementCreate(rInterpolationTypeId, incidence);
        newElementIds.push_back(newElementId);
    }

    bool showTime = mShowTime;
    mShowTime = false;

    // create element group containing the new elements
    int newElementGroup = GroupCreate(eGroupId::Elements);
    for (int newElementId : newElementIds)
        GroupAddElement(newElementGroup, newElementId);

    mShowTime = showTime;
    return newElementGroup;
}


int NuTo::Structure::BoundaryElementsCreate(int rElementGroupId, int rNodeGroupId, NodeBase* rControlNode)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    // find groups
    boost::ptr_map<int, GroupBase>::iterator itGroupElements = mGroupMap.find(rElementGroupId);
    if (itGroupElements == mGroupMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroupElements->second->GetType() != NuTo::eGroupId::Elements)
        throw Exception(__PRETTY_FUNCTION__, "Group is not an element group.");

    boost::ptr_map<int, GroupBase>::iterator itGroupBoundaryNodes = mGroupMap.find(rNodeGroupId);
    if (itGroupBoundaryNodes == mGroupMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroupBoundaryNodes->second->GetType() != NuTo::eGroupId::Nodes)
        throw Exception(__PRETTY_FUNCTION__, "Group is not a node group.");

    Group<ElementBase>& elementGroup = *(itGroupElements->second->AsGroupElement());
    Group<NodeBase>& nodeGroup = *(itGroupBoundaryNodes->second->AsGroupNode());

    // since the search is done via the id's, the surface nodes are ptr, so make another set with the node ptrs
    std::set<const NodeBase*> nodePtrSet;
    for (auto itNode : nodeGroup)
    {
        nodePtrSet.insert(itNode.second);
    }

    std::vector<int> newBoundaryElementIds;

    // loop over all elements
    for (auto itElement : elementGroup)
    {
        ElementBase* elementPtr = itElement.second;
        const InterpolationType& interpolationType = elementPtr->GetInterpolationType();

        // std::cout << typeid(*elementPtr).name() << "\n";
        // std::cout << typeid(ContinuumElement<1>).name() << std::endl;
        if (typeid(*elementPtr) != typeid(ContinuumElement<1>) && typeid(*elementPtr) != typeid(ContinuumElement<2>) &&
            typeid(*elementPtr) != typeid(ContinuumElement<3>))
            throw Exception(__PRETTY_FUNCTION__, "Element is not a ContinuumElement.");

        // loop over all surfaces
        for (int iSurface = 0; iSurface < interpolationType.GetNumSurfaces(); ++iSurface)
        {
            bool elementSurfaceNodesMatchBoundaryNodes = true;
            Eigen::VectorXi surfaceNodeIndices = interpolationType.GetSurfaceNodeIndices(iSurface);

            int numSurfaceNodes = surfaceNodeIndices.rows();
            std::vector<const NodeBase*> surfaceNodes(numSurfaceNodes);

            for (int iSurfaceNode = 0; iSurfaceNode < numSurfaceNodes; ++iSurfaceNode)
            {
                surfaceNodes[iSurfaceNode] = elementPtr->GetNode(surfaceNodeIndices(iSurfaceNode, 0));
            }

            // check, if all surface nodes are in the node group
            for (auto& surfaceNode : surfaceNodes)
            {
                if (nodePtrSet.find(surfaceNode) == nodePtrSet.end())
                {
                    // this surface has at least one node that is not in the list, continue
                    elementSurfaceNodesMatchBoundaryNodes = false;
                    break;
                }
            }
            if (not elementSurfaceNodesMatchBoundaryNodes)
                continue;

            int surfaceId = iSurface;

            int elementId = GetUnusedId(mElementMap);

            ElementBase* boundaryElement = nullptr;
            ConstitutiveBase& constitutiveLaw = elementPtr->GetConstitutiveLaw(0);

            switch (elementPtr->GetLocalDimension())
            {
            case 1:
            {
                const auto& integrationType = *GetPtrIntegrationType(eIntegrationType::IntegrationType0DBoundary);
                auto& element = dynamic_cast<ContinuumElement<1>&>(*elementPtr);
                if (rControlNode == nullptr)
                    boundaryElement = new ContinuumBoundaryElement<1>(element, integrationType, surfaceId);
                else
                    boundaryElement = new ContinuumBoundaryElementConstrainedControlNode<1>(element, integrationType,
                                                                                            surfaceId, rControlNode);
                break;
            }
            case 2:
            {
                eIntegrationType integrationTypeEnum;
                // check for 2D types
                switch (interpolationType.GetStandardIntegrationType())
                {
                case eIntegrationType::IntegrationType2D3NGauss1Ip:
                case eIntegrationType::IntegrationType2D4NGauss1Ip:
                    integrationTypeEnum = eIntegrationType::IntegrationType1D2NGauss1Ip;
                    break;
                case eIntegrationType::IntegrationType2D3NGauss3Ip:
                case eIntegrationType::IntegrationType2D4NGauss4Ip:
                    integrationTypeEnum = eIntegrationType::IntegrationType1D2NGauss2Ip;
                    break;
                case eIntegrationType::IntegrationType2D3NGauss6Ip:
                case eIntegrationType::IntegrationType2D4NGauss9Ip:
                    integrationTypeEnum = eIntegrationType::IntegrationType1D2NGauss3Ip;
                    break;
                case eIntegrationType::IntegrationType2D3NGauss12Ip:
                    integrationTypeEnum = eIntegrationType::IntegrationType1D2NGauss5Ip;
                    break;

                case eIntegrationType::IntegrationType2D4NLobatto9Ip:
                    integrationTypeEnum = eIntegrationType::IntegrationType1D2NLobatto3Ip;
                    break;

                case eIntegrationType::IntegrationType2D4NLobatto16Ip:
                    integrationTypeEnum = eIntegrationType::IntegrationType1D2NLobatto4Ip;
                    break;

                case eIntegrationType::IntegrationType2D4NLobatto25Ip:
                    integrationTypeEnum = eIntegrationType::IntegrationType1D2NLobatto5Ip;
                    break;

                default:
                    throw Exception(
                            __PRETTY_FUNCTION__,
                            "Could not automatically determine integration type of the boundary element.");
                }
                const auto& integrationType = *GetPtrIntegrationType(integrationTypeEnum);
                auto& element = dynamic_cast<ContinuumElement<2>&>(*elementPtr);
                if (rControlNode == nullptr)
                    boundaryElement = new ContinuumBoundaryElement<2>(element, integrationType, surfaceId);
                else
                    boundaryElement = new ContinuumBoundaryElementConstrainedControlNode<2>(element, integrationType,
                                                                                            surfaceId, rControlNode);

                break;
            }
            case 3:
            {
                eIntegrationType integrationTypeEnum;
                // check for 3D types
                switch (interpolationType.GetStandardIntegrationType())
                {
                case eIntegrationType::IntegrationType3D4NGauss1Ip:
                    integrationTypeEnum = eIntegrationType::IntegrationType2D3NGauss1Ip;
                    break;

                case eIntegrationType::IntegrationType3D4NGauss4Ip:
                    integrationTypeEnum = eIntegrationType::IntegrationType2D3NGauss3Ip;
                    break;

                case eIntegrationType::IntegrationType3D8NGauss1Ip:
                    integrationTypeEnum = eIntegrationType::IntegrationType2D4NGauss1Ip;
                    break;

                case eIntegrationType::IntegrationType3D8NGauss2x2x2Ip:
                    integrationTypeEnum = eIntegrationType::IntegrationType2D4NGauss4Ip;
                    break;

                case eIntegrationType::IntegrationType3D8NLobatto3x3x3Ip:
                    integrationTypeEnum = eIntegrationType::IntegrationType2D4NLobatto9Ip;
                    break;

                case eIntegrationType::IntegrationType3D8NLobatto4x4x4Ip:
                    integrationTypeEnum = eIntegrationType::IntegrationType2D4NLobatto16Ip;
                    break;

                case eIntegrationType::IntegrationType3D8NLobatto5x5x5Ip:
                    integrationTypeEnum = eIntegrationType::IntegrationType2D4NLobatto16Ip;
                    break;

                default:
                    throw Exception(
                            __PRETTY_FUNCTION__,
                            "Could not automatically determine integration type of the boundary element.");
                }

                const auto& integrationType = *GetPtrIntegrationType(integrationTypeEnum);
                auto& element = dynamic_cast<ContinuumElement<3>&>(*elementPtr);
                if (rControlNode == nullptr)
                    boundaryElement = new ContinuumBoundaryElement<3>(element, integrationType, surfaceId);
                else
                    boundaryElement = new ContinuumBoundaryElementConstrainedControlNode<3>(element, integrationType,
                                                                                            surfaceId, rControlNode);
                break;
            }
            default:
                throw Exception(__PRETTY_FUNCTION__, "Boundary element for Continuum element with dimension " +
                                                                      std::to_string(elementPtr->GetLocalDimension()) +
                                                                      "not implemented");
            }

            mElementMap.insert(elementId, boundaryElement);
            newBoundaryElementIds.push_back(elementId);

            boundaryElement->SetConstitutiveLaw(constitutiveLaw);
        }
    }

    int boundaryElementGroup = GroupCreate(eGroupId::Elements);
    for (int boundaryElementId : newBoundaryElementIds)
        GroupAddElement(boundaryElementGroup, boundaryElementId);

    return boundaryElementGroup;
}

std::pair<int, int> NuTo::Structure::InterfaceElementsCreate(int rElementGroupId, int rInterfaceInterpolationType,
                                                             int rFibreInterpolationType)
{

    // find element group
    boost::ptr_map<int, GroupBase>::iterator itGroupElements = mGroupMap.find(rElementGroupId);
    if (itGroupElements == mGroupMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroupElements->second->GetType() != NuTo::eGroupId::Elements)
        throw Exception(__PRETTY_FUNCTION__, "Group is not an element group.");


    // gets member ids from an element group. The element group must only contain truss elements
    auto elementIds = GroupGetMemberIds(rElementGroupId);

    int groupElementsInterface = GroupCreate(NuTo::eGroupId::Elements);
    int groupElementsFibre = GroupCreate(NuTo::eGroupId::Elements);

    // loop over elements in element group
    for (int elementId : elementIds)
    {
        auto nodeIds = ElementGetNodes(elementId);

        assert((nodeIds.size() == 2 or nodeIds.size() == 3) and
               "Only implemented for the 4 node and 6 node interface element");

        std::vector<int> nodeIdsFibre(nodeIds.size());
        std::vector<int> nodeIdsMatrix(nodeIds.size());

        // loop over nodes of element
        for (unsigned int k = 0; k < nodeIds.size(); ++k)
        {
            Eigen::VectorXd nodeCoordinates;
            NodeGetCoordinates(nodeIds[k], nodeCoordinates);

            int groupNodes = GroupCreate(NuTo::eGroupId::Nodes);
            GroupAddNodeRadiusRange(groupNodes, nodeCoordinates, 0.0, 1e-6);

            // create an additional node at the same position if it has not been created already
            if (GroupGetNumMembers(groupNodes) > 1)
            {
                assert(GroupGetNumMembers(groupNodes) == 2 and
                       "This group should have exactly two members. Check what went wrong!");
                auto groupNodeMemberIds = GroupGetMemberIds(groupNodes);
                if (groupNodeMemberIds[0] == nodeIds[k])
                {
                    nodeIdsFibre[k] = groupNodeMemberIds[1];
                }
                else
                {
                    nodeIdsFibre[k] = groupNodeMemberIds[0];
                }
            }
            else
            {
                std::set<NuTo::Node::eDof> dofs;
                dofs.insert(NuTo::Node::eDof::COORDINATES);
                dofs.insert(NuTo::Node::eDof::DISPLACEMENTS);

                nodeIdsFibre[k] = NodeCreate(nodeCoordinates, dofs);
            }

            nodeIdsMatrix[k] = nodeIds[k];
        }

        // create interface element
        std::vector<int> nodeIndicesInterface(2 * nodeIds.size());
        for (unsigned int iIndex = 0; iIndex < nodeIds.size(); ++iIndex)
        {
            nodeIndicesInterface[iIndex] = nodeIdsMatrix[iIndex];
            nodeIndicesInterface[nodeIndicesInterface.size() - iIndex - 1] = nodeIdsFibre[iIndex];
        }

        int newElementInterface = ElementCreate(rInterfaceInterpolationType, nodeIndicesInterface);
        GroupAddElement(groupElementsInterface, newElementInterface);

        // create new truss element with duplicated nodes
        int newElementFibre = ElementCreate(rFibreInterpolationType, nodeIdsFibre);
        GroupAddElement(groupElementsFibre, newElementFibre);

        // delete  old element
        ElementDelete(elementId);
    }
    return std::make_pair(groupElementsFibre, groupElementsInterface);
}

void NuTo::Structure::ElementGroupDelete(int rGroupNumber, bool deleteNodes)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, GroupBase>::iterator itGroup = mGroupMap.find(rGroupNumber);
    if (itGroup == mGroupMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType() != NuTo::eGroupId::Elements)
        throw Exception(__PRETTY_FUNCTION__, "Group is not an element group.");

    // the group has to be copied, since the elements are removed from this group, which invalidates the iterators
    Group<ElementBase> copyOfElementGroup = Group<ElementBase>(*(itGroup->second->AsGroupElement()));

    std::set<NodeBase*> potentialNodesToBeRemoved;
    for (auto& elementPtr : copyOfElementGroup)
    {
        // save the nodes, which are eventually to be removed
        if (deleteNodes)
        {
            for (int countNode = 0; countNode < elementPtr.second->GetNumNodes(); countNode++)
            {
                NodeBase* nodePtr = elementPtr.second->GetNode(countNode);
                potentialNodesToBeRemoved.insert(nodePtr);
            }
        }
        ElementDeleteInternal(ElementGetId(elementPtr.second));
    }

    // check all the other elements and see, if they have one of the potential Nodes To Be Removed as valid node
    for (boost::ptr_map<int, ElementBase>::iterator itElement = mElementMap.begin(); itElement != mElementMap.end();
         itElement++)
    {
        for (int countNode = 0; countNode < (itElement->second)->GetNumNodes(); countNode++)
        {
            NodeBase* nodePtr = (itElement->second)->GetNode(countNode);
            // int numRemoved = potentialNodesToBeRemoved.erase(nodePtr);
            potentialNodesToBeRemoved.erase(nodePtr);
        }
    }

    std::map<NodeBase*, int> nodeToId;
    for (auto it = mNodeMap.begin(); it != mNodeMap.end(); ++it)
        nodeToId[it->second] = it->first;

    for (auto node : potentialNodesToBeRemoved)
    {
        int nodeId = nodeToId[node];
        NodeDelete(nodeId, false);
    }
}

//! @brief Deletes an element
//! @param rItElement iterator of the map
void NuTo::Structure::ElementDeleteInternal(int rElementId)
{
    // find element
    boost::ptr_map<int, ElementBase>::iterator itElement = mElementMap.find(rElementId);
    if (itElement == this->mElementMap.end())
    {
        throw Exception(__PRETTY_FUNCTION__, "Element does not exist.");
    }
    else
    {
        // Search for elements in groups: using a loop over all groups
        for (boost::ptr_map<int, GroupBase>::iterator groupIt = mGroupMap.begin(); groupIt != mGroupMap.end();
             ++groupIt)
        {
            if (groupIt->second->GetType() == NuTo::eGroupId::Elements)
            {
                if (groupIt->second->Contain(rElementId))
                {
                    groupIt->second->RemoveMember(rElementId);
                }
            }
        }

        // delete element from map
        this->mElementMap.erase(itElement);
    }
}

// store all elements of a structure in a vector
void NuTo::Structure::GetElementsTotal(std::vector<const ElementBase*>& rElements) const
{
    rElements.reserve(mElementMap.size());
    rElements.resize(0);
    boost::ptr_map<int, ElementBase>::const_iterator ElementIter = this->mElementMap.begin();
    while (ElementIter != this->mElementMap.end())
    {
        rElements.push_back(ElementIter->second);
        ElementIter++;
    }
}

// store all elements of a structure in a vector
void NuTo::Structure::GetElementsTotal(std::vector<std::pair<int, const ElementBase*>>& rElements) const
{
    rElements.reserve(mElementMap.size());
    rElements.resize(0);
    boost::ptr_map<int, ElementBase>::const_iterator ElementIter = this->mElementMap.begin();
    while (ElementIter != this->mElementMap.end())
    {
        rElements.push_back(std::pair<int, const ElementBase*>(ElementIter->first, ElementIter->second));
        ElementIter++;
    }
}

// store all elements of a structure in a vector
void NuTo::Structure::GetElementsTotal(std::vector<ElementBase*>& rElements)
{
    rElements.reserve(mElementMap.size());
    rElements.resize(0);
    boost::ptr_map<int, ElementBase>::iterator ElementIter = this->mElementMap.begin();
    while (ElementIter != this->mElementMap.end())
    {
        rElements.push_back(ElementIter->second);
        ElementIter++;
    }
}

// store all elements of a structure in a vector
void NuTo::Structure::GetElementsTotal(std::vector<std::pair<int, ElementBase*>>& rElements)
{
    rElements.reserve(mElementMap.size());
    rElements.resize(0);
    boost::ptr_map<int, ElementBase>::iterator ElementIter = this->mElementMap.begin();
    while (ElementIter != this->mElementMap.end())
    {
        rElements.push_back(std::pair<int, ElementBase*>(ElementIter->first, ElementIter->second));
        ElementIter++;
    }
}
