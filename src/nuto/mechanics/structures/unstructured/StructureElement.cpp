// $Id$

#include <assert.h>
#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>

#include "nuto/math/FullMatrix.h"

#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/elements/IpDataEnum.h"


#include "nuto/mechanics/elements/ContinuumElement.h"
#include "nuto/mechanics/elements/ContinuumElementIGA.h"
#include "nuto/mechanics/elements/ContinuumElementIGALayer.h"
#include "nuto/mechanics/elements/ContinuumContactElement.h"
#include "nuto/mechanics/elements/ContinuumBoundaryElement.h"
#include "nuto/mechanics/elements/ContinuumBoundaryElementConstrainedControlNode.h"
#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/elements/ElementDataEnum.h"
#include "nuto/mechanics/elements/Element1DInXD.h"
#include "nuto/mechanics/elements/Element2DInterface.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "nuto/mechanics/interpolationtypes/InterpolationBase.h"
#include "nuto/mechanics/interpolationtypes/InterpolationType.h"
#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/nodes/NodeEnum.h"

#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/groups/GroupEnum.h"

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
        throw MechanicsException("[NuTo::Structure::ElementGetElementPtr] Element with identifier " + std::to_string(rIdent) + " does not exist.");
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
        throw MechanicsException("[NuTo::Structure::ElementGetElementPtr] Element with identifier " + std::to_string(rIdent) + " does not exist.");
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
    throw MechanicsException("[NuTo::Structure::GetElementId] Element does not exist.");
}

//! @brief returns a vector with the node ids of an element
//! @param identifier
//! @return vector with node ids
NuTo::FullVector<int, Eigen::Dynamic> NuTo::Structure::ElementGetNodes(int rId)
{
    NuTo::ElementBase* elementPtr = ElementGetElementPtr(rId);
    NuTo::FullVector<int, Eigen::Dynamic> nodeVector(elementPtr->GetNumNodes());
    for (int count = 0; count < elementPtr->GetNumNodes(); count++)
        nodeVector(count) = this->NodeGetId(elementPtr->GetNode(count));
    return nodeVector;
}

//! @brief info about one single element
//! @param rElement (Input) ... pointer to the element
//! @param rVerboseLevel (Input) ... level of verbosity
void NuTo::Structure::ElementInfo(const ElementBase* rElement, int rVerboseLevel) const
{
    std::cout << "element : " << rElement->ElementGetId() << std::endl;
    if (rVerboseLevel > 2)
    {
        std::cout << "\tenum::type=" << Element::ElementTypeToString(rElement->GetEnumType()) << std::endl;
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
                    std::cout << "\t\t" << iIp << ": [" << coor(0,0) << ";" << coor(1,0) << ";" << coor(2,0) << "]" << std::endl;
                }
            }
        }
    }
}

//! @brief info about the elements in the Structure
void NuTo::Structure::ElementInfo(int rVerboseLevel) const
{
    std::cout << "number of elements: " << mElementMap.size() << std::endl;
    if (rVerboseLevel > 3)
    {
        std::cout << "\t\telements :" << std::endl;
        for (boost::ptr_map<int, ElementBase>::const_iterator it = mElementMap.begin(); it != mElementMap.end(); it++)
        {
            std::cout << "\t\t" << it->first;
            if (rVerboseLevel > 4)
            {
                std::cout << "\t:";
                for (unsigned short iNode = 0; iNode < it->second->GetNumNodes(); ++iNode)
                    std::cout << "\t" << this->NodeGetId(it->second->GetNode(iNode));
            }
            std::cout << std::endl;
        }
    }
}

//! @brief changes the node structure to match the interpolation type
//! the node merge distance and the box size are calculated from the element sizes
void NuTo::Structure::ElementTotalConvertToInterpolationType()
{
    // create a group with all elements
#ifdef SHOW_TIME
    bool showTime = mShowTime;
    mShowTime = false;
#endif

    // create element group containing the new elements
    int groupNumber = GroupCreate("Elements");
    GroupBase* group = mGroupMap.find(groupNumber)->second;

    for (boost::ptr_map<int, ElementBase>::iterator itElement = mElementMap.begin(); itElement != mElementMap.end(); ++itElement)
    {
        int elementId = itElement->first;
        ElementBase* elementPtr = itElement->second;
        group->AddMember(elementId, elementPtr);
    }
#ifdef SHOW_TIME
    mShowTime = showTime;
#endif

    // convert elements
    ElementConvertToInterpolationType(groupNumber);

    // delete load from map
    mGroupMap.erase(groupNumber);
}

//! @brief changes the node structure to match the interpolation type
//! @remark The node merge distance is the smallest element length divided by 1000.
//! @remark The mesh size is the median element length divided by 2.
//! @param rGroupNumberElements group for elements (coordinates only) to be converted
void NuTo::Structure::ElementConvertToInterpolationType(int rGroupNumberElements)
{
    // determine the element sizes and calculate the mesh size and the merge distance using these values.

    boost::ptr_map<int, GroupBase>::iterator itGroup = mGroupMap.find(rGroupNumberElements);
    if (itGroup == mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::ElementConvertToInterpolationType] Group with the given identifier does not exist.");

    if (itGroup->second->GetType() != NuTo::eGroupId::Elements)
        throw MechanicsException("[NuTo::Structure::ElementConvertToInterpolationType] Group is not an element group.");

    Group<ElementBase> *elementGroup = itGroup->second->AsGroupElement();
    assert(elementGroup != 0);

    NuTo::FullVector<int, Eigen::Dynamic> elementIndices = elementGroup->GetMemberIds();
    int numElements = elementIndices.GetNumRows();

    // calculate and store the 'size' of each element
    // = volume in 3D
    // = area in 2D
    // = length in 1D
    Eigen::VectorXd elementSize(numElements);

    for (int iElement = 0; iElement < numElements; ++iElement)
    {
        int elementId = elementIndices[iElement];
        ElementBase* element = &(mElementMap.at(elementId));
        Eigen::VectorXd sizeForEachIntegrationPoint = element->GetIntegrationPointVolume();

        elementSize[iElement] = sizeForEachIntegrationPoint.sum();
    }

    std::sort(elementSize.data(), elementSize.data() + numElements, std::less<double>()); // sort in accending order

    double sizeMin = elementSize(0);
    double sizeMed = elementSize(numElements / 2);

    double lengthMin = std::pow(sizeMin, 1. / mDimension);
    double lengthMed = std::pow(sizeMed, 1. / mDimension);

    double mergeDist = lengthMin / 1000.;
    double meshSize = lengthMed / 2.;

    if (mVerboseLevel > 0)
    {
        std::cout << "[NuTo::Structure::ElementConvertToInterpolationType] total size:            " << elementSize.sum() << std::endl;
        std::cout << "[NuTo::Structure::ElementConvertToInterpolationType] min. element size:     " << sizeMin << std::endl;
        std::cout << "[NuTo::Structure::ElementConvertToInterpolationType] median element size:   " << sizeMed << std::endl;
        std::cout << "[NuTo::Structure::ElementConvertToInterpolationType] ====================== " << std::endl;
        std::cout << "[NuTo::Structure::ElementConvertToInterpolationType] min. element length:   " << lengthMin << std::endl;
        std::cout << "[NuTo::Structure::ElementConvertToInterpolationType] median element length: " << lengthMed << std::endl;
        std::cout << "[NuTo::Structure::ElementConvertToInterpolationType] ====================== " << std::endl;
        std::cout << "[NuTo::Structure::ElementConvertToInterpolationType] merge distance:        " << mergeDist << std::endl;
        std::cout << "[NuTo::Structure::ElementConvertToInterpolationType] mesh size:             " << meshSize << std::endl;
    }

    ElementConvertToInterpolationType(rGroupNumberElements, mergeDist, meshSize);

}

//! @brief changes the node structure to match the interpolation type for all elements
//! @param rNodeDistanceMerge Distance of nodes to be joined (should be significantly smaller than the node distance in the mesh)
//! @param rMeshSize approximate size of the elements
void NuTo::Structure::ElementTotalConvertToInterpolationType(double rNodeDistanceMerge, double rMeshSize)
{
    // create a group with all elements

#ifdef SHOW_TIME
    bool showTime = mShowTime;
    mShowTime = false;
#endif

    // create element group containing the new elements
    int groupNumber = GroupCreate("Elements");
    GroupBase* group = mGroupMap.find(groupNumber)->second;

    for (boost::ptr_map<int, ElementBase>::iterator itElement = mElementMap.begin(); itElement != mElementMap.end(); ++itElement)
    {
        int elementId = itElement->first;
        ElementBase* elementPtr = itElement->second;
        group->AddMember(elementId, elementPtr);
    }
#ifdef SHOW_TIME
    mShowTime = showTime;
#endif

    // convert elements
    ElementConvertToInterpolationType(groupNumber, rNodeDistanceMerge, rMeshSize);

    // delete load from map
    mGroupMap.erase(groupNumber);
}

//! @brief changes the node structure to match the interpolation type
//! @param rGroupNumberElements group for elements (coordinates only) to be converted
//! @param rNodeDistanceMerge Distance of nodes to be joined (should be significantly smaller than the node distance in the mesh)
//! @param rMeshSize approximate size of the elements
void NuTo::Structure::ElementConvertToInterpolationType(int rGroupNumberElements, double rNodeDistanceMerge, double rMeshSize)
{
// Summary:
// 1) check stuff, calculate dimensions and number of bounding boxes, suitable for 1D/2D/3D
// 2) Create a TmpNode for each node in each element and store it in the bounding box
// 3) In each bounding box: check for TmpNodes with the same global coordinates
// 4) "Unite" these nodes by creating a new Node with the appropriate dofs
//   4.1) Node already exists: Exchange its pointer in all elements containing it (keep the original node numbering!!!)
//   4.2) new Node: Resize the elements Node storage and add it to all elements containing it
// 5) Profit!

#ifdef SHOW_TIME
    std::clock_t start, end;
    start = clock();
#endif

    // checks

    boost::ptr_map<int, GroupBase>::iterator itGroup = mGroupMap.find(rGroupNumberElements);
    if (itGroup == mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::ElementConvertToInterpolationType] Group with the given identifier does not exist.");

    if (itGroup->second->GetType() != NuTo::eGroupId::Elements)
        throw MechanicsException("[NuTo::Structure::ElementConvertToInterpolationType] Group is not an element group.");

    Group<ElementBase> *elementGroup = itGroup->second->AsGroupElement();
    assert(elementGroup != 0);

    NuTo::FullVector<int, Eigen::Dynamic> elementIndices = elementGroup->GetMemberIds();
    int numElements = elementIndices.GetNumRows();

//    const InterpolationType* interpolationType = this->ElementGetElementPtr(elementIndices.at(0,0))->GetInterpolationType();

    std::vector<ElementBase*> elements(numElements);
    for (int iElement = 0; iElement < numElements; ++iElement)
    {
        int elementId = elementIndices[iElement];
        elements[iElement] = &(mElementMap.at(elementId));
//        if (interpolationType != elements[iElement]->GetInterpolationType())
//            throw MechanicsException("[NuTo::Structure::ElementConvertToInterpolationType] All elements must have the same interpolation type.");
    }




    // getting rid of the two biggest bottlenecks
    // 1) always looping trough all elements when exchanging the node pointer
    //    instead of just the ones containing the node
    //    --> build a map: key-node; value list of elements containing the key;
    // 2) finding the node id
    //    --> build a reversed node map, NodeBase* --> NodeId
    std::map<NodeBase*, std::vector<ElementBase*>> nodeToElement;
    std::map<const NodeBase*, int> nodeToId;

    for (boost::ptr_map<int,NodeBase>::const_iterator it = mNodeMap.begin(); it!= mNodeMap.end(); it++)
    {
        nodeToId[it->second] = it->first;
    }

    // find the bounding box corners = max/min x,y,z coordinates, initialize with random node coordinates
    Eigen::VectorXd boundingBoxMax = elements[0]->GetNode(0)->Get(Node::eDof::COORDINATES);
    Eigen::VectorXd boundingBoxMin = elements[0]->GetNode(0)->Get(Node::eDof::COORDINATES);

    for (ElementBase* element : elements)
    {
        // loop through nodes with coordinates
        for (int iNode = 0; iNode < element->GetNumNodes(Node::eDof::COORDINATES); ++iNode)
        {
            NodeBase* node = element->GetNode(iNode, Node::eDof::COORDINATES);
            Eigen::VectorXd nodeCoordinates = node->Get(Node::eDof::COORDINATES);
            assert(nodeCoordinates.rows() == mDimension);
            boundingBoxMax = boundingBoxMax.cwiseMax(nodeCoordinates);
            boundingBoxMin = boundingBoxMin.cwiseMin(nodeCoordinates);

            nodeToElement[node].push_back(element);
        }
    }


    Eigen::Vector3i numBoxes = Eigen::Vector3i::Ones(); // initialize as 3D array with at least one box per dimension
    Eigen::VectorXd deltaBox(mDimension);
    for (int iDim = 0; iDim < mDimension; ++iDim)
    {
        // sligtly move the boundingBoxMin
        boundingBoxMin[iDim] -= rNodeDistanceMerge;
        numBoxes[iDim] = (boundingBoxMax[iDim] - boundingBoxMin[iDim]) / rMeshSize + 1;
        deltaBox[iDim] = (boundingBoxMax[iDim] - boundingBoxMin[iDim] + rNodeDistanceMerge) / numBoxes[iDim]; // small offset to fit coordinates at boundingBoxMax in the last box
        if (deltaBox[iDim] < rNodeDistanceMerge)
            throw MechanicsException("[NuTo::Structure::ElementConvertToInterpolationType] The merge distance is larger than the mesh size, that should not happen.");
    }

    if (mVerboseLevel > 0)
    {
        mLogger << "[NuTo::Structure::ElementConvertToInterpolationType] bounding box: from  " << boundingBoxMin.transpose() << " to " << boundingBoxMax.transpose() << "\n";
        mLogger << "[NuTo::Structure::ElementConvertToInterpolationType] number of boxes: " << numBoxes.transpose() << "\n";
        mLogger << "[NuTo::Structure::ElementConvertToInterpolationType] distance  of boxes: " << deltaBox.transpose() << "\n";
    }

    if ((boundingBoxMax - boundingBoxMin).minCoeff() <= 0)
        // The group has at least one element and this element has no spacial size.
        throw MechanicsException("[NuTo::Structure::ElementConvertToInterpolationType] Bounding box with zero length. Your element definition might be messed up.");

    double nodeDistanceMerge2 = rNodeDistanceMerge * rNodeDistanceMerge;

    // data structure for storing temporary node informations
    struct TmpNode
    {
        TmpNode(const Eigen::VectorXd& rCoords, ElementBase* rElement, int rElementNodeId, double rDist) :
                coords(rCoords), elementNodeId(rElementNodeId), mDist(rDist)
        {
            element = rElement;
        }

        const Eigen::VectorXd coords;   // global coordinates, either 1D, 2D or 3D
        ElementBase* element;           // single element connected to this node
        const int elementNodeId;        // local node id within that element
        const double mDist;             // merge distance squared

        bool operator==(const TmpNode& rOther) const
        {
            return (coords - rOther.coords).dot(coords - rOther.coords) < mDist;
        }
    };

    // Use 1D storage for the bounding boxes. Each bounding box contains a vector of TmpNodes
    std::vector<std::vector<TmpNode>> boxes(numBoxes.prod());

    for (ElementBase* element : elements)
    {
        const InterpolationType* interpolationType = element->GetInterpolationType();

        for (int iNode = 0; iNode < element->GetNumNodes(); ++iNode)
        {

            Eigen::VectorXd naturalNodeCoordinates = interpolationType->GetNaturalNodeCoordinates(iNode);
            Eigen::VectorXd globalNodeCoordinates = element->InterpolateDofGlobal(naturalNodeCoordinates, Node::eDof::COORDINATES);

            Eigen::Vector3i boxIndex3D = Eigen::Vector3i::Zero(); // default access is the first box
            for (int iDim = 0; iDim < mDimension; ++iDim)
                boxIndex3D[iDim] = std::floor((globalNodeCoordinates(iDim,0) - boundingBoxMin(iDim,0)) / deltaBox(iDim,0));

            // 3D-index --> 1D-index mapping
            int boxIndex1D = boxIndex3D[0] + boxIndex3D[1] * numBoxes[0] + boxIndex3D[2] * numBoxes[0] * numBoxes[1];
            assert(boxIndex1D < (int )boxes.size());

            boxes[boxIndex1D].push_back(TmpNode(globalNodeCoordinates, element, iNode, nodeDistanceMerge2));
        }
    }

    int numEmptyBoxes = 0; // for statistics only
    int maxBoxSize = -1;   // for statistics only
    int sumBoxSize = 0;    // for statistics only

    for (std::vector<TmpNode>& box : boxes)
    {
        int boxSize = box.size();
        maxBoxSize = std::max(maxBoxSize, boxSize);
        sumBoxSize += boxSize;
        if (boxSize == 0)
            numEmptyBoxes++;

        std::vector<std::vector<TmpNode>> sameNodes; // multiple TmpNodes with the same global coordinates == the same nodes
        for (TmpNode& node : box)
        {
            unsigned int sameNodeIndex = sameNodes.size();
            for (unsigned int iUniqueNode = 0; iUniqueNode < sameNodes.size(); ++iUniqueNode)
                if (sameNodes[iUniqueNode][0] == node)
                {
                    sameNodeIndex = iUniqueNode;
                    break;
                }

            if (sameNodeIndex == sameNodes.size())
            {
                // no match with existing nodes --> create new vector
                sameNodes.push_back(std::vector<TmpNode>());
            }
            sameNodes[sameNodeIndex].push_back(node);
        }

        for (std::vector<TmpNode>& singleSameNode : sameNodes)
        {
            TmpNode& tmpNode = singleSameNode[0]; // the first one always exists

            // check if all same nodes (= nodes with the same global coordinates) have the same dofs
            std::set<Node::eDof> nodeDofs;
            for (TmpNode& singleUniqueNode : singleSameNode)
            {
                auto nodeDofsOther = singleUniqueNode.element->GetInterpolationType()->GetNodeDofs(singleUniqueNode.elementNodeId);
                for (auto dof : nodeDofsOther)
                    nodeDofs.insert(dof);
            }

            NuTo::FullVector<double, Eigen::Dynamic> nodeCoordinates(tmpNode.coords);

            if (nodeDofs.find(Node::eDof::COORDINATES) != nodeDofs.end())
            {

                // If the node is a coordinate node, it should already be in the elements mNodes, check that:
                NodeBase* oldNode = singleSameNode[0].element->GetNode(singleSameNode[0].elementNodeId);
                assert(oldNode != nullptr);
                for (TmpNode& singleUniqueNode : singleSameNode)
                    assert(singleUniqueNode.element->GetNode(singleUniqueNode.elementNodeId) == oldNode);

                NodeBase* newNode = NodePtrCreate(nodeDofs, nodeCoordinates);

                // exchange the pointer in all elements, all groups, all constraints etc
                int oldNodeIndex = nodeToId[oldNode];
                std::vector<ElementBase*> elementsToChange = nodeToElement[oldNode];
                NodeExchangePtr(oldNodeIndex, oldNode, newNode, elementsToChange); // from now on, newNode is nullPtr

            } else
            {
                // create node manually
                NodeBase* newNode = NodePtrCreate(nodeDofs, nodeCoordinates);

                int id(mNodeMap.size());
                boost::ptr_map<int, NodeBase>::iterator it = mNodeMap.find(id);
                while (it != mNodeMap.end())
                {
                    id++;
                    it = mNodeMap.find(id);
                }
                mNodeMap.insert(id, newNode);

                // set the node pointer for all suitable elements
                for (TmpNode& singleUniqueNode : singleSameNode)
                {
                    ElementBase* theElement = singleUniqueNode.element;
                    int theNodeIndex = singleUniqueNode.elementNodeId;
                    theElement->ResizeNodes(theElement->GetInterpolationType()->GetNumNodes());

                    theElement->SetNode(theNodeIndex, newNode);
                }
            }

        }
    }

    if (mVerboseLevel > 0)
    {
        mLogger << "[NuTo::Structure::ElementConvertToInterpolationType] Maximum box entries: " << maxBoxSize << "\n";
        mLogger << "[NuTo::Structure::ElementConvertToInterpolationType] number of boxes:     " << numBoxes.prod() << "\n";
        mLogger << "[NuTo::Structure::ElementConvertToInterpolationType] Average box entries: " << static_cast<double>(sumBoxSize) / numBoxes.prod() << "\n";
        mLogger << "[NuTo::Structure::ElementConvertToInterpolationType] Empty boxes:         " << static_cast<double>(numEmptyBoxes) / numBoxes.prod() * 100 << "%" << "\n";
    }

#ifdef SHOW_TIME
    end = clock();
    if (mShowTime)
        std::cout << "[NuTo::Structure::ElementConvertToInterpolationType] " << difftime(end, start) / CLOCKS_PER_SEC << "sec" << std::endl;
#endif

}

//! @brief Deletes an element
//! @param rElementIdent identifier for the element
void NuTo::Structure::ElementDelete(int rElementNumber)
{
#ifdef SHOW_TIME
    std::clock_t start, end;
    start = clock();
#endif
    ElementDeleteInternal(rElementNumber);
#ifdef SHOW_TIME
    end = clock();
    if (mShowTime)
        std::cout << "[NuTo::Structure::ElementDelete] " << difftime(end, start) / CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

int NuTo::Structure::ElementCreate(int rInterpolationTypeId,
                                   const Eigen::VectorXi &rNodeNumbers,
                                   const Eigen::MatrixXd &rKnots,
                                   const Eigen::VectorXi &rKnotIDs)
{
    return ElementCreate(rInterpolationTypeId, rNodeNumbers, rKnots, rKnotIDs,  "CONSTITUTIVELAWIP", "NOIPDATA");
}


int NuTo::Structure::ElementCreate(int rInterpolationTypeId, const NuTo::FullVector<int, Eigen::Dynamic>& rNodeNumbers)
{
    return ElementCreate(rInterpolationTypeId, rNodeNumbers, "CONSTITUTIVELAWIP", "NOIPDATA");
}

int NuTo::Structure::ElementCreate(int   rInterpolationTypeId,
                                   const Eigen::VectorXi &rNodeNumbers,
                                   const Eigen::MatrixXd &rKnots,
                                   const Eigen::VectorXi &rKnotIDs,
                                   const std::string     &rElementDataType,
                                   const std::string     &rIpDataType)
{
    //find unused integer id
    int elementNumber(mElementMap.size());
    boost::ptr_map<int, ElementBase>::iterator it = mElementMap.find(elementNumber);
    while (it != mElementMap.end())
    {
        elementNumber++;
        it = mElementMap.find(elementNumber);
    }
    // create element
    this->ElementCreate(elementNumber, rInterpolationTypeId, rNodeNumbers, rKnots, rKnotIDs, rElementDataType, rIpDataType);
    return elementNumber;
}

int NuTo::Structure::ElementCreate(int   rInterpolationTypeId,
                                   const NuTo::FullVector<int, Eigen::Dynamic>& rNodeNumbers,
                                   const std::string& rElementDataType,
                                   const std::string& rIpDataType)
{

    //find unused integer id
    int elementNumber(mElementMap.size());
    boost::ptr_map<int, ElementBase>::iterator it = mElementMap.find(elementNumber);
    while (it != mElementMap.end())
    {
        elementNumber++;
        it = mElementMap.find(elementNumber);
    }
    // create element
    this->ElementCreate(elementNumber, rInterpolationTypeId, rNodeNumbers, rElementDataType, rIpDataType);
    return elementNumber;
}

void NuTo::Structure::ElementCreate(int rElementNumber,
                                    int rInterpolationTypeId,
                                    const Eigen::VectorXi& rNodeNumbers,
                                    const Eigen::MatrixXd &rKnots,
                                    const Eigen::VectorXi &rKnotIDs,
                                    const std::string& rElementDataType,
                                    const std::string& rIpDataType)
{
    // convert node numbers to pointer
    if (rNodeNumbers.cols() != 1)
        throw MechanicsException("[NuTo::Structure::ElementCreate] Matrix with node numbers should have a single column.");
    std::vector<NodeBase*> nodeVector;
    for (int iNode = 0; iNode < rNodeNumbers.rows(); iNode++)
    {
        nodeVector.push_back(NodeGetNodePtr(rNodeNumbers(iNode)));
    }

    ElementData::eElementDataType elementDataType = ElementData::ElementDataTypeToEnum(rElementDataType);
    IpData::eIpDataType ipDataType = NuTo::IpData::IpDataTypeToEnum(rIpDataType);

    ElementCreate(rElementNumber, rInterpolationTypeId, nodeVector, rKnots, rKnotIDs, elementDataType, ipDataType);
}

void NuTo::Structure::ElementCreate(int rElementNumber,
                                    int rInterpolationTypeId,
                                    const NuTo::FullVector<int, Eigen::Dynamic>& rNodeNumbers,
                                    const std::string& rElementDataType,
                                    const std::string& rIpDataType)
{
    // convert node numbers to pointer
    if (rNodeNumbers.GetNumColumns() != 1)
        throw MechanicsException("[NuTo::Structure::ElementCreate] Matrix with node numbers should have a single column.");
    std::vector<NodeBase*> nodeVector;
    for (int iNode = 0; iNode < rNodeNumbers.GetNumRows(); iNode++)
    {
        nodeVector.push_back(NodeGetNodePtr(rNodeNumbers.GetValue(iNode)));
    }

    ElementData::eElementDataType elementDataType = ElementData::ElementDataTypeToEnum(rElementDataType);
    IpData::eIpDataType ipDataType = NuTo::IpData::IpDataTypeToEnum(rIpDataType);

    ElementCreate(rElementNumber, rInterpolationTypeId, nodeVector, elementDataType, ipDataType);
}

int NuTo::Structure::ElementCreate(int rInterpolationTypeId,
                                   const std::vector<NodeBase*>& rNodeNumbers,
                                   const Eigen::MatrixXd &rKnots,
                                   const Eigen::VectorXi &rKnotIDs,
                                   ElementData::eElementDataType rElementDataType,
                                   IpData::eIpDataType rIpDataType)
{
    //find unused integer id
    int elementNumber(mElementMap.size());
    boost::ptr_map<int, ElementBase>::iterator it = mElementMap.find(elementNumber);
    while (it != mElementMap.end())
    {
        elementNumber++;
        it = mElementMap.find(elementNumber);
    }

    // create element
    this->ElementCreate(elementNumber, rInterpolationTypeId, rNodeNumbers, rKnots, rKnotIDs, rElementDataType, rIpDataType);

    // return element number
    return elementNumber;
}

void  NuTo::Structure::ElementCreate(int rElementNumber, int rInterpolationTypeId, const std::vector<int>& rNodeIds, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType)
{
    std::vector<NodeBase*> nodeVector;
    for (const auto& nodeId : rNodeIds)
        nodeVector.push_back(NodeGetNodePtr(nodeId));

    ElementCreate(rElementNumber, rInterpolationTypeId, nodeVector, rElementDataType, rIpDataType);
}

//! @brief Creates an element
//! @param rInterpolationTypeId interpolation type id
//! @param rNodeVector pointers to the corresponding nodes
//! @param rElementDataType Element data for the elements
//! @param rIpDataType Integration point data for the elements
//! @return int rElementNumber
int NuTo::Structure::ElementCreate(int rInterpolationTypeId, const std::vector<NodeBase*>& rNodeNumbers, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType)
{
    //find unused integer id
    int elementNumber(mElementMap.size());
    boost::ptr_map<int, ElementBase>::iterator it = mElementMap.find(elementNumber);
    while (it != mElementMap.end())
    {
        elementNumber++;
        it = mElementMap.find(elementNumber);
    }

    // create element
    this->ElementCreate(elementNumber, rInterpolationTypeId, rNodeNumbers, rElementDataType, rIpDataType);

    // return element number
    return elementNumber;
}

int NuTo::Structure::ElementCreate(int rInterpolationTypeId,
                                   const Eigen::VectorXi& rNodeNumbers,
                                   const Eigen::MatrixXd &rKnots,
                                   const Eigen::VectorXi &rKnotIDs,
                                   ElementData::eElementDataType rElementDataType,
                                   IpData::eIpDataType rIpDataType)
{
    //find unused integer id
    int elementNumber(mElementMap.size());
    boost::ptr_map<int, ElementBase>::iterator it = mElementMap.find(elementNumber);
    while (it != mElementMap.end())
    {
        elementNumber++;
        it = mElementMap.find(elementNumber);
    }

    // convert node numbers to pointer
    if (rNodeNumbers.cols() != 1)
        throw MechanicsException("[NuTo::Structure::ElementCreate] Matrix with node numbers should have a single column.");
    std::vector<NodeBase*> nodeVector;
    for (int iNode = 0; iNode < rNodeNumbers.rows(); iNode++)
    {
        nodeVector.push_back(NodeGetNodePtr(rNodeNumbers(iNode)));
    }

    // create element
    this->ElementCreate(elementNumber, rInterpolationTypeId, nodeVector, rKnots, rKnotIDs, rElementDataType, rIpDataType);

    // return element number
    return elementNumber;
}

int NuTo::Structure::ElementCreate(int rInterpolationTypeId,
                                   const NuTo::FullVector<int, Eigen::Dynamic>& rNodeNumbers,
                                   ElementData::eElementDataType rElementDataType,
                                   IpData::eIpDataType rIpDataType)
{
    //find unused integer id
    int elementNumber(mElementMap.size());
    boost::ptr_map<int, ElementBase>::iterator it = mElementMap.find(elementNumber);
    while (it != mElementMap.end())
    {
        elementNumber++;
        it = mElementMap.find(elementNumber);
    }

    // convert node numbers to pointer
    if (rNodeNumbers.GetNumColumns() != 1)
        throw MechanicsException("[NuTo::Structure::ElementCreate] Matrix with node numbers should have a single column.");
    std::vector<NodeBase*> nodeVector;
    for (int iNode = 0; iNode < rNodeNumbers.GetNumRows(); iNode++)
    {
        nodeVector.push_back(NodeGetNodePtr(rNodeNumbers.GetValue(iNode)));
    }

    // create element
    this->ElementCreate(elementNumber, rInterpolationTypeId, nodeVector, rElementDataType, rIpDataType);

    // return element number
    return elementNumber;
}


void NuTo::Structure::ElementCreate(int rElementNumber,
                                    int rInterpolationTypeId,
                                    const std::vector<NodeBase*>& rNodeVector,
                                    const Eigen::MatrixXd &rKnots,
                                    const Eigen::VectorXi &rKnotIDs,
                                    ElementData::eElementDataType rElementDataType,
                                    IpData::eIpDataType rIpDataType)
{
    boost::ptr_map<int, InterpolationType>::iterator itIterator = mInterpolationTypeMap.find(rInterpolationTypeId);

    if (itIterator == mInterpolationTypeMap.end())
        throw NuTo::MechanicsException("[NuTo::Structure::ElementCreate] Interpolation type does not exist.");

    InterpolationType* interpolationType = itIterator->second;

    if (interpolationType->IsDof(Node::eDof::COORDINATES) == false)
        throw NuTo::MechanicsException("[NuTo::Structure::ElementCreate] COORDINATE interpolation required.");

    unsigned int numNodesCoordinates = interpolationType->Get(Node::eDof::COORDINATES).GetNumNodes();
    if (numNodesCoordinates != rNodeVector.size())
        throw NuTo::MechanicsException("[NuTo::Structure::ElementCreate] COORDINATE interpolation requires " + std::to_string(numNodesCoordinates) + " nodes. " + std::to_string(rNodeVector.size()) + " are provided.");

    ElementBase* ptrElement = 0;

    if (interpolationType->GetCurrentIntegrationType() == nullptr)
    {
        eIntegrationType integrationTypeEnum = interpolationType->GetStandardIntegrationType();
        IntegrationTypeBase* integrationType = GetPtrIntegrationType(integrationTypeEnum);
        interpolationType->UpdateIntegrationType(*integrationType);
    }

    try
    {
        switch (interpolationType->GetShapeType())
        {
        case NuTo::Interpolation::eShapeType::SPRING:
            throw NuTo::MechanicsException("[NuTo::Structure::ElementCreate] Element1DSpring currently not implemented.");
//            ptrElement = new Element1DSpring(this, rNodeVector, rElementDataType, rIpDataType, interpolationType);
//            ptrElement->CheckElement();
            break;
        case NuTo::Interpolation::eShapeType::TRUSS1D:
        case NuTo::Interpolation::eShapeType::TRUSSXD:
        case NuTo::Interpolation::eShapeType::TRIANGLE2D:
        case NuTo::Interpolation::eShapeType::QUAD2D:
        case NuTo::Interpolation::eShapeType::TETRAHEDRON3D:
        case NuTo::Interpolation::eShapeType::BRICK3D:
        case NuTo::Interpolation::eShapeType::INTERFACE:
            throw NuTo::MechanicsException("[NuTo::Structure::ElementCreate] Please use approriate functions for element creation, this is IGA implementation.");
            break;
        case NuTo::Interpolation::eShapeType::IGA1D:
            ptrElement = new ContinuumElementIGA<1>(this, rNodeVector, rKnots, rKnotIDs, rElementDataType, rIpDataType, interpolationType);
            ptrElement->CheckElement();
            break;
        case NuTo::Interpolation::eShapeType::IGA2D:
            ptrElement = new ContinuumElementIGA<2>(this, rNodeVector, rKnots, rKnotIDs, rElementDataType, rIpDataType, interpolationType);
            ptrElement->CheckElement();
            break;
        case NuTo::Interpolation::eShapeType::IGA1DLAYER:
            ptrElement = new ContinuumElementIGALayer<1>(this, rNodeVector, rKnots, rKnotIDs, rElementDataType, rIpDataType, interpolationType);
            ptrElement->CheckElement();
            break;
        case NuTo::Interpolation::eShapeType::IGA2DLAYER:
            ptrElement = new ContinuumElementIGALayer<2>(this, rNodeVector, rKnots, rKnotIDs, rElementDataType, rIpDataType, interpolationType);
            ptrElement->CheckElement();
            break;
        default:
            throw MechanicsException("[NuTo::Structure::ElementCreate] invalid dimension.");
        }

        mElementMap.insert(rElementNumber, ptrElement);

    } catch (NuTo::MechanicsException &e)
    {
        std::stringstream ss;
        ss << rElementNumber;
        e.AddMessage("[NuTo::Structure::ElementCreate] Error creating element " + ss.str() + ".");
        throw e;
    } catch (...)
    {
        std::stringstream ss;
        ss << rElementNumber;
        throw NuTo::MechanicsException("[NuTo::Structure::ElementCreate] Error creating element " + ss.str() + ".");
    }
}

//! @brief Creates an element
//! @param rElementNumber element number
//! @param rInterpolationTypeId interpolation type id
//! @param rNodeVector pointers to the corresponding nodes
//! @param rElementType element type
//! @param rIpDataType Integration point data for the elements
void NuTo::Structure::ElementCreate(int rElementNumber, int rInterpolationTypeId, const std::vector<NodeBase*>& rNodeVector, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType)
{
    boost::ptr_map<int, InterpolationType>::iterator itIterator = mInterpolationTypeMap.find(rInterpolationTypeId);

    if (itIterator == mInterpolationTypeMap.end())
        throw NuTo::MechanicsException("[NuTo::Structure::ElementCreate] Interpolation type does not exist.");

    InterpolationType* interpolationType = itIterator->second;

    if (interpolationType->IsDof(Node::eDof::COORDINATES) == false)
        throw NuTo::MechanicsException("[NuTo::Structure::ElementCreate] COORDINATE interpolation required.");

    unsigned int numNodesCoordinates = interpolationType->Get(Node::eDof::COORDINATES).GetNumNodes();
    if (numNodesCoordinates != rNodeVector.size())
        throw NuTo::MechanicsException("[NuTo::Structure::ElementCreate] COORDINATE interpolation requires " + std::to_string(numNodesCoordinates) + " nodes. " + std::to_string(rNodeVector.size()) + " are provided.");

    ElementBase* ptrElement = 0;

    if (interpolationType->GetCurrentIntegrationType() == nullptr)
    {
        eIntegrationType integrationTypeEnum = interpolationType->GetStandardIntegrationType();
        IntegrationTypeBase* integrationType = GetPtrIntegrationType(integrationTypeEnum);
        interpolationType->UpdateIntegrationType(*integrationType);
    }

    try
    {
        switch (interpolationType->GetShapeType())
        {
        case NuTo::Interpolation::eShapeType::SPRING:
            throw NuTo::MechanicsException("[NuTo::Structure::ElementCreate] Element1DSpring currently not implemented.");
//            ptrElement = new Element1DSpring(this, rNodeVector, rElementDataType, rIpDataType, interpolationType);
//            ptrElement->CheckElement();
            break;
        case NuTo::Interpolation::eShapeType::TRUSS1D:
            ptrElement = new ContinuumElement<1>(this, rNodeVector, rElementDataType, rIpDataType, interpolationType);
            ptrElement->CheckElement();
            break;
        case NuTo::Interpolation::eShapeType::TRUSSXD:
            ptrElement = new Element1DInXD(this, rNodeVector, rElementDataType, rIpDataType, interpolationType);
            ptrElement->CheckElement();
            break;
        case NuTo::Interpolation::eShapeType::TRIANGLE2D:
        case NuTo::Interpolation::eShapeType::QUAD2D:
            ptrElement = new ContinuumElement<2>(this, rNodeVector, rElementDataType, rIpDataType, interpolationType);
            ptrElement->CheckElement();
            break;
        case NuTo::Interpolation::eShapeType::TETRAHEDRON3D:
        case NuTo::Interpolation::eShapeType::BRICK3D:
            ptrElement = new ContinuumElement<3>(this, rNodeVector, rElementDataType, rIpDataType, interpolationType);
            ptrElement->CheckElement();
            break;
        case NuTo::Interpolation::eShapeType::INTERFACE:
            ptrElement = new Element2DInterface(this, rNodeVector, rElementDataType, rIpDataType, interpolationType);
            ptrElement->CheckElement();
            break;
        case NuTo::Interpolation::eShapeType::IGA1D:
            throw NuTo::MechanicsException("[NuTo::Structure::ElementCreate] Please use the ElementCreate function for IGA elements, where the knot parameters are provided.");
            break;
        case NuTo::Interpolation::eShapeType::IGA2D:
            throw NuTo::MechanicsException("[NuTo::Structure::ElementCreate] Please use the ElementCreate function for IGA elements, where the knot parameters are provided.");
            break;
        default:
            throw MechanicsException("[NuTo::Structure::ElementCreate] invalid dimension.");
        }

        mElementMap.insert(rElementNumber, ptrElement);

    } catch (NuTo::MechanicsException &e)
    {
        std::stringstream ss;
        ss << rElementNumber;
        e.AddMessage("[NuTo::Structure::ElementCreate] Error creating element " + ss.str() + ".");
        throw e;
    } catch (...)
    {
        std::stringstream ss;
        ss << rElementNumber;
        throw NuTo::MechanicsException("[NuTo::Structure::ElementCreate] Error creating element " + ss.str() + ".");
    }
}

//! @brief creates multiple elements and adds them to an element group
//! @param rInterpolationTypeId interpolation type id
//! @param rNodeIdents Identifier for the corresponding nodes (Incidences have to be stored column-wise)
//! @return index to the new element group
int NuTo::Structure::ElementsCreate(int rInterpolationTypeId, NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic>& rNodeNumbers)
{
    return ElementsCreate(rInterpolationTypeId, rNodeNumbers, "CONSTITUTIVELAWIP", "NOIPDATA");
}

//! @brief creates multiple elements and adds them to an element group
//! @param rInterpolationTypeId interpolation type id
//! @param rNodeIdents Identifier for the corresponding nodes (Incidences have to be stored column-wise)
//! @param rElementDataType Element data for the elements
//! @param rIpDataType Integration point data for the elements
//! @return index to the new element group
int NuTo::Structure::ElementsCreate(int rInterpolationTypeId, NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic>& rNodeNumbers, const std::string& rElementDataType, const std::string& rIpDataType)
{
    std::vector<int> newElementIds;
    /// go through the elements
    for (int iNode = 0; iNode < rNodeNumbers.GetNumColumns(); ++iNode)
    {
        const NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> incidence(rNodeNumbers.GetColumn(iNode));
        int newElementId = ElementCreate(rInterpolationTypeId, incidence, rElementDataType, rIpDataType);
        newElementIds.push_back(newElementId);
    }

#ifdef SHOW_TIME
    bool showTime = mShowTime;
    mShowTime = false;
#endif

    // create element group containing the new elements
    int newElementGroup = GroupCreate("Elements");
    for (int newElementId : newElementIds)
        GroupAddElement(newElementGroup, newElementId);

#ifdef SHOW_TIME
    mShowTime = showTime;
#endif
    return newElementGroup;
}

namespace NuTo
{

template<>
Eigen::Matrix<std::pair<const ContinuumElementIGA<1>*, int>, Eigen::Dynamic, Eigen::Dynamic> NuTo::Structure::ContactElementsCreateMaster<1>(const Eigen::Matrix<std::pair<int, int>, Eigen::Dynamic, Eigen::Dynamic> &rMasterElementsID)
{
    Eigen::Matrix<std::pair<const ContinuumElementIGA<1>*, int>, Eigen::Dynamic, Eigen::Dynamic> masterElements(rMasterElementsID.rows(), rMasterElementsID.cols());

    for(int i = 0 ; i < rMasterElementsID.rows(); i++)
    {
        for(int j = 0 ; j < rMasterElementsID.cols(); j++)
        {
            const NuTo::ElementBase* elementPtrSlave = ElementGetElementPtr(rMasterElementsID(i,j).first);

            if(elementPtrSlave->GetEnumType() == Element::eElementType::CONTINUUMELEMENTIGA)
            {
                if(elementPtrSlave->GetLocalDimension() == 1)
                {
                    masterElements(i,j) = std::make_pair(&elementPtrSlave->AsContinuumElementIGA1D(), rMasterElementsID(i,j).second);
                }
                else
                    throw MechanicsException(__PRETTY_FUNCTION__,"Wrong element dimension.");
            }
            else
                throw MechanicsException(__PRETTY_FUNCTION__,"Only ContinuumElementIGA available as master elements.");
        }
    }
    return masterElements;
}

template<>
Eigen::Matrix<std::pair<const ContinuumElementIGA<2>*, int>, Eigen::Dynamic, Eigen::Dynamic> NuTo::Structure::ContactElementsCreateMaster<2>(const Eigen::Matrix<std::pair<int, int>, Eigen::Dynamic, Eigen::Dynamic> &rMasterElementsID)
{
    Eigen::Matrix<std::pair<const ContinuumElementIGA<2>*, int>, Eigen::Dynamic, Eigen::Dynamic> masterElements(rMasterElementsID.rows(), rMasterElementsID.cols());

    for(int i = 0 ; i < rMasterElementsID.rows(); i++)
    {
        for(int j = 0 ; j < rMasterElementsID.cols(); j++)
        {
            const NuTo::ElementBase* elementPtrSlave = ElementGetElementPtr(rMasterElementsID(i,j).first);

            if(elementPtrSlave->GetEnumType() == Element::eElementType::CONTINUUMELEMENTIGA)
            {
                if(elementPtrSlave->GetLocalDimension() == 2)
                {
                    masterElements(i,j) = std::make_pair(&elementPtrSlave->AsContinuumElementIGA2D(), rMasterElementsID(i,j).second);
                }
                else
                    throw MechanicsException(__PRETTY_FUNCTION__,"Wrong element dimension.");
            }
            else
                throw MechanicsException(__PRETTY_FUNCTION__,"Only ContinuumElementIGA available as master elements.");
        }
    }
    return masterElements;
}
}
template<int TDimSlave, int TDimMaster>
void NuTo::Structure::ContactElementsCreate(int rElementsGroupIDSlave,
                                            int rNodeGroupSlaveId,
                                            const Eigen::Matrix<std::pair<int, int>,
                                            Eigen::Dynamic, Eigen::Dynamic> &rMasterElementsID,
                                            eIntegrationType rIntegrationType,
                                            double rPenalty,
                                            int rContactAlgorithm)
{
    Eigen::Matrix<std::pair<const ContinuumElementIGA<TDimMaster>*, int>, Eigen::Dynamic, Eigen::Dynamic> masterElements = ContactElementsCreateMaster<TDimMaster>(rMasterElementsID);

    //find groups
    boost::ptr_map<int, GroupBase>::iterator itGroupElementsSlave = mGroupMap.find(rElementsGroupIDSlave);
    if (itGroupElementsSlave == mGroupMap.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroupElementsSlave->second->GetType() != NuTo::eGroupId::Elements)
        throw MechanicsException(__PRETTY_FUNCTION__,"Group is not an element group.");

    Group<ElementBase>& elementGroupSlave = *(itGroupElementsSlave->second->AsGroupElement());


    boost::ptr_map<int, GroupBase>::iterator itGroupBoundaryNodes = mGroupMap.find(rNodeGroupSlaveId);
    if (itGroupBoundaryNodes == mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::BoundaryElementsCreate] Group with the given identifier does not exist.");
    if (itGroupBoundaryNodes->second->GetType() != NuTo::eGroupId::Nodes)
        throw MechanicsException("[NuTo::Structure::BoundaryElementsCreate] Group is not a node group.");

    Group<NodeBase>&    nodeGroup       = *(itGroupBoundaryNodes->second->AsGroupNode());

    // since the search is done via the id's, the surface nodes are ptr, so make another set with the node ptrs
    std::set<const NodeBase*> nodePtrSet;
    for (auto itNode : nodeGroup)
    {
        nodePtrSet.insert(itNode.second);
    }

    std::vector<int> newBoundaryElementIds;

    //loop over all elements
    for (auto itElementSlave : elementGroupSlave)
    {
        ElementBase* elementPtrSlave = itElementSlave.second;
         const InterpolationType* interpolationType = elementPtrSlave->GetInterpolationType();
        //loop over all surfaces
        for (int iSurface = 0; iSurface < interpolationType->GetNumSurfaces(); ++iSurface)
        {
            bool elementSurfaceNodesMatchBoundaryNodes = true;
            Eigen::VectorXi surfaceNodeIndices = interpolationType->GetSurfaceNodeIndices(iSurface);

            int numSurfaceNodes = surfaceNodeIndices.rows();
            std::vector<const NodeBase*> surfaceNodes(numSurfaceNodes);

            for (int iSurfaceNode = 0; iSurfaceNode < numSurfaceNodes; ++iSurfaceNode)
            {
                surfaceNodes[iSurfaceNode] = elementPtrSlave->GetNode(surfaceNodeIndices(iSurfaceNode, 0));
            }

            //check, if all surface nodes are in the node group
            for (unsigned int countNode = 0; countNode < surfaceNodes.size(); countNode++)
            {
                if (nodePtrSet.find(surfaceNodes[countNode]) == nodePtrSet.end())
                {
                    //this surface has at least one node that is not in the list, continue
                    elementSurfaceNodesMatchBoundaryNodes = false;
                    break;
                }
            }

            if (elementSurfaceNodesMatchBoundaryNodes)
            {
                int surfaceId = iSurface;

                IpData::eIpDataType ipDataType = elementPtrSlave->GetIpDataType(0);
                ConstitutiveBase* constitutiveLaw = elementPtrSlave->GetConstitutiveLaw(0);

                ElementBase* boundaryElement = nullptr;

                int localDimension = elementPtrSlave->GetLocalDimension();

                switch (elementPtrSlave->GetEnumType())
                {
                case Element::eElementType::CONTINUUMELEMENT:
                    switch (localDimension)
                    {
                    case 2:
                    {
                        boundaryElement = new ContinuumContactElement<2, TDimMaster>(&elementPtrSlave->AsContinuumElement2D(), surfaceId, masterElements, rPenalty, rContactAlgorithm);
                        break;
                    }
                    default:
                        throw MechanicsException(__PRETTY_FUNCTION__,"Contact element for Continuum element with dimension " +
                                                 std::to_string(elementPtrSlave->GetLocalDimension()) + "not implemented");
                    }
                    break;
                case Element::eElementType::CONTINUUMELEMENTIGA:
                    switch (localDimension)
                    {
                    case 1:
                    {
                        boundaryElement = new ContinuumContactElement<1, TDimMaster>(&elementPtrSlave->AsContinuumElementIGA1D(), surfaceId, masterElements, rPenalty, rContactAlgorithm);
                        break;
                    }
                    case 2:
                    {
                        boundaryElement = new ContinuumContactElement<2, TDimMaster>(&elementPtrSlave->AsContinuumElementIGA2D(), surfaceId, masterElements, rPenalty, rContactAlgorithm);
                        break;
                    }
                    default:
                        throw MechanicsException(__PRETTY_FUNCTION__,"Contact element for Continuum element with dimension " +
                                                 std::to_string(elementPtrSlave->GetLocalDimension()) + "not implemented");
                    }
                    break;
                default:
                    throw MechanicsException(__PRETTY_FUNCTION__,"Only ContinuumElement and ContinuumElementIGA available as master elements.");
                }

                //find unused integer id
                int elementId = mElementMap.size();
                boost::ptr_map<int, ElementBase>::iterator it = mElementMap.find(elementId);
                while (it != mElementMap.end())
                {
                    elementId++;
                    it = mElementMap.find(elementId);
                }

                mElementMap.insert(elementId, boundaryElement);
                newBoundaryElementIds.push_back(elementId);

                if (rIntegrationType == eIntegrationType::NotSet)
                    throw MechanicsException(__PRETTY_FUNCTION__, "Could not automatically determine integration type of the boundary element.");

                boundaryElement->SetIntegrationType(GetPtrIntegrationType(rIntegrationType), ipDataType);
                boundaryElement->SetConstitutiveLaw(constitutiveLaw);
            }
        }
    }

}

template void NuTo::Structure::ContactElementsCreate<3,2>(int rElementsGroupIDSlave,  int rNodeGroupSlaveId, const Eigen::Matrix<std::pair<int, int>, Eigen::Dynamic, Eigen::Dynamic> &rMasterElementsID, eIntegrationType rIntegrationType, double rPenalty,int rContactAlgorithm);
template void NuTo::Structure::ContactElementsCreate<2,2>(int rElementsGroupIDSlave,  int rNodeGroupSlaveId, const Eigen::Matrix<std::pair<int, int>, Eigen::Dynamic, Eigen::Dynamic> &rMasterElementsID, eIntegrationType rIntegrationType, double rPenalty,int rContactAlgorithm);
template void NuTo::Structure::ContactElementsCreate<2,1>(int rElementsGroupIDSlave,  int rNodeGroupSlaveId, const Eigen::Matrix<std::pair<int, int>, Eigen::Dynamic, Eigen::Dynamic> &rMasterElementsID, eIntegrationType rIntegrationType, double rPenalty,int rContactAlgorithm);
template void NuTo::Structure::ContactElementsCreate<1,2>(int rElementsGroupIDSlave,  int rNodeGroupSlaveId, const Eigen::Matrix<std::pair<int, int>, Eigen::Dynamic, Eigen::Dynamic> &rMasterElementsID, eIntegrationType rIntegrationType, double rPenalty,int rContactAlgorithm);
template void NuTo::Structure::ContactElementsCreate<1,1>(int rElementsGroupIDSlave,  int rNodeGroupSlaveId, const Eigen::Matrix<std::pair<int, int>, Eigen::Dynamic, Eigen::Dynamic> &rMasterElementsID, eIntegrationType rIntegrationType, double rPenalty,int rContactAlgorithm);


//! @brief creates boundary elements and add them to an element group
//! @param rElementGroupId ... group id including the base elements
//! @param rNodeGroupId ... node group id that includes the surface nodes
//! @param rControlNode if not nullptr, then a boundary element with control node will be created
//! @return ... ids of the created boundary element group
int NuTo::Structure::BoundaryElementsCreate(int rElementGroupId, int rNodeGroupId, NodeBase* rControlNode)
{
#ifdef SHOW_TIME
    std::clock_t start, end;
    start = clock();
#endif

    //find groups
    boost::ptr_map<int, GroupBase>::iterator itGroupElements = mGroupMap.find(rElementGroupId);
    if (itGroupElements == mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::BoundaryElementsCreate] Group with the given identifier does not exist.");
    if (itGroupElements->second->GetType() != NuTo::eGroupId::Elements)
        throw MechanicsException("[NuTo::Structure::BoundaryElementsCreate] Group is not an element group.");

    boost::ptr_map<int, GroupBase>::iterator itGroupBoundaryNodes = mGroupMap.find(rNodeGroupId);
    if (itGroupBoundaryNodes == mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::BoundaryElementsCreate] Group with the given identifier does not exist.");
    if (itGroupBoundaryNodes->second->GetType() != NuTo::eGroupId::Nodes)
        throw MechanicsException("[NuTo::Structure::BoundaryElementsCreate] Group is not a node group.");

    Group<ElementBase>& elementGroup    = *(itGroupElements->second->AsGroupElement());
    Group<NodeBase>&    nodeGroup       = *(itGroupBoundaryNodes->second->AsGroupNode());

    // since the search is done via the id's, the surface nodes are ptr, so make another set with the node ptrs
    std::set<const NodeBase*> nodePtrSet;
    for (auto itNode : nodeGroup)
    {
        nodePtrSet.insert(itNode.second);
    }

    std::vector<int> newBoundaryElementIds;

    //loop over all elements
    for (auto itElement : elementGroup)
    {
        ElementBase* elementPtr = itElement.second;
        const InterpolationType* interpolationType = elementPtr->GetInterpolationType();
        try
        {

            //loop over all surfaces
            for (int iSurface = 0; iSurface < interpolationType->GetNumSurfaces(); ++iSurface)
            {
                bool elementSurfaceNodesMatchBoundaryNodes = true;
                Eigen::VectorXi surfaceNodeIndices = interpolationType->GetSurfaceNodeIndices(iSurface);

                int numSurfaceNodes = surfaceNodeIndices.rows();
                std::vector<const NodeBase*> surfaceNodes(numSurfaceNodes);

                for (int iSurfaceNode = 0; iSurfaceNode < numSurfaceNodes; ++iSurfaceNode)
                {
                    surfaceNodes[iSurfaceNode] = elementPtr->GetNode(surfaceNodeIndices(iSurfaceNode, 0));
                }

                //check, if all surface nodes are in the node group
                for (unsigned int countNode = 0; countNode < surfaceNodes.size(); countNode++)
                {
                    if (nodePtrSet.find(surfaceNodes[countNode]) == nodePtrSet.end())
                    {
                        //this surface has at least one node that is not in the list, continue
                        elementSurfaceNodesMatchBoundaryNodes = false;
                        break;
                    }
                }

                if (elementSurfaceNodesMatchBoundaryNodes)
                {
                    int surfaceId = iSurface;

                    //find unused integer id
                    int elementId = mElementMap.size();
                    boost::ptr_map<int, ElementBase>::iterator it = mElementMap.find(elementId);
                    while (it != mElementMap.end())
                    {
                        elementId++;
                        it = mElementMap.find(elementId);
                    }

                    // extract data from base element
                    IpData::eIpDataType ipDataType = elementPtr->GetIpDataType(0);

                    ElementBase* boundaryElement = nullptr;
                    ConstitutiveBase* constitutiveLaw = elementPtr->GetConstitutiveLaw(0);

                    eIntegrationType integrationType = eIntegrationType::NotSet;

                    switch (elementPtr->GetEnumType())
                    {
                    case Element::eElementType::CONTINUUMELEMENT:
                        switch (elementPtr->GetLocalDimension())
                        {
                        case 1:
                        {
                            if(rControlNode == nullptr)
                            {
                                boundaryElement = new ContinuumBoundaryElement<1>(&elementPtr->AsContinuumElement1D(), surfaceId);
                            }
                            else
                            {
                                boundaryElement = new ContinuumBoundaryElementConstrainedControlNode<1>(&elementPtr->AsContinuumElement1D(), surfaceId,rControlNode);
                            }
                            integrationType = eIntegrationType::IntegrationType0DBoundary;
                            break;
                        }
                        case 2:
                        {
                            if(rControlNode == nullptr)
                            {
                                boundaryElement = new ContinuumBoundaryElement<2>(&elementPtr->AsContinuumElement2D(), surfaceId);
                            }
                            else
                            {
                                boundaryElement = new ContinuumBoundaryElementConstrainedControlNode<2>(&elementPtr->AsContinuumElement2D(), surfaceId,rControlNode);
                            }
                            // check for 2D types
                            auto it = std::find(mMappingIntEnum2String.begin(), mMappingIntEnum2String.end(), interpolationType->GetCurrentIntegrationType()->GetStrIdentifier());
                            if (it == mMappingIntEnum2String.end()) break;

                            switch ((eIntegrationType)std::distance(mMappingIntEnum2String.begin(), it)) // Oh my god. Someone please remove this "map" <--- Totally agree with Thomas (vhirtham)
                            {
                            case eIntegrationType::IntegrationType2D3NGauss1Ip:
                            case eIntegrationType::IntegrationType2D4NGauss1Ip:
                                integrationType = eIntegrationType::IntegrationType1D2NGauss1Ip;
                                break;
                            case eIntegrationType::IntegrationType2D3NGauss3Ip:
                            case eIntegrationType::IntegrationType2D4NGauss4Ip:
                                integrationType = eIntegrationType::IntegrationType1D2NGauss2Ip;
                                break;
                            case eIntegrationType::IntegrationType2D3NGauss6Ip:
                            case eIntegrationType::IntegrationType2D4NGauss9Ip:
                                integrationType = eIntegrationType::IntegrationType1D2NGauss3Ip;
                                break;
                            case eIntegrationType::IntegrationType2D3NGauss12Ip:
                                integrationType = eIntegrationType::IntegrationType1D2NGauss5Ip;
                                break;

                            case eIntegrationType::IntegrationType2D4NLobatto9Ip:
                                integrationType = eIntegrationType::IntegrationType1D2NLobatto3Ip;
                                break;

                            case eIntegrationType::IntegrationType2D4NLobatto16Ip:
                                integrationType = eIntegrationType::IntegrationType1D2NLobatto4Ip;
                                break;

                            case eIntegrationType::IntegrationType2D4NLobatto25Ip:
                                integrationType = eIntegrationType::IntegrationType1D2NLobatto5Ip;
                                break;

                                default:
                                    break;
                            }


                            break;
                        }
                        case 3:
                        {
                            if(rControlNode == nullptr)
                            {
                                boundaryElement = new ContinuumBoundaryElement<3>(&elementPtr->AsContinuumElement3D(), surfaceId);
                            }
                            else
                            {
                                boundaryElement = new ContinuumBoundaryElementConstrainedControlNode<3>(&elementPtr->AsContinuumElement3D(), surfaceId,rControlNode);
                            }

                            // check for 3D types
                            auto it = std::find(mMappingIntEnum2String.begin(), mMappingIntEnum2String.end(), interpolationType->GetCurrentIntegrationType()->GetStrIdentifier());
                            if (it == mMappingIntEnum2String.end()) break;

                            switch ((eIntegrationType)std::distance(mMappingIntEnum2String.begin(), it))
                            {
                            case eIntegrationType::IntegrationType3D4NGauss1Ip:
                                integrationType = eIntegrationType::IntegrationType2D3NGauss1Ip;
                                break;

                            case eIntegrationType::IntegrationType3D4NGauss4Ip:
                                integrationType = eIntegrationType::IntegrationType2D3NGauss3Ip;
                                break;

                            case eIntegrationType::IntegrationType3D8NGauss1Ip:
                                integrationType = eIntegrationType::IntegrationType2D4NGauss1Ip;
                                break;

                            case eIntegrationType::IntegrationType3D8NGauss2x2x2Ip:
                                integrationType = eIntegrationType::IntegrationType2D4NGauss4Ip;
                                break;

                            case eIntegrationType::IntegrationType3D8NLobatto3x3x3Ip:
                                integrationType = eIntegrationType::IntegrationType2D4NLobatto9Ip;
                                break;

                            case eIntegrationType::IntegrationType3D8NLobatto4x4x4Ip:
                                integrationType = eIntegrationType::IntegrationType2D4NLobatto16Ip;
                                break;

                            case eIntegrationType::IntegrationType3D8NLobatto5x5x5Ip:
                                integrationType = eIntegrationType::IntegrationType2D4NLobatto16Ip;
                                break;

                            default:
                                    break;
                            }

                            break;
                        }
                        default:
                            throw MechanicsException(__PRETTY_FUNCTION__,"Boundary element for Continuum element with dimension "+
                                                     std::to_string(elementPtr->GetLocalDimension()) + "not implemented");
                        }
                        break;


                    default:

                        break;
                    }

                    mElementMap.insert(elementId, boundaryElement);
                    newBoundaryElementIds.push_back(elementId);

                    if (integrationType == eIntegrationType::NotSet)
                        throw MechanicsException(__PRETTY_FUNCTION__, "Could not automatically determine integration type of the boundary element.");

                    boundaryElement->SetIntegrationType(GetPtrIntegrationType(integrationType), ipDataType);
                    boundaryElement->SetConstitutiveLaw(constitutiveLaw);
                }
            }

        } catch (NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement.second) == itElement.first);
            ss << itElement.first;
            e.AddMessage("[NuTo::Structure::BoundaryElementsCreate] Error creating boundary element for element " + ss.str() + ".");
            throw e;
        } catch (...)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement.second) == itElement.first);
            ss << itElement.first;
            throw NuTo::MechanicsException("[NuTo::Structure::BoundaryElementsCreate] Error creating boundary element for element " + ss.str() + ".");
        }

    }

    int boundaryElementGroup = GroupCreate(eGroupId::Elements);
    for (int boundaryElementId : newBoundaryElementIds)
        GroupAddElement(boundaryElementGroup, boundaryElementId);

#ifdef SHOW_TIME
    end = clock();
    if (mShowTime)
        std::cout << "[NuTo::Structure::BoundaryElementsCreate] " << difftime(end, start) / CLOCKS_PER_SEC << "sec" << std::endl;
#endif

    return boundaryElementGroup;
}

std::pair<int,int> NuTo::Structure::InterfaceElementsCreate(int rElementGroupId, int rInterfaceInterpolationType, int rFibreInterpolationType)
{

    // find element group
    boost::ptr_map<int, GroupBase>::iterator itGroupElements = mGroupMap.find(rElementGroupId);
    if (itGroupElements == mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::BoundaryElementsCreate] Group with the given identifier does not exist.");
    if (itGroupElements->second->GetType() != NuTo::eGroupId::Elements)
        throw MechanicsException("[NuTo::Structure::BoundaryElementsCreate] Group is not an element group.");


    // gets member ids from an element group. The element group must only contain truss elements
    auto elementIds = GroupGetMemberIds(rElementGroupId);

    int groupElementsInterface = GroupCreate(NuTo::eGroupId::Elements);
    int groupElementsFibre = GroupCreate(NuTo::eGroupId::Elements);

    // loop over elements in element group
    for (int i = 0; i < elementIds.size(); ++i)
    {
        auto nodeIds = ElementGetNodes(elementIds(i,0));

        assert( (nodeIds.size() == 2 or nodeIds.size() == 3) and "Only implemented for the 4 node and 6 node interface element");

        NuTo::FullVector<int, Eigen::Dynamic> nodeIdsFibre(nodeIds.size());
        NuTo::FullVector<int, Eigen::Dynamic> nodeIdsMatrix(nodeIds.size());

        // loop over nodes of element
        for (int k = 0; k < nodeIds.size(); ++k)
        {
            FullVector<double, Eigen::Dynamic> nodeCoordinates;
            NodeGetCoordinates(nodeIds(k,0), nodeCoordinates);

            int groupNodes = GroupCreate(NuTo::eGroupId::Nodes);
            GroupAddNodeRadiusRange(groupNodes, nodeCoordinates, 0.0, 1e-6);

            // create an additional node at the same position if it has not been created already
            if (GroupGetNumMembers(groupNodes) > 1)
            {
                assert(GroupGetNumMembers(groupNodes) == 2 and "This group should have exactly two members. Check what went wrong!");
                auto groupNodeMemberIds = GroupGetMemberIds(groupNodes);
                if (groupNodeMemberIds(0,0) == nodeIds(k,0))
                {
                    nodeIdsFibre[k] = groupNodeMemberIds(1,0);

                } else
                {
                    nodeIdsFibre[k] = groupNodeMemberIds(0,0);
                }
            }
            else
            {
                std::set<NuTo::Node::eDof> dofs;
                dofs.insert(NuTo::Node::eDof::COORDINATES);
                dofs.insert(NuTo::Node::eDof::DISPLACEMENTS);

                nodeIdsFibre[k] = NodeCreate(nodeCoordinates, dofs);


            }

            nodeIdsMatrix[k] = nodeIds(k,0);

        }

        // create interface element
        NuTo::FullVector<int, Eigen::Dynamic> nodeIndicesInterface(2 * nodeIds.size());
        for (int iIndex = 0; iIndex < nodeIds.size(); ++iIndex)
        {
            nodeIndicesInterface[iIndex] = nodeIdsMatrix[iIndex];
            nodeIndicesInterface[nodeIndicesInterface.size() - iIndex - 1] = nodeIdsFibre[iIndex];
        }

        int newElementInterface = ElementCreate(rInterfaceInterpolationType, nodeIndicesInterface, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
        GroupAddElement(groupElementsInterface, newElementInterface);

        // create new truss element with duplicated nodes
        int newElementFibre = ElementCreate(rFibreInterpolationType, nodeIdsFibre, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::NOIPDATA);
        GroupAddElement(groupElementsFibre, newElementFibre);

        // delete  old element
        ElementDelete(elementIds(i,0));

    }


    return std::make_pair(groupElementsFibre, groupElementsInterface);

}

//! @brief Deletes a group of elements element
//! @param rGroupNumber group number
void NuTo::Structure::ElementGroupDelete(int rGroupNumber, bool deleteNodes)
{
#ifdef SHOW_TIME
    std::clock_t start, end;
    start = clock();
#endif
    boost::ptr_map<int, GroupBase>::iterator itGroup = mGroupMap.find(rGroupNumber);
    if (itGroup == mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupDelete] Group with the given identifier does not exist.");
    if (itGroup->second->GetType() != NuTo::eGroupId::Elements)
        throw MechanicsException("[NuTo::StructureBase::ElementGroupDelete] Group is not an element group.");

//the group has to be copied, since the elements are removed from this group, which invalidates the iterators
    Group<ElementBase> copyOfElementGroup = Group<ElementBase>(*(itGroup->second->AsGroupElement()));

    std::set<NodeBase*> potentialNodesToBeRemoved;
    for (Group<ElementBase>::iterator itElement = copyOfElementGroup.begin(); itElement != copyOfElementGroup.end(); itElement++)
    {
        try
        {
            //save the nodes, which are eventually to be removed
            if (deleteNodes)
            {
                for (int countNode = 0; countNode < itElement->second->GetNumNodes(); countNode++)
                {
                    NodeBase* nodePtr = itElement->second->GetNode(countNode);
                    potentialNodesToBeRemoved.insert(nodePtr);

                }
            }
            ElementDeleteInternal(itElement->second->ElementGetId());
        } catch (NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second) == itElement->first);
            ss << itElement->first;
            e.AddMessage("[NuTo::StructureBase::ElementGroupDelete] Error deleting element " + ss.str() + ".");
            throw e;
        } catch (...)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second) == itElement->first);
            ss << itElement->first;
            throw NuTo::MechanicsException("[NuTo::StructureBase::ElementGroupDelete] Error deleting element " + ss.str() + ".");
        }
    }

//check all the other elements and see, if they have one of the potential Nodes To Be Removed as valid node
    for (boost::ptr_map<int, ElementBase>::iterator itElement = mElementMap.begin(); itElement != mElementMap.end(); itElement++)
    {
        for (int countNode = 0; countNode < (itElement->second)->GetNumNodes(); countNode++)
        {
            NodeBase* nodePtr = (itElement->second)->GetNode(countNode);
            //int numRemoved = potentialNodesToBeRemoved.erase(nodePtr);
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
#ifdef SHOW_TIME
    end = clock();
    if (mShowTime)
        std::cout << "[NuTo::StructureBase::ElementGroupDelete] " << difftime(end, start) / CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief Deletes an element
//! @param rItElement iterator of the map
void NuTo::Structure::ElementDeleteInternal(int rElementId)
{
// find element
    boost::ptr_map<int, ElementBase>::iterator itElement = mElementMap.find(rElementId);
    if (itElement == this->mElementMap.end())
    {
        throw MechanicsException("[NuTo::Structure::ElementDeleteInternal] Element does not exist.");
    } else
    {
        // Search for elements in groups: using a loop over all groups
        for (boost::ptr_map<int, GroupBase>::iterator groupIt = mGroupMap.begin(); groupIt != mGroupMap.end(); ++groupIt)
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
void NuTo::Structure::GetElementsTotal(std::vector<std::pair<int, const ElementBase*> >& rElements) const
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
void NuTo::Structure::GetElementsTotal(std::vector<std::pair<int, ElementBase*> >& rElements)
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

