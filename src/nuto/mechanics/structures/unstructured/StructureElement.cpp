// $Id$

#include <assert.h>
#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/elements/IpDataEnum.h"

#include "nuto/mechanics/elements/Element1D.h"
#include "nuto/mechanics/elements/Element2D.h"
#include "nuto/mechanics/elements/Element3D.h"

#include "nuto/mechanics/elements/BoundaryElement1D.h"
#include "nuto/mechanics/elements/BoundaryElement2D.h"

#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/nodes/NodeDof.h"

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
        std::cout << "\tenum::type=" << rElement->GetEnumType() << std::endl;
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
                    std::cout << "\t\t" << iIp << ": [" << coor[0] << ";" << coor[1] << ";" << coor[2] << "]" << std::endl;
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
    bool showTime = mShowTime;
    mShowTime = false;

    // create element group containing the new elements
    int groupNumber = GroupCreate("Elements");
    GroupBase* group = mGroupMap.find(groupNumber)->second;

    for (boost::ptr_map<int, ElementBase>::iterator itElement = mElementMap.begin(); itElement != mElementMap.end(); ++itElement)
    {
        int elementId = itElement->first;
        ElementBase* elementPtr = itElement->second;
        group->AddMember(elementId, elementPtr);
    }
    mShowTime = showTime;

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

    if (itGroup->second->GetType() != NuTo::Groups::Elements)
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

    std::sort(elementSize.data(), elementSize.data()+numElements, std::less<double>()); // sort in accending order

    double sizeMin = elementSize(0);
    double sizeMed = elementSize(numElements/2);

    double lengthMin = std::pow(sizeMin, 1./mDimension);
    double lengthMed = std::pow(sizeMed, 1./mDimension);

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

    ElementTotalConvertToInterpolationType(mergeDist, meshSize);

}

//! @brief changes the node structure to match the interpolation type for all elements
//! @param rNodeDistanceMerge Distance of nodes to be joined (should be significantly smaller than the node distance in the mesh)
//! @param rMeshSize approximate size of the elements
void NuTo::Structure::ElementTotalConvertToInterpolationType(double rNodeDistanceMerge, double rMeshSize)
{
    // create a group with all elements
    bool showTime = mShowTime;
    mShowTime = false;

    // create element group containing the new elements
    int groupNumber = GroupCreate("Elements");
    GroupBase* group = mGroupMap.find(groupNumber)->second;

    for (boost::ptr_map<int, ElementBase>::iterator itElement = mElementMap.begin(); itElement != mElementMap.end(); ++itElement)
    {
        int elementId = itElement->first;
        ElementBase* elementPtr = itElement->second;
        group->AddMember(elementId, elementPtr);
    }
    mShowTime = showTime;

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

    if (itGroup->second->GetType() != NuTo::Groups::Elements)
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

    // find the bounding box corners = max/min x,y,z coordinates, initialize with random node coordinates
    Eigen::VectorXd boundingBoxMax = elements[0]->GetNode(0)->GetCoordinates();
    Eigen::VectorXd boundingBoxMin = elements[0]->GetNode(0)->GetCoordinates();

    for (ElementBase* element : elements)
    {
        // loop through nodes with coordinates
        for (int iNode = 0; iNode < element->GetNumNodes(Node::COORDINATES); ++iNode)
        {
            NodeBase* node = element->GetNode(iNode, Node::COORDINATES);
            Eigen::VectorXd nodeCoordinates = node->GetCoordinates();
            assert(nodeCoordinates.rows() == mDimension);
            boundingBoxMax = boundingBoxMax.cwiseMax(nodeCoordinates);
            boundingBoxMin = boundingBoxMin.cwiseMin(nodeCoordinates);
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

        for (int iNode = 0; iNode < interpolationType->GetNumNodes(); ++iNode)
        {
            Eigen::VectorXd naturalNodeCoordinates = interpolationType->GetNaturalNodeCoordinates(iNode);
            Eigen::VectorXd globalNodeCoordinates = element->InterpolateDof(naturalNodeCoordinates, Node::COORDINATES);
            Eigen::Vector3i boxIndex3D = Eigen::Vector3i::Zero(); // default access is the first box
            for (int iDim = 0; iDim < mDimension; ++iDim)
                boxIndex3D[iDim] = std::floor((globalNodeCoordinates[iDim] - boundingBoxMin[iDim]) / deltaBox[iDim]);

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
            std::set<Node::eAttributes> nodeDofs;
            for (TmpNode& singleUniqueNode : singleSameNode)
            {
                auto nodeDofsOther =  singleUniqueNode.element->GetInterpolationType()->GetNodeDofs(singleUniqueNode.elementNodeId);
                for (auto dof : nodeDofsOther)
                    nodeDofs.insert(dof);
            }

            NuTo::FullVector<double, Eigen::Dynamic> nodeCoordinates(tmpNode.coords);

            if (nodeDofs.find(Node::COORDINATES) != nodeDofs.end())
            {
                // If the node is a coordinate node, it should already be in the elements mNodes, check that:
                NodeBase* oldNode = singleSameNode[0].element->GetNode(singleSameNode[0].elementNodeId);
                assert(oldNode != nullptr);
                for (TmpNode& singleUniqueNode : singleSameNode)
                    assert(singleUniqueNode.element->GetNode(singleUniqueNode.elementNodeId) == oldNode);

                NodeBase* newNode = NodePtrCreate(nodeDofs, nodeCoordinates);

                // exchange the pointer in all elements, all groups, all constraints etc
                int oldNodeIndex = NodeGetId(oldNode);
                NodeExchangePtr(oldNodeIndex, oldNode, newNode); // from now on, newNode is nullPtr

            }
            else
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

//! @brief Creates an element
//! @param rInterpolationTypeId interpolation type id
//! @param rNodeNumbers node indices
//! @return int rElementNumber
int NuTo::Structure::ElementCreate(int rInterpolationTypeId, const NuTo::FullVector<int, Eigen::Dynamic>& rNodeNumbers)
{
    return ElementCreate(rInterpolationTypeId, rNodeNumbers, "CONSTITUTIVELAWIP", "NOIPDATA");
}

//! @brief Creates an element
//! @param rInterpolationTypeId interpolation type id
//! @param rNodeNumbers node indices
//! @param rElementDataType Element data for the elements
//! @param rIpDataType Integration point data for the elements
//! @return int rElementNumber
int NuTo::Structure::ElementCreate(int rInterpolationTypeId, const NuTo::FullVector<int, Eigen::Dynamic>& rNodeNumbers, const std::string& rElementDataType, const std::string& rIpDataType)
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

//! @brief Creates an element
//! @param rElementNumber element number
//! @param rInterpolationTypeId interpolation type id
//! @param rElementType element type
//! @param rNodeIdents pointers to the corresponding nodes
void NuTo::Structure::ElementCreate(int rElementNumber, int rInterpolationTypeId, const NuTo::FullVector<int, Eigen::Dynamic>& rNodeNumbers, const std::string& rElementDataType, const std::string& rIpDataType)
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

//! @brief Creates an element
//! @param rInterpolationTypeId interpolation type id
//! @param rNodeNumbers node indices
//! @param rElementDataType Element data for the elements
//! @param rIpDataType Integration point data for the elements
//! @return int rElementNumber
int NuTo::Structure::ElementCreate(int rInterpolationTypeId, const NuTo::FullVector<int, Eigen::Dynamic>& rNodeNumbers, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType)
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

    if (interpolationType->IsDof(Node::COORDINATES) == false)
        throw NuTo::MechanicsException("[NuTo::Structure::ElementCreate] COORDINATE interpolation required.");

    unsigned int numNodesCoordinates = interpolationType->Get(Node::COORDINATES).GetNumNodes();
    if (numNodesCoordinates != rNodeVector.size())
        throw NuTo::MechanicsException("[NuTo::Structure::ElementCreate] COORDINATE interpolation requires " + std::to_string(numNodesCoordinates) + " nodes. " + std::to_string(rNodeVector.size()) + " are provided.");

    ElementBase* ptrElement = 0;

    if (interpolationType->GetCurrentIntegrationType() == nullptr)
    {
        IntegrationType::eIntegrationType integrationTypeEnum = interpolationType->GetStandardIntegrationType();
        IntegrationTypeBase* integrationType = GetPtrIntegrationType(integrationTypeEnum);
        interpolationType->UpdateIntegrationType(*integrationType);
    }

    try
    {
        switch (mDimension)
        {
        case 1:
            ptrElement = new Element1D(this, rNodeVector, rElementDataType, rIpDataType, interpolationType);
            break;
        case 2:
            ptrElement = new Element2D(this, rNodeVector, rElementDataType, rIpDataType, interpolationType);
            break;
        case 3:
            ptrElement = new Element3D(this, rNodeVector, rElementDataType, rIpDataType, interpolationType);
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

    bool showTime = mShowTime;
    mShowTime = false;

    // create element group containing the new elements
    int newElementGroup = GroupCreate("Elements");
    for (int newElementId : newElementIds)
        GroupAddElement(newElementGroup, newElementId);

    mShowTime = showTime;
    return newElementGroup;
}

//! @brief creates boundary elements
//! @param rElementGroupId ... group id including the base elements
//! @param rNodeGroupId ... node group id that includes the surface nodes
//! @return ... ids of the created boundary elements
NuTo::FullVector<int, Eigen::Dynamic> NuTo::Structure::BoundaryElementsCreate(int rElementGroupId, int rNodeGroupId)
{
#ifdef SHOW_TIME
    std::clock_t start, end;
    start = clock();
#endif

    //find groups
    boost::ptr_map<int, GroupBase>::iterator itGroupElements = mGroupMap.find(rElementGroupId);
    if (itGroupElements == mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::BoundaryElementsCreate] Group with the given identifier does not exist.");
    if (itGroupElements->second->GetType() != NuTo::Groups::Elements)
        throw MechanicsException("[NuTo::Structure::BoundaryElementsCreate] Group is not an element group.");

    boost::ptr_map<int, GroupBase>::iterator itGroupBoundaryNodes = mGroupMap.find(rNodeGroupId);
    if (itGroupBoundaryNodes == mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::BoundaryElementsCreate] Group with the given identifier does not exist.");
    if (itGroupBoundaryNodes->second->GetType() != NuTo::Groups::Nodes)
        throw MechanicsException("[NuTo::Structure::BoundaryElementsCreate] Group is not a node group.");

    Group<ElementBase>& elementGroup = *(itGroupElements->second->AsGroupElement());
    Group<NodeBase>& nodeGroup = *(itGroupBoundaryNodes->second->AsGroupNode());

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
        const InterpolationType* InterpolationType = elementPtr->GetInterpolationType();

        try
        {

            //loop over all surfaces
            for (int iSurface = 0; iSurface < InterpolationType->GetNumSurfaces(); iSurface++)
            {
                bool elementSurfaceNodesMatchBoundaryNodes = true;
                Eigen::VectorXi surfaceNodeIndices = InterpolationType->GetSurfaceNodeIndices(iSurface);

                int numSurfaceNodes = surfaceNodeIndices.rows();
                std::vector<const NodeBase*> surfaceNodes(numSurfaceNodes);

                for (int iSurfaceNode = 0; iSurfaceNode < numSurfaceNodes; ++iSurfaceNode)
                {
                    surfaceNodes[iSurfaceNode] = elementPtr->GetNode(surfaceNodeIndices.at(iSurfaceNode, 0));
                }

                //check, if all surface nodes are in the node group
                for (unsigned int countNode = 0; countNode < surfaceNodes.size(); countNode++)
                {
                    if (nodePtrSet.find(surfaceNodes[countNode]) == nodePtrSet.end())
                    {
                        //this surface has at least one node that is not in the list, continue
                        elementSurfaceNodesMatchBoundaryNodes = false;
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

                    IntegrationType::eIntegrationType integrationType;

                    switch (elementPtr->GetEnumType())
                    {
                    case Element::ELEMENT1D:
                        // create BoundaryElement1D

                        boundaryElement = new BoundaryElement1D(elementPtr, surfaceId);
                        integrationType = IntegrationType::IntegrationType0DBoundary;
                        break;
                    case Element::ELEMENT2D:
                        boundaryElement = new BoundaryElement2D(elementPtr, surfaceId);
                        integrationType = IntegrationType::IntegrationType1D2NGauss3Ip; // TODO, esp. for 3D. Maybe InterpolationType::GetSurfaceInterpolationType
                        break;
                    case Element::ELEMENT3D:
                        // create BoundaryElement3D
                        throw MechanicsException("[NuTo::Structure::BoundaryElementsCreate] not yet implemented for 3D.");
                        break;
                    default:
                        break;
                    }

                    mElementMap.insert(elementId, boundaryElement);
                    newBoundaryElementIds.push_back(elementId);

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

    NuTo::FullVector<int, Eigen::Dynamic> boundaryElementIds(newBoundaryElementIds);

#ifdef SHOW_TIME
    end = clock();
    if (mShowTime)
        std::cout << "[NuTo::Structure::BoundaryElementsCreate] " << difftime(end, start) / CLOCKS_PER_SEC << "sec" << std::endl;
#endif

    return boundaryElementIds;
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
    if (itGroup->second->GetType() != NuTo::Groups::Elements)
        throw MechanicsException("[NuTo::StructureBase::ElementGroupDelete] Group is not an element group.");

//the group has to be copied, since the elements are removed from this group, which invalidates the iterators
    Group<ElementBase> copyOfElementGroup = *(dynamic_cast<Group<ElementBase>*>(itGroup->second));

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

    for (std::set<NodeBase*>::iterator itNode = potentialNodesToBeRemoved.begin(); itNode != potentialNodesToBeRemoved.end(); itNode++)
    {
        int nodeId(NodeGetId(*itNode));
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
    }
    else
    {
        // Search for elements in groups: using a loop over all groups
        for (boost::ptr_map<int, GroupBase>::iterator groupIt = mGroupMap.begin(); groupIt != mGroupMap.end(); ++groupIt)
        {
            if (groupIt->second->GetType() == NuTo::Groups::Elements)
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

