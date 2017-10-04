#pragma once


#include "mechanics/structures/StructureBase.h"
#include <set>


namespace NuTo
{


namespace Interpolation
{
enum class eShapeType;
enum class eTypeOrder;
} // namespace Interpolation

//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard class for irregular (unstructured) structures

class Structure : public StructureBase
{

public:
    //! @brief Typedefinitions
    //! @todo check if it is useful to switch to Boost::Ptr container types
    typedef ElementBase* elementBasePtr_t;
    typedef std::vector<elementBasePtr_t> elementBasePtrVec_t;
    typedef std::set<elementBasePtr_t> elementBasePtrSet_t;
    typedef NodeBase* nodeBasePtr_t;
    typedef std::map<nodeBasePtr_t, nodeBasePtr_t> nodeBasePtrMap_t;
    typedef std::set<nodeBasePtr_t> nodeBasePtrSet_t;
    typedef std::vector<nodeBasePtr_t> nodeBasePtrVec_t;

    //! @brief constructor
    //! @param mDimension  Structural dimension (1,2 or 3)
    Structure(int mDimension);

    //! @brief destructor
    virtual ~Structure();


#ifndef SWIG

    //! @brief Calculates the initial value rates (velocities) of the system to meet equilibrium
    virtual void CalculateInitialValueRates(TimeIntegrationBase& rTimeIntegrationScheme) override;

    //! @brief ... evaluates the structure
    virtual void Evaluate(const NuTo::ConstitutiveInputMap& rInput,
                          std::map<eStructureOutput, StructureOutputBase*>& rStructureOutput) override;

#endif


    //*************************************************
    //************ Node routines        ***************
    //***  defined in structures/StructureNode.cpp  ***
    //*************************************************
    //! @brief returns the number of nodes
    //! @return number of nodes
    int GetNumNodes() const override;

    //! @brief returns a reference to a node
    //! @param identifier
    //! @return reference to a node
    NodeBase* NodeGetNodePtr(int rIdent) override;

#ifndef SWIG
    //! @brief returns a reference to a node
    //! @param identifier
    //! @return reference to a node
    const NodeBase* NodeGetNodePtr(int rIdent) const override;

    //! @brief ... store all elements connected to this node in a vector
    //! @param rNode (Input) 		... node pointer
    //! @param rElements (Output) 	... vector of element pointers
    void NodeGetElements(const NodeBase* rNodePtr, std::vector<ElementBase*>& rElements) override;

    //! @brief gives the identifier of a node
    //! @param pointer to a node
    //! @return id
    int NodeGetId(const NodeBase* rNode) const override;

    //! @brief returns const reference to mNodeMap
    //! @return mNodeMap
    const boost::ptr_map<int, NodeBase>& NodeGetNodeMap() const;
#endif // SWIG

    //! @brief ... store all elements connected to this node in a vector
    //! @param rNodeId (Input) 			... node id
    //! @param rElementNumbers (Output)	... vector of element ids
    void NodeGetElements(const int rNodeId, std::vector<int>& rElementNumbers);

    //! @brief creates a node at coordinate's origin
    //! @return node number
    int NodeCreate();

    //! @brief creates a node with coordinates only
    //! @param rCoordinates coordinates of the node
    //! @return node number
    int NodeCreate(Eigen::VectorXd rCoordinates);

    //! @brief creates a node with coordinates only
    //! @param rNodeNumber ... node number
    //! @param rCoordinates ...  node coordinates
    void NodeCreate(int rNodeNumber, Eigen::VectorXd rCoordinates);

    //! @brief creates multiple nodes with coordinates only
    //! @param rCoordinates ...  nodal coordinates (column-wise storage of each nodal coordinate)
    //! @return a std::vector<int> containing the node numbers
    std::vector<int> NodesCreate(Eigen::MatrixXd& rCoordinates);

    //! @brief creates a node with specific dofs at coordinate's origin
    //! @param rDOFs ... space separated string containing the node dofs (e.g. displacements, rotations, temperatures)
    //! @param rCoordinates ...  node coordinates
    //! @return node number
    int NodeCreateDOFs(std::string rDOFs);

    //! @brief creates a node with specific dofs
    //! @param rDOFs ... space separated string containing the node dofs (e.g. displacements, rotations, temperatures)
    //! @param rCoordinates ...  node coordinates
    //! @return node number
    int NodeCreateDOFs(std::string rDOFs, Eigen::VectorXd rCoordinates);

    //! @brief creates a node with specific dofs
    //! @param node number
    //! @param rDOFs ... space separated string containing the node dofs (e.g. displacements, rotations, temperatures)
    //! @param rCoordinates ...  node coordinates
    void NodeCreateDOFs(int rNodeNumber, std::string rDOFs, Eigen::VectorXd rCoordinates);

#ifndef SWIG
    //! @brief creates a node with specific dofs at coordinate's origin
    //! @param rDOFs ... set containing the node dof enums (e.g. displacements, rotations, temperatures)
    //! @param rCoordinates ...  node coordinates
    //! @return node number
    int NodeCreateDOFs(std::set<Node::eDof> rDOFs);

    //! @brief creates a node with specific dofs
    //! @param rDOFs ... set containing the node dof enums (e.g. displacements, rotations, temperatures)
    //! @param rCoordinates ...  node coordinates
    //! @return node number
    int NodeCreateDOFs(std::set<Node::eDof> rDOFs, Eigen::VectorXd rCoordinates);

    //! @brief creates a node with specific dofs
    //! @param node number
    //! @param rDOFs ... set containing the node dof enums (e.g. displacements, rotations, temperatures)
    //! @param rCoordinates ...  node coordinates
    void NodeCreateDOFs(int rNodeNumber, std::set<Node::eDof> rDOFs, Eigen::VectorXd rCoordinates);
#endif // SWIG

    //! @brief deletes a node
    //! @param rNodeNumber ... node number
    void NodeDelete(const int rNodeNumber) override;

    //! @brief info about the nodes in the Structure
    void NodeInfo(int mVerboseLevel) const override;

    //! @brief numbers the dofs in the structure
    //! @param rCallerName ... if the method throws it is nice to know by whom it was called.
    void NodeBuildGlobalDofs(std::string rCallerName = "") override;

    //! @brief extract dof values (e.g. displacements, temperatures to the nodes)
    //! @param rTimeDerivative time derivative (0 disp 1 vel 2 acc)
    //! @return ... StructureBlockVector containing the dofs (J and K)
    virtual NuTo::StructureOutputBlockVector NodeExtractDofValues(int rTimeDerivative) const override;

    virtual void NodeMergeDofValues(int rTimeDerivative, const NuTo::StructureOutputBlockVector& rDofValues) override
    {
        NodeMergeDofValues(rTimeDerivative, rDofValues.J, rDofValues.K);
    }

    virtual void NodeMergeDofValues(NuTo::StructureOutputBlockVector& rDofValues) override
    {
        NodeMergeDofValues(0, rDofValues);
    }

    virtual void NodeMergeDofValues(int rTimeDerivative, const NuTo::BlockFullVector<double>& rActiveDofValues,
                                    const NuTo::BlockFullVector<double>& rDependentDofValues) override;

    //! @brief calculate dependent dof values (for the zeroth time derivative)
    //! @param rActiveDofValues ... vector of independent dof values (ordering according to global dofs, size is number
    //! of active dofs)
    //! @return  ... vector of dependent  dof values (ordering according to global dofs, size is number of active dofs)
    virtual NuTo::BlockFullVector<double>
    NodeCalculateDependentDofValues(const NuTo::BlockFullVector<double>& rActiveDofValues) const override;

    Eigen::Matrix<double, 3, 3> DoubleMatrix(const Eigen::Matrix<double, 3, 3>& matrix)
    {
        return 2.0 * matrix;
    }

#ifndef SWIG


    //! @brief creates a node with rDofs degrees of freedom
    //! @param rCoordinates coordinates of the node
    //! @param rDofs degrees of freedom of the node
    //! @return node number
    int NodeCreate(Eigen::VectorXd rCoordinates, std::set<Node::eDof> rDofs);


    //! @brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
    //! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still
    //! identical
    //! @param rId ... old node id
    //! @param rOldPtr ... old node ptr
    //! @param rNewPtr ... new node ptr
    //! @param rElements (optional) ... vector of all elements that contain the node - speedup!
    void NodeExchangePtr(int rId, NuTo::NodeBase* rOldPtr, NuTo::NodeBase* rNewPtr,
                         std::vector<ElementBase*> rElements = std::vector<ElementBase*>());


#endif // SWIG

    //*************************************************
    //************ Element routines     ***************
    //**  defined in structures/StructureElement.cpp **
    //*************************************************
    //! @brief returns the number of nodes
    //! @return number of nodes
    int GetNumElements() const override;

#ifndef SWIG
    //! @brief returns a reference to an element
    //! @param identifier
    //! @return reference to an element
    ElementBase* ElementGetElementPtr(int rIdent) override;

    //! @brief returns a reference to an element
    //! @param identifier
    //! @return reference to an element
    const ElementBase* ElementGetElementPtr(int rIdent) const override;

    //! @brief gives the identifier of an element
    //! @param pointer to an element
    //! @return identifier
    int ElementGetId(const ElementBase* rElement) const override;

    //! @brief info about one single element
    //! @param rElement (Input) ... pointer to the element
    //! @param rVerboseLevel (Input) ... level of verbosity
    void ElementInfo(const ElementBase* rElement, int rVerboseLevel) const override;

    /// \brief Set an interpolation type for all elements
    /// \param rInterpolationTypeId
    void ElementTotalSetInterpolationType(const int rInterpolationTypeId);
#endif // SWIG

    //! @brief returns a vector with the node ids of an element
    //! @param identifier
    //! @return vector with node ids
    std::vector<int> ElementGetNodes(int rId);

    //! @brief info about the elements in the Structure
    void ElementInfo(int rVerboseLevel) const override;

    //***************************************************//
    //************ ElementCreate routines ***************//
    //***************************************************//

    //---------------------------------------------------------------------------------------------------------------------------
    //! @brief Creates an IGA element, where the knot indices, beside the nodes (control points), are part of the data
    //! structures
    //! @param rInterpolationTypeId ... interpolation type id
    //! @param rNodeNumbers ... node indices
    //! @param rKnots ... knots defining the element
    //! @param rKnotIDs ... starting knot ids of the element in each element direction
    //! @return element id
    int ElementCreate(int rInterpolationTypeId, const Eigen::VectorXi& rNodeNumbers, const Eigen::MatrixXd& rKnots,
                      const Eigen::VectorXi& rKnotIDs);

    //! @brief Creates an IGA element, where the knot indices, beside the nodes (control points), are part of the data
    //! structures
    //! @param rElementNumber ... element number
    //! @param rInterpolationTypeId ... interpolation type id
    //! @param rNodeNumbers  ... node indices
    //! @param rKnots ... knots defining the element
    //! @param rKnotIDs ... starting knot ids of the element in each element direction
    void ElementCreate(int rElementNumber, int rInterpolationTypeId, const Eigen::VectorXi& rNodeNumbers,
                       const Eigen::MatrixXd& rKnots, const Eigen::VectorXi& rKnotIDs);
    //---------------------------------------------------------------------------------------------------------------------------

    //! @brief Creates an element
    //! @param interpolationTypeId ID of the interpolation type
    //! @param rNodeNumbers node indices
    //! @return ID of the newly created element
    int ElementCreate(int interpolationTypeId, const std::vector<int>& rNodeNumbers);

    //! @brief Creates an element
    //! @param elementNumber Number of the element
    //! @param interpolationTypeId Id of the interpolation type
    //! @param rNodeNumbers  Vector of node numbers
    void ElementCreate(int elementNumber, int interpolationTypeId, const std::vector<int>& rNodeNumbers);
    //---------------------------------------------------------------------------------------------------------------------------

    //! @brief creates multiple elements and adds them to an element group
    //! @param rInterpolationTypeId interpolation type id
    //! @param rNodeIdents Identifier for the corresponding nodes (Incidences have to be stored column-wise)
    //! @return index to the new element group<<
    int ElementsCreate(int rInterpolationTypeId, const Eigen::MatrixXi& rNodeNumbers);

    //! @brief changes the node structure to match the interpolation type for all elements
    //! the node merge distance and the box size are calculated from the element sizes
    void ElementTotalConvertToInterpolationType();

    //! @brief changes the node structure to match the interpolation type
    //! the node merge distance and the box size are calculated from the element sizes
    //! @param rGroupNumberElements group for elements (coordinates only) to be converted
    void ElementConvertToInterpolationType(int rGroupNumberElements);

    //! @brief changes the node structure to match the interpolation type for all elements
    //! @param rNodeDistanceMerge Distance of nodes to be joined (should be significantly smaller than the node distance
    //! in the mesh)
    //! @param rMeshSize approximate size of the elements
    void ElementTotalConvertToInterpolationType(double rNodeDistanceMerge, double rMeshSize);

    //! @brief changes the node structure to match the interpolation type
    //! @param rGroupNumberElements group for elements (coordinates only) to be converted
    //! @param rNodeDistanceMerge Distance of nodes to be joined (should be significantly smaller than the node distance
    //! in the mesh)
    //! @param rMeshSize approximate size of the elements
    void ElementConvertToInterpolationType(int rGroupNumberElements, double rNodeDistanceMerge, double rMeshSize);


    //! @brief Deletes an element
    //! @param rElementNumber element number
    void ElementDelete(int rElementNumber) override;

    //! @brief Deletes a group of elements element
    //! @param rGroupNumber group number
    void ElementGroupDelete(int rGroupNumber, bool deleteNodes);

#ifndef SWIG
    //! @brief Deletes an element
    //! @param rItElement iterator of the map
    void ElementDeleteInternal(int rElementId);

    //! @brief Creates an element
    //! @param rInterpolationTypeId interpolation type id
    //! @param rNodes node vector
    int ElementCreate(int rInterpolationTypeId, std::vector<NodeBase*> rNodes);

    //! @brief Creates an element
    //! @param rElementNumber ... element number
    //! @param rInterpolationTypeId interpolation type id
    //! @param rNodes node vector
    void ElementCreate(int rElementNumber, int rInterpolationTypeId, std::vector<NodeBase*> rNodes);

#endif // SWIG

    //! @brief creates boundary elements and add them to an element group
    //! @param rElementGroupId ... group id including the base elements
    //! @param rNodeGroupId ... node group id that includes the surface nodes
    //! @param rControlNode if not nullptr, then a boundary element with control node will be created
    //! @return ... ids of the created boundary element group
    int BoundaryElementsCreate(int rElementGroupId, int rNodeGroupId, NodeBase* rControlNode = nullptr);

    //! @brief creates boundary contact elements and add them to an element group
    //! @param rElementGroupId ... group id including the base elements
    //! @param rNodeGroupId ... node group id that includes the surface nodes
    //! @return ... ids of the created boundary contact element group
    int BoundaryElementsContactCreate(int rElementGroupId, int rNodeGroupId);

    //! @brief  Creates interface elements from an element group.
    //! @param  rElementGroupId: group id including the base elements
    //! @param  rInterfaceInterpolationType: interpolation type of the interface elements
    //! @param  rFibreInterpolationType: interpolation type of the interface elements
    //! @return returns a pair with the group ids of the new fiber and interface elements
    std::pair<int, int> InterfaceElementsCreate(int rElementGroupId, int rInterfaceInterpolationType,
                                                int rFibreInterpolationType);

    //! @brief Import from gmsh
    //!        Creates groups according to gmsh's physical entities and creates an interpolation types for each group
    //! @param rFileName File name
    //! @return Vector of pair, with element.first containing the group id, and element.second the interpolation type id
    std::vector<std::pair<int, int>> ImportFromGmsh(const std::string& rFileName);

    //*************************************************
    //**      InterpolationType routines             **
    //**    defined in structures/                   **
    //**     ...StructureBaseInterpolationType.cpp   **
    //*************************************************

    //! @brief creates a new interpolation type, calls the enum method
    //! @param rShape ... element shape "TRUSS", "TRIANGLE", "QUAD", "TETRAEDER", "BRICK", etc.
    //! @return ... interpolation type id
    int InterpolationTypeCreate(const std::string& rShape);

    //! @brief creates a new interpolation type, calls the enum method
    //! @param rShape ... element shape "TRUSS", "TRIANGLE", "QUAD", "TET", "BRICK", etc.
    //! @return ... interpolation type id
    int InterpolationTypeCreate(NuTo::Interpolation::eShapeType rShape);

    //! @brief sets the integration type for a specific interpolation type
    //! @param rInterpolationTypeId ... interpolation type id
    //! @param rIntegrationType ... integration type string
    void InterpolationTypeSetIntegrationType(int rInterpolationTypeId, const std::string& rIntegrationType);

    //! @brief prints the info to the interpolation type
    //! @param rInterpolationTypeId ... interpolation type id
    void InterpolationTypeInfo(int rInterpolationTypeId) const;

    //! @brief adds a dof to a interpolation type, calls the enum method
    //! @param rInterpolationTypeId ... interpolation type id
    //! @param rDofType ... dof type
    //! @param rTypeOrder ... type and order of interpolation
    void InterpolationTypeAdd(int rInterpolationTypeId, const std::string& rDofType, const std::string& rTypeOrder);

#ifndef SWIG

    //! @brief creates a new interpolation type
    //! @param rInterpolationTypeId ... interpolation type id
    //! @param rShape ... element shape "1DTRUSS", "2DTRIANGLE", "2DQUAD", "3DTET", "3DBRICK", etc.
    void InterpolationTypeCreate(int rInterpolationTypeId, NuTo::Interpolation::eShapeType rShape);

    //! @brief sets the integration type for a specific interpolation type
    //! @param rInterpolationTypeId ... interpolation type id
    //! @param rIntegrationType ... integration type enum
    //! @param rIpDataType ... ip data type enum
    void InterpolationTypeSetIntegrationType(int rInterpolationTypeId, eIntegrationType rIntegrationType);

    //! @brief sets the integration type for a specific interpolation type
    //! @param rInterpolationTypeId ... interpolation type id
    //! @param rIntegrationType ... integration type pointer
    void InterpolationTypeSetIntegrationType(int rInterpolationTypeId, IntegrationTypeBase* rIntegrationType);

    //! @brief adds a dof to a interpolation type
    //! @param rInterpolationTypeId ... interpolation type id
    //! @param rDofType ... dof type
    //! @param rTypeOrder ... type and order of interpolation
    void InterpolationTypeAdd(int rInterpolationTypeId, NuTo::Node::eDof rDofType,
                              NuTo::Interpolation::eTypeOrder rTypeOrder);


    //! @brief adds a dof to a interpolation type
    //! @param rInterpolationTypeId ... interpolation type id
    //! @param rDofType ... dof type
    //! @param rTypeOrder ... type and order of interpolation
    void InterpolationTypeAdd(int rInterpolationTypeId, NuTo::Node::eDof rDofType, Interpolation::eTypeOrder rTypeOrder,
                              const Eigen::VectorXi& rDegree, const std::vector<Eigen::VectorXd>& rKnots,
                              const Eigen::MatrixXd& rWeights);

    //! @brief returns the interpolation type for a given interpolation type id
    //! @param rInterpolationTypeId ... interpolation type id
    //! @return InterpolationType ptr
    InterpolationType* InterpolationTypeGet(int rInterpolationTypeId);

    //! @brief returns the interpolation type for a given interpolation type id
    //! @param rInterpolationTypeId ... interpolation type id
    //! @return InterpolationType ptr
    const InterpolationType* InterpolationTypeGet(int rInterpolationTypeId) const;

    //! @brief creates a node
    //! @param rDOFs
    //! @param rCoordinates coordinates of the node
    //! @return node pointer
    NodeBase* NodePtrCreate(std::set<Node::eDof> rDOFs, Eigen::VectorXd rCoordinates);

    //! @brief copy and move the structure
    //! most of the data is kept, but e.g. nonlocal data and
    //! @param rOffset offset (dimension x 1 has to be identical with structure dimension)
    //! @param rOld2NewNodePointer ptrMap showing the new and old node pointers
    //! @param rOld2NewElementPointer ptrMap showing the new and old element pointers
    void CopyAndTranslate(Eigen::VectorXd& rOffset, std::map<NodeBase*, NodeBase*>& rOld2NewNodePointer,
                          std::map<ElementBase*, ElementBase*>& rOld2NewElementPointer);
#endif // SWIG

    //! @brief copy and move the structure
    //! most of the data is kept, but e.g. nonlocal data and
    //! @param rOffset offset (dimension x 1 has to be identical with structure dimension)
    void CopyAndTranslate(Eigen::VectorXd& rOffset);

    //! @brief ... Adds an element to an element group
    //! @param rIdentGroup identifier for the group
    //! @param rIdentElement  identifier for the element
    void GroupAddElement(int rIdentGroup, int rIdElement);

    //! @brief ... Adds all elements to an element group
    //! @param rIdentGroup identifier for the group
    void GroupAddElementsTotal(int rIdentGroup);

    //! @brief ... Adds all elements to a group based on the type
    //! @param rIdentGroup identifier for the group
    //! @param rInterpolationType  identifier for the interpolation type
    void GroupAddElementFromType(int rIdentGroup, int rInterpolationType);

    //! @brief ... Adds all nodes in rSearchIdenGroup to rIdentGroup whose coordinates are in the specified range
    //! @param rIdentNodeGroup identifier for the group
    //! @param rSearchIdentElementGroup identifier for the group
    //! @param rDirection either 0,1,2 for x,y, or z
    //! @param rMin ... minimum value
    //! @param rMax ... maximum value
    void GroupAddNodeFromElementGroupCoordinateRange(int rIdentNodeGroup, int rSearchIdentElementGroup, int rDirection,
                                                     double rMin, double rMax);

    //! @brief adds all elements to an element group and returns its id
    int GroupGetElementsTotal();

    //! @brief adds all ndoes to an element group and returns its id
    int GroupGetNodesTotal();


    //! @brief defines the serialization of this class
    //! @param rStream serialize output stream
    virtual void NuToSerializeSave(SerializeStreamOut& rStream) override;

    //! @brief defines the serialization of this class
    //! @param rStream serialize input stream
    virtual void NuToSerializeLoad(SerializeStreamIn& rStream) override;


    //*************************************************
    //************ Info routine         ***************
    //**  defined in structures/Structure.cpp *********
    //*************************************************
    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info() const override;

//*************************************************
//************ DOF routine          ***************
//**  defined in structures/Structure.cpp *********
//*************************************************
#ifndef SWIG
    bool IsActiveDofId(int rDofId, Node::eDof rDofType) const
    {
        return (rDofId < GetNumActiveDofs(rDofType));
    }
#endif // SWIG

protected:
#ifndef SWIG
#endif

#ifndef SWIG

    //! @brief ... store all elements of a structure in a vector
    //! @param rElements ... vector of element pointer
    void GetElementsTotal(std::vector<const ElementBase*>& rElements) const override;

    //! @brief ... store all elements of a structure in a vector
    //! @param rElements ... vector of element pointer
    void GetElementsTotal(std::vector<std::pair<int, const ElementBase*>>& rElements) const override;

    //! @brief ... store all elements of a structure in a vector
    //! @param rElements ... vector of element pointer
    void GetElementsTotal(std::vector<ElementBase*>& rElements) override;

    //! @brief ... store all elements of a structure in a vector
    //! @param rElements ... vector of element pointer
    void GetElementsTotal(std::vector<std::pair<int, ElementBase*>>& rElements) override;

    //! @brief ... store all nodes of a structure in a vector
    //! @param rNodes ... vector of element pointer
    void GetNodesTotal(std::vector<const NodeBase*>& rNodess) const override;

    //! @brief ... store all nodes of a structure in a vector
    //! @param rNodes ... vector of element pointer
    void GetNodesTotal(std::vector<std::pair<int, const NodeBase*>>& rNodes) const override;

    //! @brief ... store all nodes of a structure in a vector
    //! @param rNodes ... vector of element pointer
    void GetNodesTotal(std::vector<NodeBase*>& rNodes) override;

    //! @brief ... store all nodes of a structure in a vector
    //! @param rNodes ... vector of element pointer
    void GetNodesTotal(std::vector<std::pair<int, NodeBase*>>& rNodes) override;
#endif

    //! @brief deletes a node
    //! @param rNodeNumber ... node number
    //! @param checkElements ... check the elements, if set to false, make sure that the node is not part of any element
    void NodeDelete(int rNodeNumber, bool checkElements);

    boost::ptr_map<int, NodeBase> mNodeMap;
    boost::ptr_map<int, ElementBase> mElementMap;
};
} // namespace NuTo
