// $Id$

#ifndef STRUCTURE_H
#define STRUCTURE_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#endif // ENABLE_SERIALIZATION
#include <set>

#include <boost/ptr_container/ptr_map.hpp>

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/structures/StructureBaseEnum.h"
#include "nuto/mechanics/cracks/CrackBase.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/ElementDataEnum.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/nodes/NodeBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard class for irregular (unstructured) structures
class Structure : public StructureBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief Typedefinitions
    //! @todo check if it is useful to switch to Boost::Ptr container types
	typedef ElementBase* elementBasePtr_t;
	typedef std::vector<elementBasePtr_t> elementBasePtrVec_t;
	typedef std::set<elementBasePtr_t>    elementBasePtrSet_t;
	typedef NodeBase* nodeBasePtr_t;
	typedef std::map<nodeBasePtr_t , nodeBasePtr_t > nodeBasePtrMap_t;
	typedef std::set<nodeBasePtr_t>    nodeBasePtrSet_t;
	typedef std::vector<nodeBasePtr_t> nodeBasePtrVec_t;
	typedef CrackBase* crackBasePtr_t;
	typedef std::vector<crackBasePtr_t> crackBasePtrVec_t;
    typedef boost::ptr_map<int,CrackBase> crackMap_t;

    //! @brief constructor
    //! @param mDimension  Structural dimension (1,2 or 3)
    Structure(int mDimension);

#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
    //! @brief serializes the class, this is the load routine
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void load(Archive & ar, const unsigned int version);

    //! @brief serializes the class, this is the save routine
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const;

    BOOST_SERIALIZATION_SPLIT_MEMBER()

#endif// SWIG

    //! @brief ... save the object to a file
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
    void Save (const std::string &filename, std::string rType )const;

    //! @brief ... restore the object from a file
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
    void Restore (const std::string &filename, std::string rType );
#endif // ENABLE_SERIALIZATION

    //! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    std::string GetTypeId()const
    {
        return std::string("Structure");
    }

    //! @brief Builds the nonlocal data for integral type nonlocal constitutive models
    //! @param rConstitutiveId constitutive model for which the data is build
    void BuildNonlocalData(int rConstitutiveId);

#ifndef SWIG
    //! @brief Builds the nonlocal data for integral type nonlocal constitutive models
    //! @param rConstitutiveId constitutive model for which the data is build
    void BuildNonlocalData(const ConstitutiveBase* rConstitutive);
#endif //SWIG

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    NuTo::Error::eError BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureBaseEnum::eMatrixType rType, NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK);

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    //! @param rMatrixKJ ... submatrix kj (number of dependent dof x number of active dof)
    //! @param rMatrixKK ... submatrix kk (number of dependent dof x number of dependent dof)
    NuTo::Error::eError BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureBaseEnum::eMatrixType rType, NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK, NuTo::SparseMatrix<double>& rMatrixKJ, NuTo::SparseMatrix<double>& rMatrixKK);

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    NuTo::Error::eError BuildGlobalCoefficientSubMatricesSymmetric(NuTo::StructureBaseEnum::eMatrixType rType, NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK);

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    //! @param rMatrixKK ... submatrix kk (number of dependent dof x number of dependent dof)
    NuTo::Error::eError BuildGlobalCoefficientSubMatricesSymmetric(NuTo::StructureBaseEnum::eMatrixType rType, NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK, NuTo::SparseMatrix<double>& rMatrixKK);

    //! @brief ... based on the global dofs build sub-vectors of the global lumped mass
    //! @param rActiveDofVector ... global lumped mass which corresponds to the active dofs
    //! @param rDependentDofVector ... global lumped mass which corresponds to the dependent dofs
    Error::eError BuildGlobalLumpedHession2(NuTo::FullVector<double,Eigen::Dynamic>& rActiveDofVector,
    		NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofVector);

    //! @brief ... based on the global dofs build sub-vectors of the global internal potential gradient
    //! @param rActiveDofGradientVector ... global internal potential gradient which corresponds to the active dofs
    //! @param rDependentDofGradientVector ... global internal potential gradient which corresponds to the dependent dofs
    //! @param rUpdateHistoryVariables (update history variables after having calculated the response)
    NuTo::Error::eError BuildGlobalGradientInternalPotentialSubVectors(NuTo::FullVector<double,Eigen::Dynamic>& rActiveDofGradientVector,
    		NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofGradientVector, bool rUpdateHistoryVariables);

//*************************************************
//************ Node routines        ***************
//***  defined in structures/StructureNode.cpp  ***
//*************************************************
    //! @brief returns the number of nodes
    //! @return number of nodes
    int GetNumNodes() const;

#ifndef SWIG
    //! @brief returns a reference to a node
    //! @param identifier
    //! @return reference to a node
    NodeBase* NodeGetNodePtr(int rIdent);

    //! @brief returns a reference to a node
    //! @param identifier
    //! @return reference to a node
    const NodeBase* NodeGetNodePtr(int rIdent)const;

    //! @brief ... store all elements connected to this node in a vector
    //! @param rNode (Input) 		... node pointer
    //! @param rElements (Output) 	... vector of element pointers
    void NodeGetElements(const NodeBase* rNodePtr, std::vector<ElementBase*>& rElements);

    //! @brief gives the identifier of a node
    //! @param pointer to a node
    //! @return id
    int NodeGetId(const NodeBase* rNode)const;
#endif //SWIG

    //! @brief ... return the global dof number of the displacement component of a node
    //! @param rNodeId (Input) 			... node id
    //! @param rDispDof 	... local disp dof (0,1 or 2 for x,y or z)
    //! @returnrglobal dof number
    int NodeGetDofDisplacement(int rNodeId, int rDispDof);

    //! @brief ... store all elements connected to this node in a vector
    //! @param rNodeId (Input) 			... node id
    //! @param rElementNumbers (Output)	... vector of element ids
    void NodeGetElements(const int rNodeId, NuTo::FullVector<int,Eigen::Dynamic>& rElementNumbers);

    //! @brief creates a node at coordinate's origin
    //! @param rDOFs space separated string containing the attributes like
    //! DISPLACEMENTS, ROTATIONS; TEMPERATURES
    //! @return node number
    int NodeCreate(std::string rDOFs);

    //! @brief creates a node
    //! @param rDOFs space separated string containing the attributes like
    //! DISPLACEMENTS, ROTATIONS; TEMPERATURES
    //! the number of DOFs will be chosen according to the structure
    //! e.g. a 2D structure might have 2 displacements and 1 rotation
    //! @param rCoordinates coordinates of the node
    //! @return node number
    int NodeCreate(std::string rDOFs, NuTo::FullVector<double,Eigen::Dynamic>& rCoordinates);

    //! @brief creates a node
    //! @param rDOFs space separated string containing the attributes like
    //! DISPLACEMENTS, ROTATIONS; TEMPERATURES
    //! the number of DOFs will be chosen according to the structure
    //! e.g. a 2D structure might have 2 displacements and 1 rotation
    //! @param rCoordinates coordinates of the node
    //! @param rNumTimeDerivatives ...  number of time derivatives for the dofs to be stored
    //! @return node number
    int NodeCreate(std::string rDOFs, NuTo::FullVector<double,Eigen::Dynamic>& rCoordinates, int rNumTimeDerivatives);

    //! @brief creates a node
    //! @param rNodeNumber ... node number
    //! @param rDOFs ... space separated string containing the node dofs (e.g. displacements, rotations, temperatures)
    //! @param rCoordinates ...  node coordinates
    inline void NodeCreate(int rNodeNumber, std::string rDOFs, NuTo::FullVector<double,Eigen::Dynamic>& rCoordinates)
    {
    	NodeCreate(rNodeNumber,rDOFs,rCoordinates,0);
    }

    //! @brief creates a node
    //! @param rNodeNumber ... node number
    //! @param rDOFs ... space separated string containing the node dofs (e.g. displacements, rotations, temperatures)
    //! @param rCoordinates ...  node coordinates
    //! @param rNumTimeDerivatives ...  number of time derivatives for the dofs to be stored
    void NodeCreate(int rNodeNumber, std::string rDOFs, NuTo::FullVector<double,Eigen::Dynamic>& rCoordinates, int rNumTimeDerivatives);

    //! @brief creates multiple nodes
    //! @param rDOFs ... space separated string containing the nodal dofs (e.g. displacements, rotations, temperatures)
    //! @param rCoordinates ...  nodal coordinates (column-wise storage of each nodal coordinate)
    //! @return a NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> containing the node numbers
    NuTo::FullVector<int,Eigen::Dynamic> NodesCreate(std::string rDOFs, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoordinates);

    //! @brief deletes a node
    //! @param rNodeNumber ... node number
    void NodeDelete(const int rNodeNumber);

    //! @brief info about the nodes in the Structure
    void NodeInfo(int mVerboseLevel)const;

    //! @brief numbers the dofs in the structure
    void NodeBuildGlobalDofs();

    //! @brief write dof values (e.g. displacements, temperatures to the nodes)
    //! @param rTimeDerivative time derivative (0 disp 1 vel 2 acc)
    //! @param rActiveDofValues ... vector of global dof values (ordering according to global dofs, size is number of active dofs)
    virtual void NodeMergeActiveDofValues(int rTimeDerivative, const NuTo::FullVector<double,Eigen::Dynamic>& rActiveDofValues);

    using StructureBase::NodeMergeActiveDofValues;

    //! @brief write dof values (e.g. displacements, temperatures to the nodes)
    //! @param rTimeDerivative time derivative (0 disp 1 vel 2 acc)
    //! @param rActiveDofValues ... vector of independent dof values (ordering according to global dofs, size is number of active dofs)
    //! @param rDependentDofValuess ... vector of global dependent values (ordering according to global dofs, size is number of active dofs)
    virtual void NodeMergeDofValues(int rTimeDerivative, const NuTo::FullVector<double,Eigen::Dynamic>& rActiveDofValues,const NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofValuess);

    using StructureBase::NodeMergeDofValues;

    //! @brief ...merge additional dof values
    //! @param rTimeDerivative time derivative (0 disp 1 vel 2 acc)
    //! empty in the standard case, is used in StructureIp
    virtual void NodeMergeAdditionalGlobalDofValues(int rTimeDerivative, const NuTo::FullVector<double,Eigen::Dynamic>& rActiveDofValues, const NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofValues){};

    //! @brief extract dof values (e.g. displacements, temperatures to the nodes)
    //! @param rTimeDerivative time derivative (0 disp 1 vel 2 acc)
    //! @param rActiveDofValues ... vector of global active dof values (ordering according to global dofs, size is number of active dofs)
    //! @param rDependentDofValues ... vector of global dependent dof values (ordering according to (global dofs) - (number of active dofs), size is (total number of dofs) - (number of active dofs))
    virtual void NodeExtractDofValues(int rTimeDerivative, NuTo::FullVector<double,Eigen::Dynamic>& rActiveDofValues, NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofValues) const;

    using StructureBase::NodeExtractDofValues;

    //! @brief extract dof values additional dof values
    //! @param rTimeDerivative time derivative (0 disp 1 vel 2 acc)
    //! @param rActiveDofValues ... vector of global active dof values (ordering according to global dofs, size is number of active dofs)
    //! @param rDependentDofValues ... vector of global dependent dof values (ordering according to (global dofs) - (number of active dofs), size is (total number of dofs) - (number of active dofs))
    virtual void NodeExtractAdditionalGlobalDofValues(int rTimeDerivative, NuTo::FullVector<double,Eigen::Dynamic>& rActiveDofValues, NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofValues) const{};

#ifndef SWIG
    //! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
    //! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
    void NodeExchangePtr(int rId, NuTo::NodeBase* rOldPtr, NuTo::NodeBase* rNewPtr);
#endif //SWIG

//*************************************************
//************ Element routines     ***************
//**  defined in structures/StructureElement.cpp **
//*************************************************
    //! @brief returns the number of nodes
    //! @return number of nodes
    int GetNumElements() const;

#ifndef SWIG
    //! @brief returns a reference to an element
    //! @param identifier
    //! @return reference to an element
    ElementBase* ElementGetElementPtr(int rIdent);

    //! @brief returns a reference to an element
    //! @param identifier
    //! @return reference to an element
    const ElementBase* ElementGetElementPtr(int rIdent)const;

    //! @brief gives the identifier of an element
    //! @param pointer to an element
    //! @return identifier
    int ElementGetId(const ElementBase* rElement)const;

    //! @brief info about one single element
    //! @param rElement (Input) ... pointer to the element
    //! @param rVerboseLevel (Input) ... level of verbosity
    void ElementInfo(const ElementBase* rElement, int rVerboseLevel)const;
#endif //SWIG

    //! @brief returns a vector with the node ids of an element
    //! @param identifier
    //! @return vector with node ids
    NuTo::FullVector<int,Eigen::Dynamic> ElementGetNodes(int rId);

    //! @brief info about the elements in the Structure
    void ElementInfo(int rVerboseLevel)const;

    //! @brief Creates an element
    //! @param rElementType element type
    //! @param rNodeIdents Identifier for the corresponding nodes
    //! @param rElementDataType data of the element (nonlocal,ip)
    //! @param rIpDataType data of the ip (static data, nonlocal data,..)
    //! @return element number
    int ElementCreate (const std::string& rElementType,
            const NuTo::FullVector<int,Eigen::Dynamic>& rNodeNumbers, const std::string& rElementDataType, const std::string& rIpDataType);

    //! @brief Creates an element
    //! @param rElementType element type
    //! @param rNodeIdents Identifier for the corresponding nodes
    //! @return element number
    int ElementCreate (const std::string& rElementType, const NuTo::FullVector<int,Eigen::Dynamic>& rNodeNumbers);

    //! @brief Creates an element
    //! @param rElementIdent identifier for the element
    //! @param rElementType element type
    //! @param rNodeIdents Identifier for the corresponding nodes
    void ElementCreate (int rElementNumber, const std::string& rElementType,
            const NuTo::FullVector<int,Eigen::Dynamic>& rNodeNumbers);

    //! @brief Creates an element
    //! @param rElementNumber element number
    //! @param rElementType element type
    //! @param rNodeIdents Identifier for the corresponding nodes
    //! @param rElementDataType data of the element (nonlocal,ip)
    //! @param rIpDataType data of the ip (static data, nonlocal data,..)
    //! @return element number
    void ElementCreate (int rElementNumber, const std::string& rElementType,
    		const NuTo::FullVector<int,Eigen::Dynamic>& rNodeNumbers, const std::string& rElementDataType, const std::string& rIpDataType);

    //! @brief creates multiple elements
    //! @param rElementType element type
    //! @param rNodeIdents Identifier for the corresponding nodes (Incidences have to be stored column-wise)
    //! @return a NuTo::FullVector<int,Eigen::Dynamic> containing the element numbers
    NuTo::FullVector<int,Eigen::Dynamic> ElementsCreate (const std::string& rElementType, NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic>& rNodeNumbers);

    //! @brief creates multiple elements
    //! @param rElementType element type
    //! @param rNodeIdents Identifier for the corresponding nodes (Incidences have to be stored column-wise)
    //! @param rElementDataType Element data for the elements
    //! @param rIpDataType Integration point data for the elements
    //! @return a NuTo::FullVector<int,Eigen::Dynamic> containing the element numbers
    NuTo::FullVector<int,Eigen::Dynamic> ElementsCreate (const std::string& rElementType, NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic>& rNodeNumbers, const std::string& rElementDataType, const std::string& rIpDataType);

    //! @param rGroupNumberElements group for elements on the real boundary
    //! @param rGroupNumberBoundaryNodes nodes on the boundary
    //! @param int rOrder number of additional boundary nodes
    //! @param rVirtualBoundary thickness of the virtual boundary element (related to the nonlocal radius)
    //! @param rElementType element type
    //! @param rNodeIdents Identifier for the corresponding nodes
    void BoundaryElementsCreate (const std::string& rElementType,
    		int rGroupNumberElements, int rGroupNumberBoundaryNodes,
    		int rOrder, double rVirtualBoundary,
    		const std::string& rElementDataType, const std::string& rIpDataType);

    //! @param rGroupNumberElements group for elements (Plane2D4N) to be converted
    //! @param rOder, order of the elements (2,3 or 4 is implemented)
    //! @param rNodeDistanceMerge Distance of nodes to be joined (should be significantly smaller than the node distance in the mesh)
    //! @param approximate size of the elements
    void ElementConvertPlane2D4NToPlane2D4NSpectral (int rGroupNumberElements,
    		int rOrder, double rNodeDistanceMerge, double meshSize);

    //! @param rGroupNumberElements group for elements (Plane2D4N) to be converted
    //! @param rElementTypeStr new element type (must be a triangular type)
    //! @param rNodeDistanceMerge distance between nodes for merging
    //! @param rMeshSize for creation of a more efficient neighbor search
    void ElementConvertPlane2D3N (int rGroupNumberElements,
    		std::string rElementTypeStr, double rNodeDistanceMerge, double rMeshSize);

#ifndef SWIG
    //! @brief Create boundary elements defined by all boundary elements and the nodes characterizing the edges
    //! @param rGroupNumberElements group for elements on the real boundary
    //! @param rGroupNumberBoundaryNodes nodes on the boundary
    //! @param int rOrder number of additional boundary nodes
    //! @param rVirtualBoundary thickness of the virtual boundary element (related to the nonlocal radius)
    //! @param rElementType element type
    //! @param rNodeIdents Identifier for the corresponding nodes
    void BoundaryElementsCreate (NuTo::Element::eElementType rType,
    		const Group<ElementBase>* rGroupElements, const Group<NodeBase>* rGroupBoundaryNodes,
    		int rOrder, double rVirtualBoundary,
    		ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType);

#endif //SWIG


    //! @brief Deletes an element
    //! @param rElementNumber element number
    void ElementDelete (int rElementNumber);

    //! @brief Deletes a group of elements element
    //! @param rGroupNumber group number
    void ElementGroupDelete (int rGroupNumber, bool deleteNodes);

#ifndef SWIG
   //! @brief Deletes an element
    //! @param rItElement iterator of the map
    void ElementDeleteInternal(int rElementId);

    //! @brief Creates an element
    //! @param rElementNumber element number
    //! @param rElementType element type
    //! @param rNodeIdents pointers to the corresponding nodes
    //! @return int rElementNumber
    int ElementCreate(Element::eElementType rType, const std::vector<NodeBase*>& rNodeVector,
    		ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType);

    //! @brief Creates an element
    //! @param rElementNumber element number
    //! @param rElementType element type
    //! @param rNodeIdents pointers to the corresponding nodes
    void ElementCreate(int rElementNumber, Element::eElementType rType, const std::vector<NodeBase*>& rNodeVector,
    		ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType);

    //! @brief Returns the internal enum (number) for the element types
    //! @param Element name in Nuto
    //! @return enum
    NuTo::Element::eElementType ElementTypeGetEnum(const std::string& rElementType);

#endif //SWIG

    //! @brief import from gmsh
    //! @param rFileName .. file name
    //! @param rNumTimeDerivatives .. number of time derivatives (0,1, or 2)
    //! @param rDOFs .. degrees of freedom for the nodes
    //! @param rElementData .. element data for the elements to be created
    //! @param rIPData .. ip data for the integration points to be created
    void ImportFromGmsh (const std::string& rFileName,
    		int rNumTimeDerivatives,
    		const std::string& rDOFs, const std::string& rElementData, const std::string& rIPData);

    //! @brief import from gmsh
    //! @param rFileName .. file name
    //! @param rNumTimeDerivatives .. number of time derivatives (0,1, or 2)
    //! @param rDOFs .. degrees of freedom for the nodes
    //! @param rElementData .. element data for the elements to be created
    //! @param rIPData .. ip data for the integration points to be created
    //! @param vector with the created groupes
    void ImportFromGmsh (const std::string& rFileName,
    		int rNumTimeDerivatives,
    		const std::string& rDOFs, const std::string& rElementData, const std::string& rIPData,
    		NuTo::FullVector<int,Eigen::Dynamic>& rElementGroupIds);

    //***********************************************************
    //************         Mesh routines        *****************
    //**  defined in structures/unstructured/StructureMesh.cpp **
    //***********************************************************
    //! @brief creates a lattice mesh from the positions of the circles
    //! @parameters rTypeOfSpecimen 0 box, 1 dogbone
    //! @parameters rBoundingBox box for the spheres (3*2 matrix)
    //! @parameters rCircles (coordinates x,y and radius)
    //! @parameters rTriangles (triangles connecting the circle centers)
   void MeshCreateLattice2D(int rTypeOfSpecimen, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rBoundingBox, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCircles, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rTriangles);

    //! @brief creates a lattice mesh from the positions of the spheres and the bounding box
    //! @parameters rTypeOfSpecimen 0 box, 1 dogbone
    //! @parameters rBoundingBox box for the spheres (3*2 matrix)
    //! @parameters rBoundingBox (min and max for x and y)
    //! @parameters rSpheres (coordinates x,y,z and radius)
    void MeshCreateLattice3D(int rTypeOfSpecimen, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rBoundingBox, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rSpheres, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rTetraeders);

    //*************************************************************
    //************         Crack routines        ******************
    //**  defined in structures/unstructured/StructureCrack.cpp **
    //*************************************************************
    //! @brief returns the number of cracks
    //! @return number of cracks
    unsigned int GetNumCracks() const;

    //! @brief ... Info routine that prints general information about the cracks
    //! @param ... rVerboseLevel describes how detailed the information is
    void CrackInfo(int rVerboseLevel)const;

    //! @brief ... create a new crack
    //! @param rIdent ... crack identifier
    int CrackCreate();

    //! @brief ... create a new crack with given node-Id's
    //! @param rNodes ... vector of node-Id's
    int CrackCreate(NuTo::FullVector<int,Eigen::Dynamic>& rNodes);

    //! @brief ... delete an existing crack
    //! @param rIdent ... crack identifier
    void CrackDelete(int rIdent);

    //! @brief ... extends an existing crack
    //! @param rIdent (Input) ... crack identifier
    //! @param rNode (Input) ... node Id to be attended to the crack
    void CrackPushBack(const int rIdent, const int rNodeNumber);
    void CrackPushFront(const int rIdent, const int rNodeNumber);

    //! @brief ... shortens an existing crack
    //! @param rIdent ... crack identifier
    void CrackPopBack(int rIdent);
    void CrackPopFront(int rIdent);

    //! @brief ... merge all cracks to the existing structure
    void InitiateCracks();
    //! @brief ... merge specified crack to the existing structure
    //! @param rIdent ... crack identifier
    void InitiateCrack(const int rIdent);

    //! @brief ... Initialize the PhantomNodeMethod to the structure
    void InitiatePhantomNodeMethod();

    //! @brief ... take cracked elements and initiate PhantomNodeMethod
    //! @param rNumIp (Input) ... number of integration points for the new (cracked) elements
    //! @return  ... id vector of cracked elements
    NuTo::FullVector<int,Eigen::Dynamic> InitiatePhantomNodeMethod(int rNumIp);

    //! @brief ... take cracked elements and initiate PhantomNodeMethod
    //! @param rNumIp (Input) ... number of integration points for the new (cracked) elements
    //! @return  ... id vector of cracked elements
    NuTo::FullVector<int,Eigen::Dynamic> InitiatePhantomNodeMethodTriangle(int rNumIp);
#ifndef SWIG
    //! @brief ... extends an existing crack
    //! @param rIdent ... crack identifier
    //! @param rNode ... pointer to the node to be attended to the crack
    void CrackPushBack(int rIdent, NodeBase* rNode);
    void CrackPushFront(int rIdent, NodeBase* rNode);

    //! @brief ... merge all cracks to the existing structure
    //! @param rCrackedElems (Output) ... vector of cracked elements
    void InitiateCracks(elementBasePtrSet_t & rCrackedElems);

    //! @brief ... merge specified crack to the existing structure
    //! @param rIdent (Input) ... crack identifier
    //! @param rCrackedElems (Output) ... vector of cracked elements
    void InitiateCrack(const int rIdent, elementBasePtrSet_t & rCrackedElems);

    //! @brief ... Initialize the PhantomNodeMethod to the structure: rCrackedElems will be doubled after initialization
    //! @param rCrackedElems (Input) ... vector of cracked elements
    void InitiatePhantomNodeMethod(elementBasePtrSet_t & rCrackedElems);

    //! @brief returns a reference to a crack
    //! @param identifier
    //! @return reference to a crack
    CrackBase* CrackGetCrackPtr(int rIdent);

    //! @brief returns a const reference to a crack
    //! @param identifier
    //! @return const reference to a crack
    const CrackBase* CrackGetCrackPtr(int rIdent)const;

    //! @brief gives the identifier of a crack
    //! @param pointer to a crack
    //! @return identifier
    int CrackGetId(const CrackBase* rCrack)const;

    //! @brief copy and move the structure
    //! most of the data is kept, but e.g. nonlocal data and
    //! @param rOffset offset (dimension x 1 has to be identical with structure dimension)
    //! @param rOld2NewNodePointer ptrMap showing the new and old node pointers
    //! @param rOld2NewElementPointer ptrMap showing the new and old element pointers
    void CopyAndTranslate(NuTo::FullVector<double,Eigen::Dynamic>& rOffset, std::map<NodeBase*, NodeBase* >& rOld2NewNodePointer, std::map<ElementBase*, ElementBase* >& rOld2NewElementPointer);
#endif //SWIG

    //! @brief copy and move the structure
    //! most of the data is kept, but e.g. nonlocal data and
    //! @param rOffset offset (dimension x 1 has to be identical with structure dimension)
    void CopyAndTranslate(NuTo::FullVector<double,Eigen::Dynamic>& rOffset);

    //! @brief ... Adds an element to an element group
    //! @param ... rIdentGroup identifier for the group
    //! @param ... rIdentElement  identifier for the element
    void GroupAddElement(int rIdentGroup, int rIdElement);

    //! @brief ... Adds all elements to a group based on the type
    //! @param ... rIdentGroup identifier for the group
    //! @param ... rElemTypeStr  identifier for the element type
    void GroupAddElementFromType(int rIdentGroup, std::string rElemTypeStr);


    //*************************************************
    //************ Info routine         ***************
    //**  defined in structures/Structure.cpp *********
    //*************************************************
    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info()const;

protected:
    //! @brief ... standard constructor just for the serialization routine
    Structure()
    {}

#ifndef SWIG
    //! @brief ... store all elements of a structure in a vector
    //! @param rElements ... vector of element pointer
    void GetElementsTotal(std::vector<const ElementBase*>& rElements) const;

    //! @brief ... store all elements of a structure in a vector
    //! @param rElements ... vector of element pointer
    void GetElementsTotal(std::vector<std::pair<int,const ElementBase*> >& rElements) const;

    //! @brief ... store all elements of a structure in a vector
    //! @param rElements ... vector of element pointer
    void GetElementsTotal(std::vector<ElementBase*>& rElements);

    //! @brief ... store all elements of a structure in a vector
    //! @param rElements ... vector of element pointer
    void GetElementsTotal(std::vector<std::pair<int, ElementBase*> >& rElements);

    //! @brief ... store all nodes of a structure in a vector
    //! @param rNodes ... vector of element pointer
    void GetNodesTotal(std::vector<const NodeBase*>& rNodess) const;

    //! @brief ... store all nodes of a structure in a vector
    //! @param rNodes ... vector of element pointer
    void GetNodesTotal(std::vector<std::pair<int,const NodeBase*> >& rNodes) const;

    //! @brief ... store all nodes of a structure in a vector
    //! @param rNodes ... vector of element pointer
    void GetNodesTotal(std::vector<NodeBase*>& rNodes);

    //! @brief ... store all nodes of a structure in a vector
    //! @param rNodes ... vector of element pointer
    void GetNodesTotal(std::vector<std::pair<int, NodeBase*> >& rNodes);

    //! @brief ... store all cracks of a structure in a vector
    //! @param rElements ... vector of const crack pointer
    void GetCracksTotal(std::vector<const CrackBase*>& rCracks) const;

    //! @brief ... store all cracks of a structure in a vector
    //! @param rElements ... vector of crack pointer
    void GetCracksTotal(std::vector<CrackBase*>& rCracks);
#endif

    //! @brief deletes a node
    //! @param rNodeNumber ... node number
    //! @param checkElements ... check the elements, if set to false, make sure that the node is not part of any element
    void NodeDelete(int rNodeNumber, bool checkElements);

#ifndef SWIG
    //! @brief import from gmsh, do the actual work
    //! @param rFileName .. file name
    //! @param rNumTimeDerivatives number of time derivatives (0,1 or 2)
    //! @param rDOFs .. degrees of freedom for the nodes
    //! @param rElementData .. element data for the elements to be created
    //! @param rIPData .. ip data for the integration points to be created
    //! @param vector with the created groupes
    void ImportFromGmshAux (const std::string& rFileName,
    		int rNumTimeDerivatives,
    		const std::string& rDOFs, const std::string& rElementData, const std::string& rIPData,
    		bool rAddGroups, std::set<int>& rElementGroupIds);
#endif //SWIG

    boost::ptr_map<int,NodeBase> mNodeMap;
    boost::ptr_map<int,ElementBase> mElementMap;


    //! @brief ... map storing the cracks and a pointer to the objects
    //! @sa CrackBase
    crackMap_t mCrackMap;


};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::Structure)
#endif // SWIG
#endif // ENABLE_SERIALIZATION


#endif // STRUCTURE_H
