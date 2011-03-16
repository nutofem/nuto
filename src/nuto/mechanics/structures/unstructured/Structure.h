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

    //! @brief ... store all elements connected to this node in a vector
    //! @param rNodeId (Input) 			... node id
    //! @param rElementNumbers (Output)	... vector of element ids
    void NodeGetElements(const int rNodeId, NuTo::FullMatrix<int>& rElementNumbers);

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
    int NodeCreate(std::string rDOFs, NuTo::FullMatrix<double>& rCoordinates);

    //! @brief creates a node
    //! @param rNodeNumber ... node number
    //! @param rDOFs ... space separated string containing the node dofs (e.g. displacements, rotations, temperatures)
    //! @param rCoordinates ...  node coordinates
    void NodeCreate(int rNodeNumber, std::string rDOFs, NuTo::FullMatrix<double>& rCoordinates);

    //! @brief creates multiple nodes
    //! @param rDOFs ... space separated string containing the nodal dofs (e.g. displacements, rotations, temperatures)
    //! @param rCoordinates ...  nodal coordinates (column-wise storage of each nodal coordinate)
    //! @return a NuTo::FullMatrix<int> containing the node numbers
    NuTo::FullMatrix<int> NodesCreate(std::string rDOFs, NuTo::FullMatrix<double>& rCoordinates);

    //! @brief deletes a node
    //! @param rNodeNumber ... node number
    void NodeDelete(const int rNodeNumber);

    //! @brief info about the nodes in the Structure
    void NodeInfo(int mVerboseLevel)const;

    //! @brief numbers the dofs in the structure
    void NodeBuildGlobalDofs();

    //! @brief renumbers the global dofs in the structure after
    //! @only relevant for structureip with global dofs, otherwise just empty
    virtual void ReNumberAdditionalGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
    {}

    //! @brief numbers non standard DOFs' e.g. in StructureIp, for standard structures this routine is empty
    virtual void NumberAdditionalGlobalDofs(){};

    //! @brief write dof values (e.g. displacements, temperatures to the nodes)
    //! @param rActiveDofValues ... vector of global dof values (ordering according to global dofs, size is number of active dofs)
    virtual void NodeMergeActiveDofValues(const NuTo::FullMatrix<double>& rActiveDofValues);

    //! @brief ...merge additional dof values
    //! empty in the standard case, is used in StructureIp
    virtual void NodeMergeAdditionalGlobalDofValues(const NuTo::FullMatrix<double>& rActiveDofValues, const NuTo::FullMatrix<double>& rDependentDofValues){};

    //! @brief extract dof values (e.g. displacements, temperatures to the nodes)
    //! @param rActiveDofValues ... vector of global active dof values (ordering according to global dofs, size is number of active dofs)
    //! @param rDependentDofValues ... vector of global dependent dof values (ordering according to (global dofs) - (number of active dofs), size is (total number of dofs) - (number of active dofs))
    virtual void NodeExtractDofValues(NuTo::FullMatrix<double>& rActiveDofValues, NuTo::FullMatrix<double>& rDependentDofValues) const;

    //! @brief extract dof values additional dof values
    //! @param rActiveDofValues ... vector of global active dof values (ordering according to global dofs, size is number of active dofs)
    //! @param rDependentDofValues ... vector of global dependent dof values (ordering according to (global dofs) - (number of active dofs), size is (total number of dofs) - (number of active dofs))
    virtual void NodeExtractAdditionalGlobalDofValues(NuTo::FullMatrix<double>& rActiveDofValues, NuTo::FullMatrix<double>& rDependentDofValues) const{};

    //! @brief write first time derivative of dof values (e.g. velocities) to the nodes
    //! @param rActiveDofValues ... global vector of the first time derivatives of dof values (ordering according to global dofs, size is number of active dofs)
    void NodeMergeActiveDofFirstTimeDerivativeValues(const NuTo::FullMatrix<double>& rActiveDofValues);

    //! @brief extract first time derivatives of dof values (e.g. velocities) from the nodes
    //! @param rActiveDofValues ... vector of global active dof values (ordering according to global dofs, size is number of active dofs)
    //! @param rDependentDofValues ... vector of global dependent dof values (ordering according to (global dofs) - (number of active dofs), size is (total number of dofs) - (number of active dofs))
    void NodeExtractDofFirstTimeDerivativeValues(NuTo::FullMatrix<double>& rActiveDofValues, NuTo::FullMatrix<double>& rDependentDofValues) const;

    //! @brief write second time derivative of dof values (e.g. accelerations) to the nodes
    //! @param rActiveDofValues ... global vector of the first time derivatives of dof values (ordering according to global dofs, size is number of active dofs)
    void NodeMergeActiveDofSecondTimeDerivativeValues(const NuTo::FullMatrix<double>& rActiveDofValues);

    //! @brief extract second time derivatives of dof values (e.g. accelerations) from the nodes
    //! @param rActiveDofValues ... vector of global active dof values (ordering according to global dofs, size is number of active dofs)
    //! @param rDependentDofValues ... vector of global dependent dof values (ordering according to (global dofs) - (number of active dofs), size is (total number of dofs) - (number of active dofs))
    void NodeExtractDofSecondTimeDerivativeValues(NuTo::FullMatrix<double>& rActiveDofValues, NuTo::FullMatrix<double>& rDependentDofValues) const;

    //! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
    //! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
    void NodeExchangePtr(NuTo::NodeBase* rOldPtr, NuTo::NodeBase* rNewPtr);

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
#endif //SWIG

    //! @brief info about the elements in the Structure
    void ElementInfo(int mVerboseLevel)const;

    //! @brief Creates an element
    //! @param rElementType element type
    //! @param rNodeIdents Identifier for the corresponding nodes
    //! @param rElementDataType data of the element (nonlocal,ip)
    //! @param rIpDataType data of the ip (static data, nonlocal data,..)
    //! @return element number
    int ElementCreate (const std::string& rElementType,
            const NuTo::FullMatrix<int>& rNodeNumbers, const std::string& rElementDataType, const std::string& rIpDataType);

    //! @brief Creates an element
    //! @param rElementType element type
    //! @param rNodeIdents Identifier for the corresponding nodes
    //! @return element number
    int ElementCreate (const std::string& rElementType, const NuTo::FullMatrix<int> &rNodeNumbers);

    //! @brief Creates an element
    //! @param rElementIdent identifier for the element
    //! @param rElementType element type
    //! @param rNodeIdents Identifier for the corresponding nodes
    void ElementCreate (int rElementNumber, const std::string& rElementType,
            const NuTo::FullMatrix<int> &rNodeNumbers);

    //! @brief Creates an element
    //! @param rElementNumber element number
    //! @param rElementType element type
    //! @param rNodeIdents Identifier for the corresponding nodes
    //! @param rElementDataType data of the element (nonlocal,ip)
    //! @param rIpDataType data of the ip (static data, nonlocal data,..)
    //! @return element number
    void ElementCreate (int rElementNumber, const std::string& rElementType,
    		const NuTo::FullMatrix<int> &rNodeNumbers, const std::string& rElementDataType, const std::string& rIpDataType);

    //! @brief creates multiple elements
    //! @param rElementType element type
    //! @param rNodeIdents Identifier for the corresponding nodes (Incidences have to be stored column-wise)
    //! @return a NuTo::FullMatrix<int> containing the element numbers
    NuTo::FullMatrix<int> ElementsCreate (const std::string& rElementType, NuTo::FullMatrix<int> &rNodeNumbers);

    //! @brief creates multiple elements
    //! @param rElementType element type
    //! @param rNodeIdents Identifier for the corresponding nodes (Incidences have to be stored column-wise)
    //! @param rElementDataType Element data for the elements
    //! @param rIpDataType Integration point data for the elements
    //! @return a NuTo::FullMatrix<int> containing the element numbers
    NuTo::FullMatrix<int> ElementsCreate (const std::string& rElementType, NuTo::FullMatrix<int> & rNodeNumbers, const std::string& rElementDataType, const std::string& rIpDataType);

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
    int ElementCreate(Element::eElementType rType, std::vector<NodeBase*> rNodeVector,
    		ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType);

    //! @brief Creates an element
    //! @param rElementNumber element number
    //! @param rElementType element type
    //! @param rNodeIdents pointers to the corresponding nodes
    void ElementCreate(int rElementNumber, Element::eElementType rType, std::vector<NodeBase*> rNodeVector,
    		ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType);

#endif //SWIG

    //! @brief import from gmsh
    //! @param rFileName .. file name
    //! @param rDOFs .. degrees of freedom for the nodes
    //! @param rElementData .. element data for the elements to be created
    //! @param rIPData .. ip data for the integration points to be created
    void ImportFromGmsh (const std::string& rFileName,
    		const std::string& rDOFs, const std::string& rElementData, const std::string& rIPData);

    //! @brief import from gmsh
    //! @param rFileName .. file name
    //! @param rDOFs .. degrees of freedom for the nodes
    //! @param rElementData .. element data for the elements to be created
    //! @param rIPData .. ip data for the integration points to be created
    //! @param vector with the created groupes
    void ImportFromGmsh (const std::string& rFileName,
    		const std::string& rDOFs, const std::string& rElementData, const std::string& rIPData,
    		NuTo::FullMatrix<int>& rElementGroupIds);

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

    //! @brief ... delete an existing crack
    //! @param rIdent ... crack identifier
    void CrackDelete(int rIdent);

    //! @brief ... extends an existing crack
    //! @param rIdent ... crack identifier
    //! @param rNode ... pointer to the node to be attended to the crack
    void CrackPushBack(int rIdent, NodeBase* rNode);
    void CrackPushFront(int rIdent, NodeBase* rNode);

    //! @brief ... shortens an existing crack
    //! @param rIdent ... crack identifier
    void CrackPopBack(int rIdent);
    void CrackPopFront(int rIdent);

    //! @brief ... merge all cracks to the existing structure
    void InitiateCracks();
    //! @brief ... merge all cracks to the existing structure
    //! @param rCrackedElems (Output) ... vector of cracked elements
    void InitiateCracks(elementBasePtrSet_t & rCrackedElems);

    //! @brief ... merge specified crack to the existing structure
    //! @param rIdent ... crack identifier
    void InitiateCrack(const int rIdent);
    //! @brief ... merge specified crack to the existing structure
    //! @param rIdent (Input) ... crack identifier
    //! @param rCrackedElems (Output) ... vector of cracked elements
    void InitiateCrack(const int rIdent, elementBasePtrSet_t & rCrackedElems);

    //! @brief ... Initialize the PhantomNodeMethod to the structure: rCrackedElems will be doubled after initialization
    //! @param rCrackedElems (Input) ... vector of cracked elements
    void InitiatePhantomNodeMethod(elementBasePtrSet_t & rCrackedElems);

#ifndef SWIG
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
#endif //SWIG

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

    //! @brief ... store all elements of a structure in a vector
    //! @param rElements ... vector of element pointer
    void GetElementsTotal(std::vector<const ElementBase*>& rElements) const;

    //! @brief ... store all elements of a structure in a vector
    //! @param rElements ... vector of element pointer
    void GetElementsTotal(std::vector<ElementBase*>& rElements);

    //! @brief ... store all nodes of a structure in a vector
    //! @param rNodes ... vector of element pointer
    void GetNodesTotal(std::vector<const NodeBase*>& rNodess) const;

    //! @brief ... store all nodes of a structure in a vector
    //! @param rNodes ... vector of element pointer
    void GetNodesTotal(std::vector<NodeBase*>& rNodes);

    //! @brief ... store all cracks of a structure in a vector
    //! @param rElements ... vector of const crack pointer
    void GetCracksTotal(std::vector<const CrackBase*>& rCracks) const;

    //! @brief ... store all cracks of a structure in a vector
    //! @param rElements ... vector of crack pointer
    void GetCracksTotal(std::vector<CrackBase*>& rCracks);

    //! @brief deletes a node
    //! @param rNodeNumber ... node number
    //! @param checkElements ... check the elements, if set to false, make sure that the node is not part of any element
    void NodeDelete(int rNodeNumber, bool checkElements);

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    virtual void BuildGlobalCoefficientSubMatrices0General(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK) const;

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    //! @param rMatrixKJ ... submatrix kj (number of dependent dof x number of active dof)
    //! @param rMatrixKK ... submatrix kk (number of dependent dof x number of dependent dof)
    virtual void BuildGlobalCoefficientSubMatrices0General(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK, NuTo::SparseMatrix<double>& rMatrixKJ, NuTo::SparseMatrix<double>& rMatrixKK) const;

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    virtual void BuildGlobalCoefficientSubMatrices0Symmetric(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK) const;

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    //! @param rMatrixKK ... submatrix kk (number of dependent dof x number of dependent dof)
    virtual void BuildGlobalCoefficientSubMatrices0Symmetric(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK, NuTo::SparseMatrix<double>& rMatrixKK) const;

    //! @brief ... based on the global dofs build sub-vectors of the global internal potential gradient
    //! @param rActiveDofGradientVector ... global internal potential gradient which corresponds to the active dofs
    //! @param rDependentDofGradientVector ... global internal potential gradient which corresponds to the dependent dofs
    virtual void BuildGlobalGradientInternalPotentialSubVectors(NuTo::FullMatrix<double>& rActiveDofGradientVector, NuTo::FullMatrix<double>& rDependentDofGradientVector) const;

#ifndef SWIG
    //! @brief import from gmsh, do the actual work
    //! @param rFileName .. file name
    //! @param rDOFs .. degrees of freedom for the nodes
    //! @param rElementData .. element data for the elements to be created
    //! @param rIPData .. ip data for the integration points to be created
    //! @param vector with the created groupes
    void ImportFromGmshAux (const std::string& rFileName,
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
