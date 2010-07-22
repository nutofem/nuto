#ifndef STRUCTUREGRID_H
#define STRUCTUREGRID_H

#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/elements/ElementDataEnum.h"
#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/MechanicsException.h"

#ifdef ENABLE_SERIALIZATION
#include <boost/ptr_container/serialize_ptr_vector.hpp>
#else
#include <boost/ptr_container/ptr_vector.hpp>
#endif //ENABLE_SERIALIZATION

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... regular structure e.g. from pixel/voxel data
class StructureGrid :  public StructureBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    //! @param mDimension  Structural dimension (1,2 or 3)
    StructureGrid(int rDimension);

    typedef NuTo::SparseMatrixCSRGeneral<double> SparseMat ;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

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
        return std::string("StructureGrid");
    }

    virtual void ImportFromVtkASCIIFileHeader(const char* rFileName);

    //! @brief returns number of Voxels
    //! @return number of Voxels
    int GetNumVoxels() const;

    //! @brief returns  VoxelSpacing
    //! @return VoxelSpacing
    const double* GetVoxelSpacing() const;

    //! @brief returns GridOrigin
    //! @return GridOrigin
    const double* GetGridOrigin() const;

     //! @brief returns GridDimension
     //! @return GridDimension
    const int* GetGridDimension() const;

     //! @brief Get NumMaterials
     //! @return NumMaterial
     const int GetNumMaterials() const;


    //! @brief Get LocalCoefficientMatrix0
    //! @param NumLocalCoefficientMatrix0 number of stiffness matrix
    const SparseMat* GetLocalCoefficientMatrix0(int rNumLocalCoefficientMatrix0) const;

//*************************************************
//************ Node routines        ***************
//***  defined in structures/StructureGridNode.cpp  ***
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

    //! @brief a reference to a node
    //! @param node GridNum
    //! @return reference to a node
    NodeBase* NodeGetNodePtrFromGridNum(int rNodeGridNum);

    //! @brief a reference to a node
    //! @param node GridNum
    //! @return reference to a node
    const NodeBase* NodeGetNodePtrFromGridNum(int rNodeGridNum) const;

    //! @brief a reference to a node
    //! @param node ID
    //! @return reference to a node
    int NodeGetNodeNumberFromId(int rNodeId);

    //! @brief a reference to a node
    //! @param node ID
    //! @return reference to a node
    const int NodeGetNodeNumberFromId(int rNodeId) const;



    //! @brief gives the identifier of a node
    //! @param reference to a node
    //! @return identifier
    int NodeGetId(const NodeBase* rNode)const;

    //! @brief a identifier of a node
    //! @param node GridNum
    //! @return identifier
    const int NodeGetIdFromGridNum(int rNodeGridNum) const;

#endif //SWIG
    //! @param rNodeNumber ... node number
    void NodeDelete(const int rNodeNumber);

    //! @brief info about the nodes in the Structure
    virtual void NodeInfo(int mVerboseLevel) const;

    void NodeCreate(int rNodeNumber, int rNodeID, std::string rDOFs);

    void CreateNodeGrid(std::string rDOFs);

    typedef std::vector<int> TCoincidentVoxelList;
    TCoincidentVoxelList GetCoincidenceVoxelIDs(int rNodeID);
    //! @brief numbers the dofs in the structure
    void NodeBuildGlobalDofs();


    //! @brief extract dof values (e.g. displacements, temperatures to the nodes)
    //! @param rActiveDofValues ... vector of global active dof values (ordering according to global dofs, size is number of active dofs)
    //! @param rDependentDofValues ... vector of global dependent dof values (ordering according to (global dofs) - (number of active dofs), size is (total number of dofs) - (number of active dofs))
    void NodeExtractDofValues(NuTo::FullMatrix<double>& rActiveDofValues, NuTo::FullMatrix<double>& rDependentDofValues) const;

    //! @brief merge dof values
    void NodeMergeActiveDofValues(const NuTo::FullMatrix<double>& rActiveDofValues);
/*
    //! @brief get internal forces
    void NodeGetInternalForce(const NodeBase* rNode, NuTo::FullMatrix<double>& rNodeForce)const;

    //! @brief calculates the internal force vector for a given node
    //! @param rNodeId node id
    //! @param rNodeForce return value
    void NodeGetInternalForce(int rNodeId, NuTo::FullMatrix<double>& rNodeForce)const;
*/

//*************************************************
//************ Element routines     ***************
//**  defined in structures/StructureGridElement.cpp **
//***********************************************"elementVec",**
    //! @brief returns the number of elements
    //! @return number of elements
    int GetNumElements() const;

#ifndef SWIG
    //! @brief returns a reference to an element
    //! @param identifier
    //! @return reference to an element
    ElementBase* ElementGetElementPtr(int rIdent);

    //! @brief returns a reference to an element
    //! @param identifier
    //! @return reference to an element
    const ElementBase* ElementGetElementPtr(int rIdent) const;

    //! @brief gives the identifier of an element
    //! @param reference to an element
    //! @return identifier
    int ElementGetId(const ElementBase* rElement) const;
#endif //SWIG

    //! @brief info about the elements in the Structure
    virtual void ElementInfo(int mVerboseLevel) const;

    void CreateElementGrid(NuTo::SparseMatrixCSRGeneral<double>& rBaseCoefficientMatrix0,
            const NuTo::FullMatrix<double>& rColorToMaterialData,const std::string& rElementType);

    //! @brief Creates an element
    //! @param rElementID identifier for the element
    //! @param rElementNumber number of the element
    //! @param rElementType element type
    void ElementCreate (int rNumCoefficientMatrix0,int rElementNumber, int rElementID,
    		const std::string& rElementType);

    //! @brief Creates an element
    //! @param rElementID identifier for the element
    //! @param rElementNumber number of the element
    //! @param rElementType element type
    //! @param rIpDataType ip type
    void ElementCreate (int rNumCoefficientMatrix0,int rElementNumber, int rElementID,
    		const std::string& rElementType, const std::string& rElementDataType, const std::string& rIpDataType);

#ifndef SWIG
    //! @brief Creates an element
    //! @param rElementID identifier for the element
    //! @param rElementNumber number of the element
    //! @param rElementType element type
    //! @param rIpDataType ip type
    void ElementCreate (int rNumCoefficientMatrix0,int rElementNumber, int rElementID,
    		Element::eElementType rElementType, NuTo::ElementData::eElementDataType rElementDataType, NuTo::IpData::eIpDataType rIpDataType);
#endif //SWIG

    //! @brief Deletes an element
    //! @param rElementNumber element number
    void ElementDelete (const int rElementNumber);



protected:
    int mNumVoxel;  //number of voxels
//! @TODO length of list in function of real dimension
    double mVoxelSpacing[3]; //spacing between center of neighbor voxels / dimension of each voxel
    int mGridDimension[3]; //dimension of the voxel model
    double mGridOrigin[3];// origin of the model , in the center of the first voxel
    boost::ptr_vector<NodeBase> mNodeVec;
    boost::ptr_vector<ElementBase> mElementVec;
    const char* mImageDataFile;
    int mNumMaterials;
    std::vector<SparseMat> mLocalCoefficientMatrix0;

    //! @brief ... store all elements of a structure in a vector
    //! @param rElements ... vector of element pointer
    void GetElementsTotal(std::vector<ElementBase*>& rElements);

    //! @brief ... store all elements of a structure in a vector
    //! @param rElements ... vector of element pointer
    void GetElementsTotal(std::vector<const ElementBase*>& rElements) const;

    //! @brief ... store all nodes of a structure in a vector
    //! @param rNodes ... vector of element pointer
    void GetNodesTotal(std::vector<const NodeBase*>& rNodess) const;

    //! @brief ... store all nodes of a structure in a vector
    //! @param rNodes ... vector of element pointer
    void GetNodesTotal(std::vector<NodeBase*>& rNodes);

    //! @brief ... create local coefficient matrix 0 for a voxel and save a pointer to the matrix
    void BuildLocalCoefficientMatrix0() const;

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    void BuildGlobalCoefficientSubMatrices0General(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK) const;

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    //! @param rMatrixKJ ... submatrix kj (number of dependent dof x number of active dof)
    //! @param rMatrixKK ... submatrix kk (number of dependent dof x number of dependent dof)
    void BuildGlobalCoefficientSubMatrices0General(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK, NuTo::SparseMatrix<double>& rMatrixKJ, NuTo::SparseMatrix<double>& rMatrixKK) const;

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    void BuildGlobalCoefficientSubMatrices0Symmetric(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK) const;

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    //! @param rMatrixKK ... submatrix kk (number of dependent dof x number of dependent dof)
    void BuildGlobalCoefficientSubMatrices0Symmetric(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK, NuTo::SparseMatrix<double>& rMatrixKK) const;

    //! @brief ... based on the global dofs build sub-vectors of the global internal potential gradient
    //! @param rActiveDofGradientVector ... global internal potential gradient which corresponds to the active dofs
    //! @param rDependentDofGradientVector ... global internal potential gradient which corresponds to the dependent dofs
    void BuildGlobalGradientInternalPotentialSubVectors(NuTo::FullMatrix<double>& rActiveDofGradientVector, NuTo::FullMatrix<double>& rDependentDofGradientVector) const;

};
} //namespace NuTo
#endif // STRUCTUREGRID_H
