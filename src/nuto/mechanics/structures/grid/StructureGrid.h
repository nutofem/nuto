// $Id$
#ifndef STRUCTUREGRID_H
#define STRUCTUREGRID_H
#include "nuto/base/NuToObject.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/math/FullMatrix.h"
#include <boost/dynamic_bitset.hpp>

namespace NuTo
{

//! @brief forward declaration to speed up compilation time
class NodeGrid3D;
class NodeGridDisplacements3D;
class Voxel8N;

//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... regular structure e.g. from pixel/voxel data
class StructureGrid: public NuToObject
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    //! @param mDimension  Structural dimension (1,2 or 3)
    StructureGrid(int rDimension);

    ~StructureGrid();

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

    //! @brief import routine for basic grid data without StructureGrid data space
    void ImportFromVtkASCIIFileHeader(const char* rFileName,int *rGridDimension,double *rVoxelSpacing,double *rGridOrigin, int &rNumVoxel);

    //! @brief import routine for basic grid data with StructureGrid data space
    void ImportFromVtkASCIIFileHeader(const char* rFileName);

	//! @brief ... imports Data from a Vtk ASCII File
	//! @param fileName ... file name
	void ImportFromVtkASCIIFile(const char* fileName,std::vector<int> &rData);

   //! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    std::string GetTypeId()const
    {
        return std::string("StructureGrid");
    }

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
    FullMatrix<double>* GetLocalCoefficientMatrix0(int rNumLocalCoefficientMatrix0);

    //! @brief Get VoxeLNumAndLocMatrix
    //! @return FullMatrix columns elements, rows voxel number and number in x, y,z direction
    FullMatrix<int>* GetVoxelNumAndLocMatrix();

   //! @brief Calculate ElementVoxelLocMatrix
   void CalculateVoxelLocations();

    //! @brief Get voxels corner numbers from bottom to top counter-clockwise
    //! @return array of number of corners with corner numbers
    void GetCornersOfVoxel(int rElementNumber,int *rVoxLoc,int *corners);

    //! @brief Set NodeIds for all nodes at all elements
    void SetAllNodeIds();

    //! @brief Set NodeIds for all nodes at all nodes
    void SetAllNodeIdsAtNode();

    //! @brief Set ElementIds for the elements at each nodes
    void SetAllElementIds();

    //! @brief GetConstraintSwitch
    //! @return switch field for constraint
    bool* GetConstraintSwitch();

    //! @brief create node data without StructureGrid
    void CreateGrid(int rThresholdMaterialValue, std::vector<int> &imageValues ,const std::vector<double>& rColorToMaterialData,int* rGridDimension,boost::dynamic_bitset<> &rNodeExist,boost::dynamic_bitset<> &rElemExist,std::vector<double>& youngsModulus,std::vector<int>& materialOfElem);


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
    NodeGrid3D* NodeGridGetNodePtr(int rIdent);

    //! @brief returns a reference to a node
    //! @param identifier
    //! @return reference to a node
    const NodeGrid3D* NodeGridGetNodePtr(int rIdent)const;

    //! @brief gives the identifier of a node
    //! @param reference to a node
    //! @return identifier
    int NodeGetId(const NodeGrid3D* rNode)const;

    //! @brief gets the displacements of a node
     //! @param rNode node identifier
     //! @param rDisplacements matrix (one column) with the displacements
     void NodeGetDisplacements(int rNode, NuTo::FullMatrix<double>& rDisplacements)const;

     //! @brief sets the displacements of a node
     //! @param rIdent node identifier
     //! @param rDisplacements matrix (one column) with the displacements
     void NodeSetDisplacements(int rId,const NuTo::FullMatrix<double>& rDisplacements);


#endif //SWIG
    //! @param rNodeNumber ... node number
    void NodeDelete(const int rNodeNumber);

    //! @brief info about the nodes in the Structure
    virtual void NodeInfo(int mVerboseLevel) const;

    void NodeCreate(bool flag, int rNodeID, std::string rDOFs);

    void CreateNodeGrid(std::string rDOFs,int rThresholdMaterialValue);

    int* GetCoincidenceVoxelIDs(int rNodeID);

    //! @brief NodeSetConstraintSwitch
    //! @param rGridNodeNum, rDirection, rConstraint
    void NodeSetConstraintSwitch(int rGridNodeNum, int rDirection, bool rConstraint);

    //! @brief NodeGetConstraintSwitch
    //! @param rGlobalDof
    //! @return switch for constraint
    bool NodeGetConstraintSwitch(int rGlobalDof);

    //! @brief Set partCoefficientmatrix for all nodes
    void SetAllPartCoefficientMatrix0();

/*
    //! @brief numbers the dofs in the structure
    void NodeBuildGlobalDofs();

    //! @brief extract dof values (e.g. displacements, temperatures to the nodes)
    //! @param rActiveDofValues ... vector of global active dof values (ordering according to global dofs, size is number of active dofs)
    //! @param rDependentDofValues ... vector of global dependent dof values (ordering according to (global dofs) - (number of active dofs), size is (total number of dofs) - (number of active dofs))
    void NodeExtractDofValues(NuTo::FullMatrix<double>& rActiveDofValues, NuTo::FullMatrix<double>& rDependentDofValues) const;

    //! @brief merge dof values
    void NodeMergeActiveDofValues(const NuTo::FullMatrix<double>& rActiveDofValues);
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
    Voxel8N* ElementGetPtr(int rIdent);

    //! @brief returns a reference to an element
    //! @param identifier
    //! @return reference to an element
    const Voxel8N* ElementGetPtr(int rIdent) const;

    //! @brief gives the identifier of an element
    //! @param reference to an element
    //! @return identifier
    int ElementGetId(Voxel8N* rElement);
#endif //SWIG

    //! @brief info about the elements in the Structure
    virtual void ElementInfo(int mVerboseLevel) const;

    void CreateElementGrid(NuTo::FullMatrix<double>& rBaseCoefficientMatrix0,
            const NuTo::FullMatrix<double>& rColorToMaterialData,const std::string& rElementType);

    //! @brief Creates an element
    //! @param number of local coefficient matrix
    //! @param rElementID identifier for the element
     void ElementCreate (bool flag, int rNumCoefficientMatrix0, int rElementID);

    //! @brief Creates an element
     //! @param number of local coefficient matrix
    //! @param rElementID identifier for the element
    //! @param rElementType element type
    void ElementCreate (bool flag, int rNumCoefficientMatrix0,int rElementID,const std::string& rElementType);

    //! @brief Deletes an element
    //! @param rElementNumber element number
    void ElementDelete (int rElementNumber);

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info()const;


protected:
    int mDimension;
    int mNumVoxel;  //number of voxels
    //! @todo length of list in function of real dimension
    double mVoxelSpacing[3]; //spacing between center of neighbor voxels / dimension of each voxel
    int mGridDimension[3]; //dimension of the voxel model
    double mGridOrigin[3];// origin of the model , in the center of the first voxel
    std::vector<NodeGrid3D*> mNodeVec;
    std::vector<Voxel8N*> mElementVec;
    //boost::ptr_vector<ElementBase> mElementVec;
    const char* mImageDataFile;
    int mNumMaterials;
    typedef NuTo::FullMatrix<double> FullMat;
    std::vector<FullMat > mLocalCoefficientMatrix0;
    NuTo::FullMatrix<int>* mVoxelLocation;
    bool* mDofIsNotConstraint; //field of bool for all dofs, length 3xnumGridNodes, constraint dof = 0 = false
    bool mCalcVoxelLocation; //true != 0 when already calculated
};
} //namespace NuTo
#endif // STRUCTUREGRID_H
