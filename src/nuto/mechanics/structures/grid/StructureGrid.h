// $Id$
#ifndef STRUCTUREGRID_H
#define STRUCTUREGRID_H
//#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/optimize/CallbackHandlerGrid.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/math/FullMatrix.h"
#include <boost/dynamic_bitset.hpp>

namespace NuTo
{

//! @brief forward declaration to speed up compilation time
//class NodeGrid3D;
//class NodeGridDisplacements3D;
//class Voxel8N;
//class MultiGridStructure;
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... regular structure e.g. from pixel/voxel data
class StructureGrid: public virtual CallbackHandlerGrid
{
friend class MultiGridStructure;
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION
public:
    //! @brief constructor
    //! @param mDimension  Structural dimension (1,2 or 3)
    StructureGrid(int rDimension);

    StructureGrid(){}


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
    void ImportFromVtkASCIIFileHeader(std::string rFileName,size_t *rGridDimension,double *rVoxelSpacing,double *rGridOrigin, size_t rNumVoxel);

    //! @brief import routine for basic grid data with StructureGrid data space
    void ImportFromVtkASCIIFileHeader(std::string rFileName);

	//! @brief ... imports Data from a Vtk ASCII File
	//! @param fileName ... file name
	void ImportFromVtkASCIIFile(const std::string fileName,std::vector<int> &rData);

    //! @brief ... get the name of the input data file
	std::string GetInputDataFile()
	{
		return mImageDataFile;
	}

	//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    std::string GetTypeId()const
    {
        return std::string("StructureGrid");
    }
    //! @brief returns dimension
     //! @return diension
     int GetDimension() const
     {
    	 return mDimension;
     }


    //! @brief returns number of Voxels
    //! @return number of Voxels
    int GetNumVoxels() const;

    //! @brief returns  VoxelSpacing
    //! @return VoxelSpacing
    const std::vector<double> GetVoxelSpacing() const;

    //! @brief returns GridOrigin
    //! @return GridOrigin
    const std::vector<double> GetGridOrigin() const;

     //! @brief returns GridDimension
     //! @return GridDimension
    const std::vector<size_t> GetGridDimension() const;

    //! @brief returns GridDimension
     //! @return GridDimension
    void  SetGridDimension(std::vector<size_t> &rGridDimension);

     //! @brief Get NumMaterials
     //! @return NumMaterial
     const int GetNumMaterials() const;

     //! @brief Get NumBasisMaterials
     //! @return NumBasisMaterial
     const int GetNumBasisMaterials() const;

     //! @brief Set NumBasisMaterials
     void SetNumBasisMaterials(int rNumBasisMaterials);

     //! @brief get id of node from grid node number
     //! @param grid node number
     //! @return node id
     size_t GetNodeId(size_t rEdgeNumber)
     {
    	 return mNodeId[rEdgeNumber];
     }

	//! @brief Get CurrentGridNumber
	//! @return rCurrentGridNumber
	const int GetCurrentGridNumber() const;

	//! @brief Set CurrentGridNumber
	void SetCurrentGridNumber(int rCurrentGridNumber);

     //! @brief Get number of Constraints
     //! @return NumConstraints
     const size_t GetNumConstraints() const;

     void SetMatrixFreeMethod(bool rMatrixFreeMethod);

     bool GetMatrixFreeMethod();

    //! @brief set basis element stiffness
     //! @param rPoissonsRatio
     //! @param rMaterialNum ... number of material with different Poisson Ratio
     void SetBasisElementStiffnessMatrix(double rPoissonsRation,int rMaterialNum);

    //! @brief Get LocalCoefficientMatrix0
    //! @param NumLocalCoefficientMatrix0 number of stiffness matrix
    const std::vector<double>* GetLocalCoefficientMatrix0(int rNumLocalCoefficientMatrix0) const;

    //! @brief set basis edge stiffnesses
    //! @param rBasisMaterialNum ... number of material,
    void SetBasisEdgeStiffnessMatrices(int rBasisMaterialNum);

    //! @brief set general neighbor nodes for node-edge based routines
    //! @brief in function of grid dimension
    void SetNeighborNodesNE();

    //! @brief ..calculate node numbers at one element
    //! @param elementNumber ... number of the element
    //! @param nodeNumbers ... nodes at element
    void CalculateNodesAtElement(size_t elementNumber,std::vector<size_t>& nodeNumbers)const;

    //! @brief set material number at edges for node-edge based routines
    void SetMaterialNumberForEdges();

    //! @brief get DisplacementConstaints
    //! @return dynamic_bitset of constraints
    const boost::dynamic_bitset<> GetDisplacementConstaints();

    //! @brief all nodes at one plane perpendicular to direction are constrained
    //! @param rDirection ... 0 = x,1 = y,2 = z
    //! @param rGridLocation ... region of constrained nodes xmin,xmax,ymin,ymax,zmin,zmax
    //! @param rValue ... value of boundary displacement
    //! @param rDisplVector ... output initial displacement vector
    void SetDisplacementConstraints(size_t rDirection,size_t *rGridLocation,double rValue,std::vector<double> &rDisplVector);

    //! @brief create grid data
    //! @param rThresholdMaterialValue ... threshold between material one and two
    //! @param fileName ... file containing image data
    //! @param rColorToMaterialData ... vector of material data (Young's Modulus) mapped to color points
   void CreateGrid(int rThresholdMaterialValue,std::string fileName, std::vector<double>& rColorToMaterialData);

   //! @brief get weighting factor for preconditioner
   //! @return rWeight ... weighting factor
   double GetWeightingFactor()const;

   //! @brief set weighting factor for preconditioner
   //! @param rWeight ... weighting factor
   void  SetWeightingFactor(double rWeight);

   //! @brief ... calculate matrix-vector-product in element-by-element way
   //! @brief ... with local vectors
   //! @param ... u - parameters input, r - gradient output
   void CalculateMatrixVectorProductEBE(std::vector<double>& u,std::vector<double>& r)const;

   //! @brief ... calculate matrix-vector-product in node-by-node way
   //! @brief ... with local vectors
   //! @param ... u - prarmeters input, r - gradient output
   void CalculateMatrixVectorProductNBN(std::vector<double> &u,std::vector<double> &r)const;

	//! @brief ... calculate reaction forces in element-by-element way
	//! @param ... u - prarmeters input, f - forces output
	void CalculateReactionForcesEBE(std::vector<double> &u,std::vector<double> &f)
	{
	       throw MechanicsException ( "[NuTo::StructureGrid::CalculateReactionForcesEBE] Routine is not implemented." );
	}

	void GetEngineeringStrain(const std::vector<double>& rDisplacements, std::vector<double>& rEngineeringStrain)const;

	void GetEngineeringStress( std::vector<double>& rEngineeringStrain, std::vector<double>& rEngineeringStress)const;

    //! @brief ... export to Vtk datafile
    //! @param rFilename ... filename
    void ExportVTKStructuredDataFile(const std::string& rFilename) const;

	//! @brief MultiGrid routines
	//! @brief initialize on coarser grid
	//! @param StructureGrid reference
	void SetCoarserGridLevel(NuTo::StructureGrid *rCoarseGrid);

	//! @brief MultiGrid routines
	//! @brief get pointer to coarser grid
	NuTo::StructureGrid& GetCoarseGridPtr()
	{
		return mpCoarseGrid[0];
	}

	//! @brief MultiGrid routines
	//! @brief Get pointer to finer grid
	NuTo::StructureGrid& GetFineGridPtr()
	{
		return mpFineGrid[0];
	}

	//! @brief MultiGrid routines
	//! @brief set pointer to coarser grid
	void SetCoarseGridPtr(NuTo::StructureGrid* rCoarseGrid);

	//! @brief MultiGrid routines
	//! @brief set pointer to finer grid
	void SetFineGridPtr(NuTo::StructureGrid* rFineGrid);

	double GetElementYoungsModulus(size_t rElementNumber)
	{
		return mYoungsModulus[rElementNumber];
	}


 //*************************************************
//************ Node routines        ***************
//***  defined in structures/StructureGridNode.cpp  ***
//*************************************************
    //! @brief returns the number of nodes
    //! @return number of nodes
    size_t GetNumNodes() const;

//
    //! @brief info about the nodes in the Structure
    void NodeInfo(int mVerboseLevel) const;
//
//    void NodeCreate(int rNodeID, std::string rDOFs);
//
    void CreateNodeGrid(std::string rDOFs,int rThresholdMaterialValue)
    {
        throw MechanicsException ( "[NuTo::StructureGrid::CreateNodeGrid] Routine is not implemented." );
    }

    int* GetCoincidenceVoxelIDs(int rNodeID);

    //! @brief NodeGetConstraintSwitch
    //! @param rGlobalDof
    //! @return switch for constraint
    bool NodeGetConstraintSwitch(int rGlobalDof);


//*************************************************
//************ Element routines     ***************
//**  defined in structures/StructureGridElement.cpp **
//***********************************************"elementVec",**
    //! @brief returns the number of elements
    //! @return number of elements
    int GetNumElements() const;

    //! @brief info about the elements in the Structure
    void ElementInfo(int mVerboseLevel) const;

    void CreateElementGrid(NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rBaseCoefficientMatrix0,
            const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rColorToMaterialData,const std::string& rElementType);

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info()const;


	// provide functions for iterative solution process
	void Gradient(std::vector<double>& rValues,std::vector<double>& rGradient);
    void Hessian(std::vector<double>&  rDiagHessian);
    void HessianDiag(std::vector<double>&  rDiagHessian);
    void CalcScalingFactors(std::vector<double> &p);

     void SetUseDiagHessian (bool rUseDiagHessian)
     {
     	mUseDiagHessian=rUseDiagHessian;
     }
     bool GetUseDiagHessian ()
     {
     	return mUseDiagHessian;
     }
     void SetMisesWielandt (bool rMisesWielandt)
     {
     	mUseMisesWielandt=rMisesWielandt;
     }

     bool GetMisesWielandt ()
     {
     	return mUseMisesWielandt;
     }
    std::vector<double>&  GetParameters();
    std::vector<double>&  GetResidual();
    std::vector<double>&  GetExtForces();
    void SetParameters(std::vector<double>& rParameters);
    void SetResidual(std::vector<double>& rResidual);
    void SetExtForces(std::vector<double>& rExtForces);
    void CalculateMultiGridCorrolations(std::string restrictionType,std::vector<double>& rRestriction,std::vector<double>& rProlongation);
    void Restriction(std::vector<double> &rRestrictionFactor);
    void Prolongation(std::vector<double> &rProlongationFactor);
    void AnsysInput(std::vector<double >&rDisplVector) const;
    void LSDynaInput() const;
protected:


    int mDimension;
    size_t mNumVoxel;  //number of voxels
    std::vector<double> mVoxelSpacing; //spacing between center of neighbor voxels / dimension of each voxel
    std::vector<size_t> mGridDimension; //dimension of the voxel model
    std::vector<double> mGridOrigin;// origin of the model , in the center of the first voxel
//    boost::ptr_map<int,NodeBase> mNodeMap;
//    boost::ptr_map<int,ElementBase> mElementMap;
    std::string mImageDataFile;
	size_t mNumConstraintDofs;
    std::vector<std::vector<double> > mLocalCoefficientMatrix0;
    std::vector<std::vector<double> > mBasisEdgeCoefficientMatrix0;
    std::vector<double> mLocalDerivativeShapeFunctions;
    double mC11,mC12,mC44; //stiffness tensor coefficients with E=1
    std::vector<int> mNeighborNodesNE;
	std::vector<size_t> mEdgeId;		//saves the id of the edge of this fe node
	std::vector<size_t> mNodeId;		//saves the id of the fe node of this edge
	std::vector<size_t> mVoxelId;		//saves the id of the voxel of this element
	boost::dynamic_bitset<> mDofIsConstraint; //0 = false, all 0 here
	std::vector<double> mDisplacements;
    int mNumMaterials;
    std::vector<double> mResidual;
	bool mMatrixFreeMethod;
    std::vector<double> mExtForces;
    int mNumBasisMaterials;
	std::vector<double> mYoungsModulus;
	std::vector<double> mLinearElasticEngineeringStrains;
	std::vector<double> mLinearElasticEngineeringStresses;
	bool mUseDiagHessian;
	bool mUseMisesWielandt;

	//MultiGrid
	size_t mCurrentGridNumber;
	double mWeightingFactor;
	std::vector<size_t> mFineEdgeId;
	StructureGrid* mpFineGrid;
	StructureGrid*  mpCoarseGrid;
};
} //namespace NuTo
#endif // STRUCTUREGRID_H
