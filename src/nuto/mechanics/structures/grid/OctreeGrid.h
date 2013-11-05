// $Id $
#ifndef OCTREEGRID_H
#define OCTREEGRID_H
#include "nuto/optimize/CallbackHandlerGrid.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/math/FullMatrix.h"
#include <boost/dynamic_bitset.hpp>
#include <stdint.h>
#include <map>

namespace NuTo
{
struct data
{
	uint32_t id;
	uint32_t level;
	uint32_t constraint=0;
	double weight;
};

//! @brief forward declaration to speed up compilation time
//class NodeGrid3D;
//class NodeGridDisplacements3D;
//class Voxel8N;
//class MultiGridStructure;
//! @author Andrea Keßler, ISM
//! @date Juni 2013
//! @brief ... regular structure e.g. from pixel/voxel data
class OctreeGrid: public virtual CallbackHandlerGrid
{
friend class MultiGridStructure;
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION
public:
    //! @brief constructor
    //! @param mDimension  Structural dimension (1,2 or 3)
    OctreeGrid(int rDimension);

    OctreeGrid(){}


    ~OctreeGrid();

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
        return std::string("OctreeGrid");
    }
    //! @brief returns dimension
     //! @return diension
     int GetDimension() const
     {
    	 return mDimension;
     }


    //! @brief returns number of Voxels
    //! @return number of Voxels
    size_t GetNumVoxels() const;

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

    //! @brief all nodes at one plane perpendicular to direction are constrained
    //! @param rDirection ... 0 = x,1 = y,2 = z
    //! @param rGridLocation ... region of constrained nodes xmin,xmax,ymin,ymax,zmin,zmax
    //! @param rValue ... value of boundary displacement
    //! @param rDisplVector ... output initial displacement vector
    void SetDisplacementConstraints(size_t rDirection,size_t *rGridLocation,double rValue,std::vector<double> &rDisplVector);

   //! @brief get weighting factor for preconditioner
   //! @return rWeight ... weighting factor
   double GetWeightingFactor()const;

   //! @brief set weighting factor for preconditioner
   //! @param rWeight ... weighting factor
   void  SetWeightingFactor(double rWeight);

   //! @brief set maximal number of octree levels
   //! @param rNumLevels ... number of levels
   inline void  SetMaxOctreeLevels(uint32_t rNumLevels)
   {
	   mNumLevels=rNumLevels;
   }
   //! @brief ... calculate matrix-vector-product in element-by-element way
   //! @brief ... with local vectors
   //! @param ... u - parameters input, r - gradient output
   void CalculateMatrixVectorProductEBE(std::vector<double>& u,std::vector<double>& r)const;

   //! @brief ... calculate reaction forces in element-by-element way
	//! @param ... u - prarmeters input, f - forces output
	void CalculateReactionForcesEBE(std::vector<double> &u,std::vector<double> &f)
	{
		throw MechanicsException ( "[NuTo::OctreeGrid::CalculateReactionForcesEBE] Routine is not implemented." );
	}

	void BuildGlobalCoefficientMatrix(std::vector<double>& rKglob,std::vector<double>& rVector)const;

	void GetEngineeringStrain(const std::vector<double>& rDisplacements, std::vector<double>& rEngineeringStrain)const;

	void GetEngineeringStress( std::vector<double>& rEngineeringStrain, std::vector<double>& rEngineeringStress)const;

    //! @brief ... export to Vtk rectilinear grid datafile
    //! @param rFilename ... filename
    void ExportVTKUnstructuredGridDataFile(const std::string& rFilename) const;

	//! @brief MultiGrid routines
	//! @brief initialize on coarser grid
	//! @param OctreeGrid reference
	void SetCoarserGridLevel(NuTo::OctreeGrid *rCoarseGrid);

	//! @brief MultiGrid routines
	//! @brief get pointer to coarser grid
	NuTo::OctreeGrid& GetCoarseGridPtr()
	{
		return mpCoarseGrid[0];
	}

	//! @brief MultiGrid routines
	//! @brief Get pointer to finer grid
	NuTo::OctreeGrid& GetFineGridPtr()
	{
		return mpFineGrid[0];
	}

	//! @brief MultiGrid routines
	//! @brief set pointer to coarser grid
	void SetCoarseGridPtr(NuTo::OctreeGrid* rCoarseGrid);

	//! @brief MultiGrid routines
	//! @brief set pointer to finer grid
	void SetFineGridPtr(NuTo::OctreeGrid* rFineGrid);

	double GetElementYoungsModulus(size_t rElementNumber)
	{
		return mData[rElementNumber].weight;
	}

    //! @brief returns the number of nodes
    //! @return number of nodes
    size_t GetNumNodes() const;

    //! @brief returns the number of elements
    //! @return number of elements
    int GetNumElements() const;

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info()const;

    //! @brief Approximate condition number through max and min diagonal singular values
    //! @return approximate condition number
    double ApproximateSystemConditionNumber();

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
    void CreateOctree(int rThresholdMaterialValue, std::string fileName,
   		 std::vector<double>& rColorToMaterialData);
    void HangingNodesSearch();
	//! @brief correct solution for hanging nodes
	//! @param displacement solution
	void HangingNodesCorrection(std::vector<double>& u)const;
	void AnsysInput(std::vector<double> &rDisplVector) const;


protected:


    int mDimension;
    size_t mNumVoxel;  //number of voxels
    size_t mNumElements;  //number of elements
    std::vector<double> mVoxelSpacing; //spacing between center of neighbor voxels / dimension of each voxel
    std::vector<size_t> mGridDimension; //dimension of the voxel model
    std::vector<double> mGridOrigin;// origin of the model , in the center of the first voxel
    std::vector<size_t> mOctreeDimension; //dimension of the voxel model
    std::string mImageDataFile;
    bool mMortonOrderIsTrue;
    std::map<uint32_t,data> mData;
    uint32_t mNumLevels;
    std::map<uint32_t,size_t> mKey; //map:  maps bitmap key to number of existing data
	size_t mNumConstraintDofs;
    std::vector<std::vector<double> > mLocalCoefficientMatrix0;
    std::vector<std::vector<double> > mBasisEdgeCoefficientMatrix0;
    std::vector<double> mLocalDerivativeShapeFunctions;
    double mC11,mC12,mC44; //stiffness tensor coefficients with E=1
 	boost::dynamic_bitset<> mDofIsConstraint; //0 = false, all 0 here
	std::vector<double> mDisplacements;
    int mNumMaterials;
    std::vector<double> mResidual;
	bool mMatrixFreeMethod;
    std::vector<double> mExtForces;
    int mNumBasisMaterials;
    std::vector<double> mHessianDiag;
    double mPoissonsRatio;
	std::vector<double> mLinearElasticEngineeringStrains;
	std::vector<double> mLinearElasticEngineeringStresses;
	bool mUseDiagHessian;
	bool mUseMisesWielandt;

	//MultiGrid
	size_t mCurrentGridNumber;
	double mWeightingFactor;
	std::vector<size_t> mFineEdgeId;
	OctreeGrid* mpFineGrid;
	OctreeGrid*  mpCoarseGrid;
};
} //namespace NuTo
#endif // OCTREEGRID_H
