//$Id$
//! @author Andrea Ke�ler, ISM
//! @date March 2013
//! @brief ... standard class for multigrid structure

#ifndef MULTIGRIDSTRUCTURE_H
#define MULTIGRIDSTRUCTURE_H

#include "nuto/base/NuToObject.h"

#include <boost/ptr_container/ptr_vector.hpp>
#ifdef ENABLE_OPTIMIZE
#include "nuto/optimize/CallbackHandlerGrid.h"
#include "nuto/optimize/ConjugateGradientGrid.h"
#include "nuto/optimize/Jacobi.h"
#endif //ENABLE_OPTIMIZE
#include "nuto/mechanics/structures/grid/StructureGrid.h"
#include "nuto/mechanics/structures/grid/OctreeGrid.h"

class StructureGrid;
namespace NuTo
{
#ifdef ENABLE_OPTIMIZE
class MultiGridStructure : public virtual CallbackHandlerGrid
#else //ENABLE_OPTIMIZE
class MultiGridStructure : public virtual NuToObject
#endif //ENABLE_OPTIMIZE
{
public:
   typedef enum
	{
		MAXCYCLES,   //maximum number of function calls is reached
		NORMGRADIENT,       //norm of gradient is smaller than a prescribed value
	} optimization_return_attributes;

	 //! @brief constructor
	MultiGridStructure()
    {
		mNumGrids=1;
        mCurrentGridDimension=0;
        mCurrentGridNumber=0;
       	mRestrictionType="TWENTYSEVENPOINTS";
        mNumPreSmoothingSteps=2;
        mNumPostSmoothingSteps=1;
        mMaxCycles=1;
        mMaxGrids=10;

    }

    //! @brief deconstructor
    virtual ~MultiGridStructure()
//    ~MultiGrid()
    {}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NuToObject)
    	   & BOOST_SERIALIZATION_NVP(mNumGrids)
    	   & BOOST_SERIALIZATION_NVP(mCurrentGridDimension)
    	   & BOOST_SERIALIZATION_NVP(mCurrentGridNumber)
    	   & BOOST_SERIALIZATION_NVP(mRestrictionType)
    	   & BOOST_SERIALIZATION_NVP(mNumPreSmoothingSteps)
    	   & BOOST_SERIALIZATION_NVP(mMaxCycles);
      }
#endif  // ENABLE_SERIALIZATION

    void Info()const;
    void SetStructure(NuTo::StructureGrid* rpStructureHandler);
 	int Initialize();
 	int MultiGridSolve(std::vector<double>& rSolution,std::vector<double>& rRightHandSide);
 	void GridLevelSolve(int rGridLevel,std::vector<double>& rSolution,std::vector<double>& rRightHandSide);
 	StructureGrid* GetGridPtr(int rIdent);

   #ifdef ENABLE_SERIALIZATION
    void Save (const std::string &filename, std::string rType )const{}
    void Restore (const std::string &filename, std::string rType ){}
    #endif // ENABLE_SERIALIZATION

    std::string GetTypeId()const
    {
    	return std::string("MultiGrid");
    }

    //! @brief gives number of grids
    //! @return number of grids
    int GetNumGrids();

    //! @brief set current grid number
    //! @param rGridNumber ... grid number
    void SetCurrentGridNumber(int rGridNumber);

	void Gradient(std::vector<double>& rValues,std::vector<double>& rGradient);
    void Hessian(std::vector<double>&  rDiagHessian);

    void HangingNodesCorrection(std::vector<double>& u)
    {
 	   //not needed her, only for octree
    }

    void SetUseMultiGridAsPreconditoner(bool rUseAsPreconditioner);
    bool GetUseMultiGridAsPreconditoner();

    void SetMaxCycle(int rMaxCycles);
    int GetMaxCycle();

    void SetMaxGrids(int rMaxGrids);
    int GetMaxGrids();
    
    std::vector<double>&  GetParameters();
    std::vector<double>&  GetRightHandSide();
    void SetParameters(std::vector<double>& rParameters);
    void SetRightHandSide(std::vector<double>& rRightHandSide);

    //! @brief gives number of relaxation iterations - smoothing on fine grid
    //! @return number of relaxation iterations
    int GetNumPreSmoothingSteps();

    //! @brief set number of relaxation iterations - smoothing on fine grid
    //! @param numIterations ... number of relaxation iterations
    void SetNumPreSmoothingSteps(int numIterations);

    //! @brief gives number of relaxation iterations - smoothing on fine grid
    //! @return number of relaxation iterations
    int GetNumPostSmoothingSteps();

    //! @brief set number of relaxation iterations - smoothing on fine grid
    //! @param numIterations ... number of relaxation iterations
    void SetNumPostSmoothingSteps(int numIterations);

    // export to Vtk Datafile
    void ExportVTKStructuredDataFile(int rGridLevel,const std::string& rFilename);

protected:


private:
    StructureGrid* mpStructureHandler; // for unregular structures an extra pointer has to be added
    OctreeGrid* mpOctreeHandler; // for unregular structures an extra pointer has to be added
    boost::ptr_vector<StructureGrid> mpGridPtr; // for unregular structures an extra pointer has to be added
    boost::ptr_vector<OctreeGrid> mpOctreePtr; // for unregular structures an extra pointer has to be added
    //! @brief number of grids
    int mNumGrids;
    //! @brief current grid dimensions (elements) in each direction
    int *mCurrentGridDimension;
    //! @brief stores current grid number - gives actual level information
    //! @brief finest grid = 0,
    int mCurrentGridNumber;
	//! @brief restriction type
	std::string mRestrictionType;
	//! @brief number of relaxation iterations
	int mNumPreSmoothingSteps;
	int mNumPostSmoothingSteps;
    std::vector<double> mRestriction; // is only for one direction
    std::vector<double> mProlongation; // is only for one direction
    std::vector<int> mCycle; // save gridnumber, and 1 for restriction and -1 for prolongation
    bool mUseAsPreconditioner;
    int mMaxCycles;
    int mMaxGrids;
};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::MultiGridStructure)
#endif // SWIG
#endif // ENABLE_SERIALIZATION

#endif // MULTIGRID_H
