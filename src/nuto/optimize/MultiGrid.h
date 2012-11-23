//$Id$
//! @author Andrea Keßler, ISM
//! @date September 2012
//! @brief ... standard class for multigrid routines

#ifndef MULTIGRID_H
#define MULTIGRID_H

#include "nuto/base/NuToObject.h"

#include <boost/ptr_container/ptr_vector.hpp>
#include "nuto/optimize/ConjugateGradientGrid.h"
#include "nuto/mechanics/structures/grid/StructureGrid.h"

//#include "nuto/optimize/OptimizeException.h"
namespace NuTo
{

class MultiGrid : public virtual NuToObject
{
	friend class ConjugateGradientGrid;
public:
    //! @brief constructor
    MultiGrid():NuToObject()
    {
    	mNumGrids=1;
        mRestrictionType="TWENTYSEVENPOINTS";
        mNumRelaxIterations=5;
        mCurrentGridDimension=0;
        mCurrentGridNumber=0;
    }

    //! @brief deconstructor
    virtual ~MultiGrid()
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
    	   & BOOST_SERIALIZATION_NVP(mNumRelaxIterations)
    	   & BOOST_SERIALIZATION_NVP(mCurrentGridDimension)
    	   & BOOST_SERIALIZATION_NVP(mCurrentGridNumber);
    }
#endif  // ENABLE_SERIALIZATION

    void Info()const;
    void SetStructure(NuTo::StructureGrid* rpStructureHandler);
 	int Initialize();
 	NuTo::StructureGrid* GetGridPtr(int rIdent);
 	int Optimize();



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
    //! @brief gives number of relaxation iterations - smoothing on fine grid
    //! @return number of relaxation iterations
    int GetNumRelaxIterations();

protected:
    //! @brief set number of relaxation iterations - smoothing on fine grid
    //! @param numIterations ... number of relaxation iterations
    void SetNumRelaxIterations(int numIterations);

    //! @brief set number of relaxation iterations - smoothing on fine grid


private:
    StructureGrid* mpStructureHandler; // for unregular structures an extra pointer has to be added
    boost::ptr_vector<StructureGrid> mpGridPtr; // for unregular structures an extra pointer has to be added
    //! @brief number of grids
    int mNumGrids;
    //! @brief restriction type
    std::string mRestrictionType;
    //! @brief number of relaxation iterations
    int mNumRelaxIterations;
    //! @brief current grid dimensions (elements) in each direction
    int *mCurrentGridDimension;
    //! @brief stores current grid number - gives actual level information
    //! @brief finest grid = 0,
    int mCurrentGridNumber;
    std::vector<double> mRestriction; // is only for one direction
    std::vector<double> mProlongation; // is only for one direction

};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::MultiGrid)
#endif // SWIG
#endif // ENABLE_SERIALIZATION

#endif // MULTIGRID_H
