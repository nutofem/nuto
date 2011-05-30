// $Id$
#ifndef CONJUGATE_GRADIENT_GRID_H
#define CONJUGATE_GRADIENT_GRID_H

#include <vector>
#include "nuto/math/FullMatrix.h"

#include "nuto/optimize/Optimizer.h"
#include "nuto/optimize/OptimizeException.h"
#include "nuto/optimize/CallbackHandlerGrid.h"

namespace NuTo
{
//! @author Andrea Keszler, ISM
//! @date July 2010
//! @brief ... standard class for element-by-element based conjugate gradient method
#ifdef ENABLE_MECHANICS
	class StructureGrid;
#endif // ENABLE_MECHANICS

class ConjugateGradientGrid : public Optimizer
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    ConjugateGradientGrid(unsigned int rNumParameters) : Optimizer(rNumParameters,(unsigned int)0,(unsigned int) 0)
    {
        mAccuracyGradient = 1e-6;
        mMinDeltaObjBetweenRestarts = 1e-6;
        mMaxGradientCalls = INT_MAX,
        mMaxHessianCalls = INT_MAX,
        mMaxIterations = INT_MAX;
        mShowSteps = 100;
       	mUseDiagHessian =true;
       	mUseMultiGrid= false;
	}

#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
    	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Optimizer)
    	   & BOOST_SERIALIZATION_NVP(mAccuracyGradient)
           & BOOST_SERIALIZATION_NVP(mMinDeltaObjBetweenRestarts)
           & BOOST_SERIALIZATION_NVP(mMaxGradientCalls)
           & BOOST_SERIALIZATION_NVP(mMaxHessianCalls)
           & BOOST_SERIALIZATION_NVP(mMaxIterations)
           & BOOST_SERIALIZATION_NVP(mShowSteps)
           & BOOST_SERIALIZATION_NVP(mUseDiagHessian)
           & BOOST_SERIALIZATION_NVP(mUseMultiGrid);
    }
#endif // SWIG
#endif // ENABLE_SERIALIZATION

	int Optimize();

    inline void SetMaxGradientCalls(int rMaxGradientCalls)
    {
        mMaxGradientCalls = rMaxGradientCalls;
    }

    inline void SetMaxHessianCalls(int rMaxHessianCalls)
    {
        mMaxHessianCalls = rMaxHessianCalls;
    }

    inline void SetMaxIterations(int rMaxIterations)
    {
        mMaxIterations = rMaxIterations;
    }

    inline void SetAccuracyGradient(double rAccuracyGradient)
    {
        mAccuracyGradient = rAccuracyGradient;
    }

    inline void SetMinDeltaObjBetweenRestarts(double rMinDeltaObjBetweenRestarts)
    {
        mMinDeltaObjBetweenRestarts = rMinDeltaObjBetweenRestarts;
    }

    inline void SetShowSteps(int rShowSteps)
    {
        mShowSteps = rShowSteps;
    }

    inline void SetUseMultiGrid(bool rUseMultiGrid)
    {
    	mUseMultiGrid = rUseMultiGrid;
    }

#ifdef ENABLE_MECHANICS
    inline void SetGridStructure(NuTo::StructureGrid* rpGrid)
    {
			mpGrid = rpGrid;
    }
#else
    inline void SetGridStructure()
    {
		throw OptimizeException ( "[ConjugateGradientGrid::SetGridStructure] Modul Mechanics is not loaded." );
    }
#endif // ENABLE_MECHANICS

    void Hessian(NuTo::FullMatrix<double>&  rHessian)const;
    void HessianDiag(NuTo::FullMatrix<double>&  rHessian)const;

#ifdef ENABLE_SERIALIZATION
    //! @brief ... save the object to a file
    //! @param filename ... filename
    //! @param rType ... type of file, either BINARY, XML or TEXT
    void Save ( const std::string &filename, std::string rType)const;


    //! @brief ... restore the object from a file
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
    void Restore ( const std::string &filename,  std::string rType);
#endif // ENABLE_SERIALIZATION

    //! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    virtual std::string GetTypeId()const;

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
	virtual void Info()const;


protected:
	void CalcScalingFactors(int& numHessianCalls, Eigen::VectorXd& scaleFactorsInv);

	//! @brief ... calculate start gradient in element-by-element way
	void CalculateStartGradient(NuTo::FullMatrix<double> &gradientOrig);

	//! @brief ... calculate start gradient in element-by-element way
	void CalculateStartGradientNodeByNode(NuTo::FullMatrix<double> &gradientOrig);

	//! @brief ... calculate matix-vector product in element-by-element way
	void CalculateMatrixVectorEBE(bool startSolution, NuTo::FullMatrix<double> &returnVector);

	//! @brief ... calculate scaled search direction multiplied with stiffness matrix in element-by-element way for each step
	void CalculateScaledSearchDirection(Eigen::VectorXd& searchDirectionScaled);

	//! @brief ... calculate search direction in node-by-node way
	void CalculateScaledSearchDirectionNodeByNode(Eigen::VectorXd& searchDirectionScaled);

	#ifdef ENABLE_MECHANICS
		StructureGrid *mpGrid;
	#endif // ENABLE_MECHANICS

	double mAccuracyGradient;
	double mMinDeltaObjBetweenRestarts;
	int    mMaxGradientCalls;
	int    mMaxHessianCalls;
	int    mMaxIterations;
	int    mShowSteps;
    bool   mUseDiagHessian;
    bool   mUseMultiGrid;

};
} // namespace NuTo
#endif // CONJUGATE_GRADIENT_GRID_H
