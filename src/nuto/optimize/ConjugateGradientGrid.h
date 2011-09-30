// $Id $
#ifndef CONJUGATE_GRADIENT_GRID_H
#define CONJUGATE_GRADIENT_GRID_H

#ifdef ENABLE_EIGEN
//#include <eigen2/Eigen/Core>
//#include <eigen3/Eigen/Dense>
#endif


#include <vector>
#include <boost/dynamic_bitset.hpp>
#include <memory>

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
#ifdef ENABLE_EIGEN
typedef Eigen::VectorXd myType;
#else
typedef std::vector<double> myType;
#endif
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

    void Initialize(int rNumParameters,boost::dynamic_bitset<>& elemExist, boost::dynamic_bitset<>& nodeExist,boost::dynamic_bitset<> &rDofIsConstraint,std::vector<double>& youngsModulus,std::vector<double>&  baseStiffness,std::vector<int>& materialOfElem,std::vector<int>& allNodesAtVoxel,std::vector<double>& parameters)
    {
    	mNumParameters=rNumParameters;
        mElemExist=elemExist;
       	mNodeExist=nodeExist;
       	mDofIsConstraint=rDofIsConstraint;
       	mYoungsModulus=youngsModulus;
       	mBaseStiffness=baseStiffness;
       	mMaterialOfElem=materialOfElem;
        mNodeIds=allNodesAtVoxel;
        mParameters=parameters;
    }

    void AnsysInput(int rNumParameters,boost::dynamic_bitset<>& elemExist, boost::dynamic_bitset<>& nodeExist,boost::dynamic_bitset<> &rDofIsConstraint,std::vector<double>& youngsModulus,int* rGridDimension,double* rVoxelSpacing,std::vector<int>& materialOfElem,std::vector<int>& allNodesAtVoxel,std::vector<double>& parameters);

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

    void GetParameters(std::vector<double>& rParameters)
    {
		rParameters.assign(mParameters.begin(),mParameters.end());
    }

#ifdef ENABLE_EIGEN
    void GetParameters(myType& rParameters)
    {
    	rParameters.Map(mParameters.data(),mNumParameters);
    }

    void SetParameters(myType& rParameters)
    {
    	mParameters.assign(mNumParameters,rParameters.coeffRef(0));
    }
#endif

    void SetParameters(std::vector<double>& rParameters)
    {
    	mParameters.assign(rParameters.begin(),rParameters.end());
    }

    void HessianDiag(myType&  rHessianDiag)const;

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
	void CalcScalingFactors(int& numHessianCalls, myType &p);

	//! @brief ... calculate start gradient in element-by-element way
	//! @param ... u - prarmeters input, r - gradient output
	void CalculateMatrixVectorProductEBE(myType &u,myType &r);

	//! @brief ... calculate start gradient in node-by-node way
	void CalculateStartGradientNodeByNode(myType &r);


	double mAccuracyGradient;
	double mMinDeltaObjBetweenRestarts;
	int    mMaxGradientCalls;
	int    mMaxHessianCalls;
	int    mMaxIterations;
	int    mShowSteps;
    bool   mUseDiagHessian;
    bool   mUseMultiGrid;

   	int mNumParameters;
   	boost::dynamic_bitset<> mElemExist;
   	boost::dynamic_bitset<> mNodeExist;
   	boost::dynamic_bitset<> mDofIsConstraint;
   	std::vector<double> mYoungsModulus;
   	std::vector<double> mBaseStiffness;
   	std::vector<int> mMaterialOfElem;
    std::vector<int> mNodeIds;
    std::vector<double> mParameters;
};
} // namespace NuTo
#endif // CONJUGATE_GRADIENT_GRID_H
