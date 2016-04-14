
#ifndef NONLINEARSOLVERBASE_H
#define NONLINEARSOLVERBASE_H


#include <ctime>
#include <array>

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include <boost/ptr_container/ptr_map.hpp>
#include <boost/function.hpp>


#include "nuto/base/NuToObject.h"
#include "nuto/base/ErrorEnum.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/FullVector.h"

//#include "nuto/mechanics/constitutive/mechanics/GradientDamageEngineeringStressFatigue.h"
//#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataGradientDamage1DFatigue.h"



namespace NuTo
{
//! @author Kindrachuk
//! @date December 2015
//! @brief ... standard abstract class for all solvers of nonlinear systems of equations
class NonlinearSolverBase : public NuToObject
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief constructor
    NonlinearSolverBase();

    //! @brief deconstructor
    virtual ~NonlinearSolverBase()
    {}

    //! @brief perform solving
    //! @param rUnknown ... unknown vector
    virtual NuTo::Error::eError Solve(NuTo::FullVector<double,Eigen::Dynamic> &rUnknown)=0;

    //! @brief sets the pointer to the residual function
    //! @param rParam ... parameters necessary to evaluate the residual
    //! @param rUnknown ... unknown vector
//    void SetResidualFunction(
//    		NuTo::FullVector<double,Eigen::Dynamic> (*rResidualFunction)(
//    		const NuTo::FullVector<double,Eigen::Dynamic>&,
//    		NuTo::FullVector<double,Eigen::Dynamic>))
//    {
//    	mResidualFunction = rResidualFunction;
//    }

    //! @brief sets the pointer to the residual function
    //! @param rParam ... parameters necessary to evaluate the residual
    //! @param rUnknown ... unknown vector
    void SetResidualFunction(
    		boost::function<NuTo::FullVector<double,Eigen::Dynamic>
    (const NuTo::FullVector<double,Eigen::Dynamic>&,NuTo::FullVector<double,Eigen::Dynamic>)> rResidualFunction)
    {
    	mResidualFunctionBoost = rResidualFunction;
    	mAssignResidual = true;
    }

    //! @brief sets the tolerance for the residual vector, (the components of the residual should be of the same magnitude of order)
    void SetResidualTolerance(double rTolResidual)
    {
        mTolResidual = rTolResidual;
    }

    //! @brief returns the tolerance for the residual vector
    double GetResidualTolerance()const
    {
    	return mTolResidual;
    }

    //! @brief sets the tolerance for the solution vector, (the components of the rUnknown vector should be of the same magnitude of order)
    void SetSolutionTolerance(double rTolSolution)
    {
        mTolSolution = rTolSolution;
    }

    //! @brief returns the tolerance for the solution rUnknown
    double GetSolutionTolerance()const
    {
    	return mTolSolution;
    }

    //! @brief sets the  maximal number of iterations
    void SetMaxIterations(int rMaxIterations)
    {
        mMaxIterationsNumber = rMaxIterations;
    }

    //! @brief returns the tolerance for the residual vector
    double GetMaxIterations()const
    {
    	return mMaxIterationsNumber;
    }

    //! @brief sets the list of parameters mParameter necessary to evaluate mResidualFunction
    void SetParameters(NuTo::FullVector<double,Eigen::Dynamic> &rParameter)
    {
    	mParameter = rParameter;
    }

    //! @brief ... numerical differentiation of the residual function mResidualFunction (numerical Jacobi matrix)
    //! @param rUnknown ... position at which the derivative is taken
    //! @param rFvec ... residual vector at the position rUnknown, rFvec = mResidualFunction(rParameter,rUnknown)
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> DResidualNum(NuTo::FullVector<double,Eigen::Dynamic> rUnknown,
    		NuTo::FullVector<double,Eigen::Dynamic> &rFvec) const;

    //! @brief ... calculates 0.5*rFvec^2 and updates rFvec = mResidualFunction(mParameter,rUnknown)
    double Fmin(NuTo::FullVector<double,Eigen::Dynamic> rUnknown, NuTo::FullVector<double,Eigen::Dynamic> &rFvec) const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    virtual void Info()const;

protected:
    //pointer to the residual function
    NuTo::FullVector<double,Eigen::Dynamic> (*mResidualFunction)
	(const NuTo::FullVector<double,Eigen::Dynamic>&,NuTo::FullVector<double,Eigen::Dynamic>);

    boost::function<NuTo::FullVector<double,Eigen::Dynamic> (const NuTo::FullVector<double,Eigen::Dynamic>&,NuTo::FullVector<double,Eigen::Dynamic>)> mResidualFunctionBoost;

    //tolerance for the residual vector
    double mTolResidual;

    //tolerance for the solution vector
    double mTolSolution;

    //maximal allowed number of iterations
    int mMaxIterationsNumber;

    //a boolean variable which gets true if a residual function is assigned
    bool mAssignResidual;

    //list of parameters necessary for evaluation of the resuduum
    NuTo::FullVector<double,Eigen::Dynamic> mParameter;

};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::NonlinearSolverBase)
#endif // SWIG
#endif  // ENABLE_SERIALIZATION
#endif // NONLINEARSOLVERBASE_H
