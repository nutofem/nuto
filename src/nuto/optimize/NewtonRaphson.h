// $Id$

#ifndef NEWTONRAPHSON_H
#define NEWTONRAPHSON_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/optimize/NonlinearSolverBase.h"

namespace NuTo
{
//! @author Vitaliy Kindrachuk
//! @date December 2015
//! @brief ... standard class for solving system of nonlinear equations
class NewtonRaphson : public NonlinearSolverBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NewtonRaphson();

//    //! @brief sets the pointer to the analytic derivative of the residual function
//    //! @param rParam ... parameters necessary to evaluate the residual
//    //! @param rUnknown ... unknown vector
//    void SetResidualDerivativeFunction(
//    		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>(*rResidualDerivativeFunction)
//			(const NuTo::FullVector<double,Eigen::Dynamic>&,NuTo::FullVector<double,Eigen::Dynamic>))
//    {
//    	mResidualDerivativeFunction = rResidualDerivativeFunction;
//    }

    //! @brief sets the pointer to the analytic derivative of the residual function
    //! @param rParam ... parameters necessary to evaluate the residual
    //! @param rUnknown ... unknown vector
    void SetResidualDerivativeFunction(
    		boost::function<NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>
    		    	(const NuTo::FullVector<double,Eigen::Dynamic>&,NuTo::FullVector<double,Eigen::Dynamic>)> rResidualDerivativeFunction)
    {
    	mResidualDerivativeFunctionBoost = rResidualDerivativeFunction;
    	mAssignResidualResidualDerivative = true;
    }

    //! @brief returns mCheckNewtonRaphson, specifying iterator exit
    double GetCheckNewtonRaphson()const
    {
    	return mCheckNewtonRaphson;
    }

#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif// SWIG
#endif // ENABLE_SERIALIZATION

    //! @brief perform iteration
    //! @param rUnknown ... unknown vector
    NuTo::Error::eError Solve(NuTo::FullVector<double,Eigen::Dynamic> &rUnknown);

    //! @brief ... the routine performs line search correction of the Newton step
    void LineSearch(NuTo::FullVector<double,Eigen::Dynamic> &rXold, const double rFold, NuTo::FullVector<double,Eigen::Dynamic> &rG, NuTo::FullVector<double,Eigen::Dynamic> &rP,
    		NuTo::FullVector<double,Eigen::Dynamic> &rX, double &rF, const double rStpmax, bool &rCheck, NuTo::FullVector<double,Eigen::Dynamic> &rFvec) const;

    //! @brief ... the routine performs Newton-Raphson integration
    NuTo::Error::eError NewtonRaphsonIterator(NuTo::FullVector<double,Eigen::Dynamic> &x, bool &check) const;

    //! @brief ... Return the name of the class, this is important for the serialize routines
    //! @return    class name
    virtual std::string GetTypeId()const;

protected:
    //pointer to the analytical derivative of the residual function (analytic Jacobi matrix)
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> (*mResidualDerivativeFunction)
    	(const NuTo::FullVector<double,Eigen::Dynamic>&,NuTo::FullVector<double,Eigen::Dynamic>);

    boost::function<NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>
    	(const NuTo::FullVector<double,Eigen::Dynamic>&,NuTo::FullVector<double,Eigen::Dynamic>)> mResidualDerivativeFunctionBoost;

    //a boolean variable which gets true if a jacobi function is assigned
    bool mAssignResidualResidualDerivative;

    //! @brief ... specify the exit: false on a normal exit; true if this is a local minimum of Fmin (that is the resudual function is not necessary zeroed)
    bool mCheckNewtonRaphson;
};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::NewtonRaphson)
#endif // SWIG
#endif // ENABLE_SERIALIZATION



#endif // NewtonRaphson_H
