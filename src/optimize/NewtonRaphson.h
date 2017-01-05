// $Id$

#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#endif // ENABLE_SERIALIZATION

#include "optimize/NonlinearSolverBase.h"

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
//			(const Eigen::VectorXd&,Eigen::VectorXd))
//    {
//    	mResidualDerivativeFunction = rResidualDerivativeFunction;
//    }

    //! @brief sets the pointer to the analytic derivative of the residual function
    //! @param rParam ... parameters necessary to evaluate the residual
    //! @param rUnknown ... unknown vector
    void SetResidualDerivativeFunction(
    		boost::function<NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>
    		    	(const Eigen::VectorXd&,Eigen::VectorXd)> rResidualDerivativeFunction)
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
    void serialize(Archive & ar, const unsigned int version){}
#endif// SWIG

    //////! @brief ... save the object to a file
    //////! @param filename ... filename
    //////! @param rType ... type of file, either BINARY, XML or TEXT
    void Save ( const std::string &filename, std::string rType)const
    {

    }
    //////! @brief ... restore the object from a file
    //////! @param filename ... filename
    //////! @param aType ... type of file, either BINARY, XML or TEXT
    void Restore ( const std::string &filename,  std::string rType)
    {

    }

#endif // ENABLE_SERIALIZATION


    //! @brief perform iteration
    //! @param rUnknown ... unknown vector
    NuTo::eError Solve(Eigen::VectorXd &rUnknown);

    //! @brief ... the routine performs line search correction of the Newton step
    void LineSearch(Eigen::VectorXd &rXold, const double rFold, Eigen::VectorXd &rG, Eigen::VectorXd &rP,
    		Eigen::VectorXd &rX, double &rF, const double rStpmax, bool &rCheck, Eigen::VectorXd &rFvec) const;

    //! @brief ... the routine performs Newton-Raphson integration
    NuTo::eError NewtonRaphsonIterator(Eigen::VectorXd &x, bool &check) const;

    //! @brief ... Return the name of the class, this is important for the serialize routines
    //! @return    class name
    virtual std::string GetTypeId()const;

protected:
    //pointer to the analytical derivative of the residual function (analytic Jacobi matrix)
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> (*mResidualDerivativeFunction)
    	(const Eigen::VectorXd&,Eigen::VectorXd);

    boost::function<NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>
    	(const Eigen::VectorXd&,Eigen::VectorXd)> mResidualDerivativeFunctionBoost;

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



