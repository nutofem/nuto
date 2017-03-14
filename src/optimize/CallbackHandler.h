#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include <vector>
#include "optimize/OptimizeException.h"
#include <eigen3/Eigen/Core>

namespace NuTo
{

//! @author Joerg F. Unger, ISM
//! @date September 2009
//! @brief Abstract class to handle callback routines
class CallbackHandler
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    CallbackHandler() = default;
    virtual ~CallbackHandler() = default;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template <class Archive> void serialize(Archive& ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

    virtual void SetParameters(const Eigen::MatrixXd& rParameters)
    {
        throw OptimizeException(
                __PRETTY_FUNCTION__, "SetParameters function not implemented in CallbackHandler object.");
    }

    virtual double Objective() const
    {
        throw OptimizeException(__PRETTY_FUNCTION__, "Objective function not implemented in CallbackHandler object.");
    }

    virtual void Gradient(Eigen::MatrixXd& rGradient) const
    {
        throw OptimizeException(__PRETTY_FUNCTION__, "Gradient function not implemented in CallbackHandler object.");
    }

    virtual void Gradient(std::vector<double>& rValues, std::vector<double>& rGradient) const
    {
        throw OptimizeException(__PRETTY_FUNCTION__, "Gradient function not implemented in CallbackHandler object.");
    }

    virtual void Hessian(Eigen::MatrixXd& rHessian) const
    {
        throw OptimizeException(__PRETTY_FUNCTION__, "Hessian function not implemented in CallbackHandler object.");
    }

    virtual void Hessian(std::vector<double>& rDiagHessian) const
    {
        throw OptimizeException(__PRETTY_FUNCTION__, "Hessian function not implemented in CallbackHandler object.");
    }

    virtual void Info() const = 0;
#ifdef ENABLE_SERIALIZATION
    virtual void Save(const std::string& filename, std::string rType) const {}
    virtual void Restore(const std::string& filename, std::string rType) {}
#endif // ENABLE_SERIALIZATION

};
} // namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::CallbackHandler)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::CallbackHandler)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
