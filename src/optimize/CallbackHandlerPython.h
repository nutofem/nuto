#pragma once

#include "optimize/CallbackHandler.h"

#include <Python.h>

namespace NuTo
{

//! @author Joerg F. Unger, ISM
//! @date September 2009
//! @brief Abstract class to handle callback routines
class CallbackHandlerPython : public CallbackHandler
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    CallbackHandlerPython() : CallbackHandler()
    {
        mCallbackSetParameters = nullptr;
        mCallbackObjective = nullptr;
        mCallbackGradient = nullptr;
        mCallbackHessian = nullptr;
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template <class Archive> void serialize(Archive& ar, const unsigned int version)
    {
        ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(CallbackHandler);
        mCallbackSetParameters = 0;
        mCallbackObjective = 0;
        mCallbackGradient = 0;
        mCallbackHessian = 0;
        /* the serializaion of the function pointers in python makes no sense
         * since in the restoring phase in another pointer session the python function might
         * be located at a totally different position

           & BOOST_SERIALIZATION_NVP(mCallbackSetParameters)
           & BOOST_SERIALIZATION_NVP(mCallbackObjective)
           & BOOST_SERIALIZATION_NVP(mCallbackGradient)
           & BOOST_SERIALIZATION_NVP(mCallbackHessian);
         */
    }
#endif // ENABLE_SERIALIZATION

    void SetCallbackFunctions(
            PyObject* args_parameters, PyObject* args_objective, PyObject* args_gradient, PyObject* args_hessian);

    void SetParameters(const Eigen::MatrixXd& rParameters) override;

    double Objective() const override;

    void Gradient(Eigen::MatrixXd& rGradient) const override;

    void Hessian(Eigen::MatrixXd& rHessian) const override;

    void Info() const override;

#ifdef ENABLE_SERIALIZATION
    //! @brief Restore the object from a file
    //! @param filename Filename
    //! @param aType Type of file, either BINARY, XML or TEXT
    //! @brief Save the object to a file
    virtual void Restore(const std::string& filename, std::string rType) {}

    //  @brief This routine has to be implemented in the final derived classes, which are no longer abstract
    //! @param filename Filename
    //! @param aType Type of file, either BINARY, XML or TEXT
    virtual void Save(const std::string& filename, std::string rType) const {}
#endif // ENABLE_SERIALIZATION

private:
    PyObject* mCallbackSetParameters;
    PyObject* mCallbackObjective;
    PyObject* mCallbackGradient;
    PyObject* mCallbackHessian;
};
} // namespace NuTo
