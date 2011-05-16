#ifndef NUTO_MECHANICS_EXCEPTION_H
#define NUTO_MECHANICS_EXCEPTION_H

#include "nuto/base/Exception.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date September 2009
//! @brief ... class for all exceptions thrown in metamodel module of NuTo
class MechanicsException : public NuTo::Exception
{
public:
	enum eError
	{
	    NOTSPECIFIED=0,
	    NOCONVERGENCE                      //!< constitutive or structure is not converging
	};
    //! @brief ...constructor
    //! @param message_ ...error message
    explicit MechanicsException(const std::string& message_, eError rError=NOTSPECIFIED) : NuTo::Exception(message_)
    {
    	mError=rError;
    }

    //! @brief ...constructor
    //! @param message_ ...error message
    //! @param FatalFlag_ ...flag to decide, if the error is fatal or not (exit the program or able to continue in a consistent way)
    explicit MechanicsException(const std::string& message_, bool FatalFlag_, eError rError=NOTSPECIFIED) : NuTo::Exception(message_, FatalFlag_)
    {
    	mError=rError;
    }

    //! @brief ...constructor
    //! @param message_ ...error message
    //! @param FatalFlag_ ...flag to decide, if the error is fatal or not (exit the program or able to continue in a consistent way)
    explicit MechanicsException(const std::string& message_, bool FatalFlag_) : NuTo::Exception(message_, FatalFlag_)
    {}

    //! @brief ...destructor
    virtual ~MechanicsException() throw()
    {}

    //! @brief ... return error message of the exception
    //! @return ... error message
    virtual std::string ErrorMessage() const throw();

    eError GetError()const
    {
    	return mError;
    }

    //! @brief ... clone the exception (important, if called from the base class)
    //! @return ... a copy of the exception
    Exception* Clone();
private:
    eError mError;
};
} //namespace NuTo
#endif //NUTO_MECHANICS_EXCEPTION_H
