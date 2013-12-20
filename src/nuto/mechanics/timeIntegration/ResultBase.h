// $Id: $

#ifndef RESULTBASE_H
#define RESULTBASE_H

#include <ctime>
#include <array>

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/base/NuToObject.h"
#include "nuto/base/ErrorEnum.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/FullVector.h"
#include "nuto/mechanics/timeIntegration/TimeIntegrationEnum.h"
namespace NuTo
{
class ResultNodeDof;
class ResultReactionNodeDofGroup;
class ResultTime;
//! @author Jörg F. Unger, BAM
//! @date December 2013
//! @brief ... standard abstract class for all results
class ResultBase : public NuToObject
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief constructor
    ResultBase(const std::string& rIdent);

    //! @brief deconstructor
    virtual ~ResultBase();

    void SetIdent(const std::string& rIdent);

    std::string GetIdent()const;

    bool IsCalculated() const;

    void SetCalculated(bool rCalculated);

    virtual void Resize(int rNumResultSteps, bool rInitialize)=0;

    virtual void WriteToFile(const std::string& rResultDir, int rNumTimeSteps)const=0;

    virtual NuTo::TimeIntegration::eResultType GetResultType()const = 0;

    virtual ResultNodeDof* AsResultNodeDof()
    {
    	throw MechanicsException("[NutO::ResultBase::AsResultNodeDof] object is not of this type.");
    }

    virtual ResultTime* AsResultTime()
    {
    	throw MechanicsException("[NutO::ResultBase::AsResultTime] object is not of this type.");
    }

/*    virtual ResultReactionNodeDofGroup* AsResultReactionNodeDofGroup()
    {
    	throw MechanicsException("[NutO::ResultBase::AsResultNodalDof] object is not of this type.");
    }
*/

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info()const
    {

    }

protected:
    std::string mIdent;
    bool mCalculated;
};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ResultBase)
#endif // SWIG
#endif  // ENABLE_SERIALIZATION
#endif // RESULTBASE_H
