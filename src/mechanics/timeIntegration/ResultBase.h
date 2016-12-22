// $Id: $

#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION


#include "math/FullMatrix_Def.h"
#include "mechanics/MechanicsException.h"

namespace NuTo
{
class ResultNodeDof;
class ResultGroupNodeDof;
class ResultTime;
class ResultElementIpData;
class StructureBase;
enum class eTimeIntegrationResultType;

//! @author Jörg F. Unger, BAM
//! @date December 2013
//! @brief ... standard abstract class for all results
class ResultBase
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

    void Resize(const StructureBase& rStructure, int rNumResultSteps, bool rInitialize);

    void WriteToFile(const std::string& rResultDir, int rNumTimeSteps)const;

    //! @brief number of data points per time step (e.g. number of displacement components of a node)
    virtual int GetNumData(const StructureBase& rStructure)const=0;

    virtual NuTo::eTimeIntegrationResultType GetResultType()const = 0;

    virtual ResultNodeDof* AsResultNodeDof()
    {
    	throw MechanicsException(std::string(__PRETTY_FUNCTION__) +"\t: object is not of this type.");
    }

    virtual ResultTime* AsResultTime()
    {
    	throw MechanicsException(std::string(__PRETTY_FUNCTION__) +"\t: object is not of this type.");
    }

    virtual ResultGroupNodeDof* AsResultGroupNodeDof()
    {
    	throw MechanicsException(std::string(__PRETTY_FUNCTION__) +"\t: object is not of this type.");
    }


    virtual ResultElementIpData* AsResultElementIpData()
    {
        throw MechanicsException(std::string(__PRETTY_FUNCTION__) +"\t: object is not of this type.");
    }

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
    FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> mData;
};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ResultBase)
#endif // SWIG
#endif  // ENABLE_SERIALIZATION
