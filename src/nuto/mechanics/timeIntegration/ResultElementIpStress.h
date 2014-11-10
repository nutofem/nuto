// $Id: $

#ifndef ResultElementIpStress_H
#define ResultElementIpStress_H

#include <ctime>
#include <array>

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/timeIntegration/ResultElementIpBase.h"
#include "nuto/base/ErrorEnum.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/FullVector.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all results
class ResultElementIpStress : public ResultElementIpBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief constructor
    ResultElementIpStress(const std::string& rIdent, int rElementId);

    //! @brief calculate the relevant nodal dofs
    void CalculateValues(const StructureBase& rStructure, NuTo::FullMatrix<double, 1, Eigen::Dynamic>& rValues)const override;

    //! @brief number of data points per time step (e.g. number of displacement components of a node
    int GetNumData(const StructureBase& rStructure)const;

    NuTo::TimeIntegration::eResultType GetResultType()const
    {
    	return NuTo::TimeIntegration::ELEMENT_IP_STRESS;
    }

    std::string GetTypeId() const
    {
    	return std::string("ResultElementIpStress");
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
};
}

//namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ResultElementIpStress)
#endif // SWIG
#endif  // ENABLE_SERIALIZATION
#endif // ResultElementIpStress_H
