// $Id: $

#ifndef ResultNodeDisp_H
#define ResultNodeDisp_H

#include <ctime>
#include <array>

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/timeIntegration/ResultNodeDof.h"
#include "nuto/base/ErrorEnum.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/FullVector.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all results
class ResultNodeDisp : public ResultNodeDof
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief constructor
    ResultNodeDisp(const std::string& rIdent, const NodeBase* rNodePtr);

    //! @brief calculate the relevant nodal dofs
    void CalculateValues(NuTo::FullMatrix<double, 1, Eigen::Dynamic>& rValues)override;

    //! @brief number of dofs (e.g. number of displacement components of a node
    int GetNumDofs()const;

    std::string GetTypeId() const
    {
    	return std::string("ResultNodeDisp");
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
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ResultNodeDisp)
#endif // SWIG
#endif  // ENABLE_SERIALIZATION
#endif // ResultNodeDisp_H
