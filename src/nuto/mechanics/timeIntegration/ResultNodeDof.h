// $Id: $

#ifndef ResultNodeDof_H
#define ResultNodeDof_H

#include <ctime>
#include <array>

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/timeIntegration/ResultBase.h"
#include "nuto/base/ErrorEnum.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/FullVector.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all results
class NodeBase;
class ResultNodeDof : public ResultBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:

    //! @brief constructor
    ResultNodeDof(const std::string& rIdent, const NodeBase* rNodePtr);

    //! @brief calculate the relevant nodal dofs and add to the internal routine
    void CalculateAndAddValues(int rTimeStepPlot);

    void WriteToFile(const std::string& rResultDir, int rTimeStepPlot)const override;

    //! @brief resize the data (to avoid reallocation every time step, a certain number of
    //! @param rNumTimeSteps ... number of time steps to allocate the matrix for
    //! @param rInitValues ... true set all the values to zero, false - leave the values as is
    void Resize(int rNumTimeSteps, bool rInitValues);

    //! @brief calculate the relevant nodal dofs
    virtual void CalculateValues(NuTo::FullMatrix<double, 1, Eigen::Dynamic>& rValues)=0;

    //! @brief number of dofs (e.g. number of displacement components of a node
    virtual int GetNumDofs()const=0;

    NuTo::TimeIntegration::eResultType GetResultType()const
    {
    	return NuTo::TimeIntegration::NODE_DOF;
    }

    ResultNodeDof* AsResultNodeDof()
    {
    	return this;
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info()const;

protected:
    const NodeBase* mNodePtr;
    FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> mData;
};
}

//namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ResultNodeDof)
#endif // SWIG
#endif  // ENABLE_SERIALIZATION
#endif // ResultNodeDof_H
