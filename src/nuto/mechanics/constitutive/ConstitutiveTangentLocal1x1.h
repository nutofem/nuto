// $Id: ConstitutiveTangentLocal.h 102 2009-11-11 10:47:23Z eckardt4 $

#ifndef CONSTITUTIVETANGENTLOCAL_1x1_H
#define CONSTITUTIVETANGENTLOCAL_1x1_H

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/array.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveTangentLocalBase.h"

namespace NuTo
{
class LinearElastic;
class ConstitutiveMisesPlasticity;
//! @brief ... tangent matrix for local constitutive formulations
//! @author Jörg F. Unger, ISM
//! @date November 2009
class ConstitutiveTangentLocal1x1: public NuTo::ConstitutiveTangentLocalBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class LinearElastic;
    friend class ConstitutiveMisesPlasticity;

public:
    //! @brief ... constructor
    ConstitutiveTangentLocal1x1();

    //! @brief ... destructor
    ~ConstitutiveTangentLocal1x1();

    //! @brief ... get the number of rows of the tangent matrix
    //! @return ... number of rows
    unsigned int GetNumberOfRows() const;

    //! @brief ... get the number of columns of the tangent matrix
    //! @return ... number of columns
    unsigned int GetNumberOfColumns() const;

    //! @brief ... get the tangent matrix
    //! @brief ... pointer to the tangent matrix (column major storage)
    const double* GetData() const;

    //! @brief ... set tangent matrix
    //! @param rTangentMatrix ... tangent matrix (column major storage)
    void SetData(const double* rTangentMatrix);

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveTangentLocalBase)
           & BOOST_SERIALIZATION_NVP(mTangent);
    }
#endif // ENABLE_SERIALIZATION
private:
    //! @brief ... tangent matrix
    double mTangent;
};

}

#endif // CONSTITUTIVETANGENTLOCAL_1x1_H
