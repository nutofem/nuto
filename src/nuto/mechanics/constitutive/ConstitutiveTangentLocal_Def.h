// $Id: $

#ifndef CONSTITUTIVETANGENTLOCAL_DEF_H
#define CONSTITUTIVETANGENTLOCAL_DEF_H

#ifdef ENABLE_SERIALIZATION
// serialize
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/archive_exception.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveTangentBase.h"
#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{
//! @brief ... tangent matrix for local constitutive formulations
//! @author Jörg F. Unger, BAM
//! @date July 2012
template <int TNumRows, int TNumColumns>
class ConstitutiveTangentLocal: public NuTo::ConstitutiveTangentBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class LinearElasticEngineeringStress;
    friend class MisesPlasticityEngineeringStress;
    friend class NonlocalDamagePlasticityEngineeringStress;
    friend class LinearHeatFlux;
    friend class GradientDamagePlasticityEngineeringStress;

public:
    //! @brief ... constructor
    ConstitutiveTangentLocal();

    //! @brief ... destructor
    ~ConstitutiveTangentLocal();

    //! @brief ... get the number of rows of the tangent matrix
    //! @return ... number of rows
    unsigned int GetNumberOfRows() const;

    //! @brief ... get the number of columns of the tangent matrix
    //! @return ... number of columns
    unsigned int GetNumberOfColumns() const;

    //! @brief ... get the tangent matrix
    //! @brief ... pointer to the tangent matrix (column major storage)
    const double* GetData() const;

    //! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
    NuTo::ConstitutiveTangentLocal<1,1>& AsConstitutiveTangentLocal_1x1() override;

    //! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
    NuTo::ConstitutiveTangentLocal<2,1>& AsConstitutiveTangentLocal_2x1() override;

    //! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
    NuTo::ConstitutiveTangentLocal<2,2>& AsConstitutiveTangentLocal_2x2() override;

    //! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
    NuTo::ConstitutiveTangentLocal<3,3>& AsConstitutiveTangentLocal_3x3() override;

    //! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
    NuTo::ConstitutiveTangentLocal<6,1>& AsConstitutiveTangentLocal_6x1() override;

    //! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
    NuTo::ConstitutiveTangentLocal<6,6>& AsConstitutiveTangentLocal_6x6() override;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template <class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
    #ifdef DEBUG_SERIALIZATION
        std::cout << "start serialize ConstitutiveTangentLocal" << std::endl;
    #endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveTangentBase)
           & BOOST_SERIALIZATION_NVP(mTangent);
    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialize ConstitutiveTangentLocal" << std::endl;
    #endif
    }

#endif // ENABLE_SERIALIZATION
private:
     //! @brief ... tangent matrix
    double mTangent[TNumRows*TNumColumns];
};

}

//#ifdef ENABLE_SERIALIZATION
//BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveTangentLocal)
//#endif // ENABLE_SERIALIZATION

#endif // CONSTITUTIVETANGENTLOCAL_DEF_H
