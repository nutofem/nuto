// $Id: ConstitutiveTangentLocal6x6.h 102 2009-11-11 10:47:23Z eckardt4 $

#ifndef CONSTITUTIVETANGENTLOCAL_6x6_H
#define CONSTITUTIVETANGENTLOCAL_6x6_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveTangentBase.h"
#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{
class LinearElastic;
class LinearElasticEngineeringStress;
class ConstitutiveMisesPlasticity;
//! @brief ... tangent matrix for local constitutive formulations
//! @author Jörg F. Unger, ISM
//! @date November 2009
class ConstitutiveTangentLocal6x6: public NuTo::ConstitutiveTangentBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class LinearElastic;
    friend class LinearElasticEngineeringStress;
    friend class ConstitutiveMisesPlasticity;
    friend class NonlocalDamagePlasticity;

public:
    //! @brief ... constructor
    ConstitutiveTangentLocal6x6();

    //! @brief ... destructor
    ~ConstitutiveTangentLocal6x6();

    //! @brief ... get the number of rows of the tangent matrix
    //! @return ... number of rows
    unsigned int GetNumberOfRows() const;

    //! @brief ... get the number of columns of the tangent matrix
    //! @return ... number of columns
    unsigned int GetNumberOfColumns() const;

    //! @brief ... get the tangent matrix
    //! @brief ... pointer to the tangent matrix (column major storage)
    const double* GetData() const;

    ConstitutiveTangentLocal6x6& GetConstitutiveTangentLocal6x6()
    {
    	return *this;
    }

    //! @brief reinterpret as ConstitutiveTangentNonlocal3x3, otherwise throw an exception
    ConstitutiveTangentNonlocal3x3* AsConstitutiveTangentNonlocal3x3()
	{
         throw MechanicsException("[NuTo::ConstitutiveTangentLocal6x6::AsConstitutiveTangentNonlocal3x3] data types can not be cast.");
	}

    //! @brief reinterpret as ConstitutiveTangentLocal1x1, otherwise throw an exception
    ConstitutiveTangentLocal1x1* AsConstitutiveTangentLocal1x1()
	{
        throw MechanicsException("[NuTo::ConstitutiveTangentLocal6x6::AsConstitutiveTangentLocal1x1] data types can not be cast.");
 	}

    //! @brief reinterpret as ConstitutiveTangentLocal3x3, otherwise throw an exception
    ConstitutiveTangentLocal2x2* AsConstitutiveTangentLocal2x2()
	{
        throw MechanicsException("[NuTo::ConstitutiveTangentLocal3x3::AsConstitutiveTangentLocal1x1] data types can not be cast.");
	}

    //! @brief reinterpret as ConstitutiveTangentLocal3x3, otherwise throw an exception
    ConstitutiveTangentLocal3x3* AsConstitutiveTangentLocal3x3()
	{
        throw MechanicsException("[NuTo::ConstitutiveTangentLocal6x6::AsConstitutiveTangentLocal3x3] data types can not be cast.");
	}

    //! @brief reinterpret as ConstitutiveTangentLocal6x6, otherwise throw an exception
    ConstitutiveTangentLocal6x6* AsConstitutiveTangentLocal6x6()
	{
        return this;
	}
#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
private:
    //! @brief ... tangent matrix
    double mTangent[36];
};

}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveTangentLocal6x6)
#endif // ENABLE_SERIALIZATION

#endif // CONSTITUTIVETANGENTLOCAL_6x6_H
