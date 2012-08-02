// $Id: ConstitutiveTangentLocal3x3.h 102 2009-11-11 10:47:23Z eckardt4 $

#ifndef CONSTITUTIVETANGENTLOCAL_3x3_H
#define CONSTITUTIVETANGENTLOCAL_3x3_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveTangentBase.h"
#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{
class LinearElastic;
class ConstitutiveMisesPlasticity;
class ConstitutiveLatticeConcrete;
class LinearHeatFlux;
template <class T> class FullMatrix;
//! @brief ... tangent matrix for local constitutive formulations
//! @author Jörg F. Unger, ISM
//! @date November 2009
class ConstitutiveTangentLocal3x3: public NuTo::ConstitutiveTangentBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class LinearElastic;
    friend class ConstitutiveMisesPlasticity;
    friend class Multiscale;
    friend class NonlocalDamagePlasticity;
    friend class ConstitutiveLatticeConcrete;
    friend class LinearHeatFlux;

public:
    //! @brief ... constructor
    ConstitutiveTangentLocal3x3();

    //! @brief ... copy constructor from matrix
    ConstitutiveTangentLocal3x3& operator= (const NuTo::FullMatrix<double>& rOtherMatrix);

    //! @brief ... destructor
    ~ConstitutiveTangentLocal3x3();

    //! @brief ... get the number of rows of the tangent matrix
    //! @return ... number of rows
    unsigned int GetNumberOfRows() const;

    //! @brief ... get the number of columns of the tangent matrix
    //! @return ... number of columns
    unsigned int GetNumberOfColumns() const;

    //! @brief ... get the tangent matrix
    //! @brief ... pointer to the tangent matrix (column major storage)
    const double* GetData() const;

    NuTo::ConstitutiveTangentLocal3x3& GetConstitutiveTangentLocal3x3()
    {
    	return *this;
    }

    //! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
    ConstitutiveTangentNonlocal3x3* AsConstitutiveTangentNonlocal3x3()
	{
         throw MechanicsException("[NuTo::ConstitutiveTangentLocal3x3::AsConstitutiveTangentNonlocal3x3] data types can not be cast.");
	}

    //! @brief reinterpret as ConstitutiveTangentLocal1x1, otherwise throw an exception
    ConstitutiveTangentLocal1x1* AsConstitutiveTangentLocal1x1()
	{
        throw MechanicsException("[NuTo::ConstitutiveTangentLocal3x3::AsConstitutiveTangentLocal1x1] data types can not be cast.");
 	}

    //! @brief reinterpret as ConstitutiveTangentLocal3x3, otherwise throw an exception
    ConstitutiveTangentLocal2x2* AsConstitutiveTangentLocal2x2()
	{
        throw MechanicsException("[NuTo::ConstitutiveTangentLocal3x3::AsConstitutiveTangentLocal1x1] data types can not be cast.");
	}


    //! @brief reinterpret as ConstitutiveTangentLocal3x3, otherwise throw an exception
    ConstitutiveTangentLocal3x3* AsConstitutiveTangentLocal3x3()
	{
        return this;
	}

    //! @brief reinterpret as ConstitutiveTangentLocal6x6, otherwise throw an exception
    ConstitutiveTangentLocal6x6* AsConstitutiveTangentLocal6x6()
	{
         throw MechanicsException("[NuTo::ConstitutiveTangentLocal3x3::AsConstitutiveTangentLocal6x6] data types can not be cast.");
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
    double mTangent[9];
};

}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveTangentLocal3x3)
#endif // ENABLE_SERIALIZATION

#endif // CONSTITUTIVETANGENTLOCAL_3x3_H
