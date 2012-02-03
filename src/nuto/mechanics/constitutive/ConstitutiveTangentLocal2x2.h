// $Id:$

#ifndef CONSTITUTIVETANGENTLOCAL_2x2_H
#define CONSTITUTIVETANGENTLOCAL_2x2_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveTangentBase.h"
#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{
class ConstitutiveLatticeConcrete;
template <class T> class FullMatrix;
//! @brief ... tangent matrix for local constitutive formulations (interface 2D)
//! @author Jörg F. Unger, ISM
//! @date Jan 2012
class ConstitutiveTangentLocal2x2: public NuTo::ConstitutiveTangentBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class ConstitutiveLatticeConcrete;
public:
    //! @brief ... constructor
    ConstitutiveTangentLocal2x2();

    //! @brief ... copy constructor from matrix
    ConstitutiveTangentLocal2x2& operator= (const NuTo::FullMatrix<double>& rOtherMatrix);

    //! @brief ... destructor
    ~ConstitutiveTangentLocal2x2();

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
    ConstitutiveTangentNonlocal3x3* AsConstitutiveTangentNonlocal3x3()
	{
         throw MechanicsException("[NuTo::ConstitutiveTangentLocal2x2::AsConstitutiveTangentNonlocal3x3] data types can not be cast.");
	}

    //! @brief reinterpret as ConstitutiveTangentLocal1x1, otherwise throw an exception
    ConstitutiveTangentLocal1x1* AsConstitutiveTangentLocal1x1()
	{
        throw MechanicsException("[NuTo::ConstitutiveTangentLocal2x2::AsConstitutiveTangentLocal1x1] data types can not be cast.");
 	}

    //! @brief reinterpret as ConstitutiveTangentLocal2x2, otherwise throw an exception
    ConstitutiveTangentLocal2x2* AsConstitutiveTangentLocal2x2()
	{
        return this;
	}

    //! @brief reinterpret as ConstitutiveTangentLocal3x3, otherwise throw an exception
    ConstitutiveTangentLocal3x3* AsConstitutiveTangentLocal3x3()
	{
        throw MechanicsException("[NuTo::ConstitutiveTangentLocal2x2::AsConstitutiveTangentLocal1x1] data types can not be cast.");
	}

    //! @brief reinterpret as ConstitutiveTangentLocal6x6, otherwise throw an exception
    ConstitutiveTangentLocal6x6* AsConstitutiveTangentLocal6x6()
	{
         throw MechanicsException("[NuTo::ConstitutiveTangentLocal2x2::AsConstitutiveTangentLocal6x6] data types can not be cast.");
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
    double mTangent[4];
};

}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveTangentLocal2x2)
#endif // ENABLE_SERIALIZATION

#endif // CONSTITUTIVETANGENTLOCAL_2x2_H
