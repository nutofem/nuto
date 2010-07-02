// $Id: ConstitutiveTangentNonlocal3x3.h 102 2009-11-11 10:47:23Z eckardt4 $

#ifndef CONSTITUTIVETANGENTNONLOCAL3x3_H
#define CONSTITUTIVETANGENTNONLOCAL3x3_H

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#else
#include <vector>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveTangentBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal3x3.h"
#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{
//! @brief ... tangent matrix for nonlocal local constitutive formulations
//! @author Jörg F. Unger, ISM
//! @date November 2009
class ConstitutiveTangentNonlocal3x3: public NuTo::ConstitutiveTangentBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class NonlocalDamagePlasticity;

public:
    //! @brief ... constructor
    ConstitutiveTangentNonlocal3x3(int rNumNonlocalElements);

    //! @brief ... destructor
    ~ConstitutiveTangentNonlocal3x3();

    //! @brief ... get the number of rows of the tangent matrix
    //! @return ... number of rows
    unsigned int GetNumberOfRows() const;

    //! @brief ... get the number of columns of the tangent matrix
    //! @return ... number of columns
    unsigned int GetNumberOfColumns() const;

    //! @brief ... get the tangent matrix
    //! @brief ... pointer to the tangent matrix (column major storage)
    const double* GetData() const;

    //! @brief ... get a local submatrix
    //! @param ... rSubMatrix number of the submatrix
    //! @return ... pointer to the tangent submatrix 3x3 matrix (column major storage)
    ConstitutiveTangentLocal3x3* GetSubMatrix(int rSubMatrix);

    //! @brief ... get a local submatrix
    //! @param ... rSubMatrix number of the submatrix
    //! @return ... pointer to the tangent submatrix 3x3 matrix (column major storage)
    const ConstitutiveTangentLocal3x3* GetSubMatrix(int rSubMatrix)const;

    //! @brief ... get a the number of local submatrices
    //! @return ... number of submatrices
    int GetNumSubMatrices()const;

    //! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
    ConstitutiveTangentNonlocal3x3* AsConstitutiveTangentNonlocal3x3()
	{
         return this;
	}

    //! @brief reinterpret as ConstitutiveTangentLocal1x1, otherwise throw an exception
    ConstitutiveTangentLocal1x1* AsConstitutiveTangentLocal1x1()
	{
         throw MechanicsException("[NuTo::ConstitutiveTangentDynamic::AsConstitutiveTangentLocal1x1] data types can not be cast.");
	}

    //! @brief reinterpret as ConstitutiveTangentLocal3x3, otherwise throw an exception
    ConstitutiveTangentLocal3x3* AsConstitutiveTangentLocal3x3()
	{
         throw MechanicsException("[NuTo::ConstitutiveTangentDynamic::AsConstitutiveTangentLocal3x3] data types can not be cast.");
	}

    //! @brief reinterpret as ConstitutiveTangentLocal6x6, otherwise throw an exception
    ConstitutiveTangentLocal6x6* AsConstitutiveTangentLocal6x6()
	{
         throw MechanicsException("[NuTo::ConstitutiveTangentDynamic::AsConstitutiveTangentLocal6x6] data types can not be cast.");
	}

    //! @brief returns the mIsNonlocal information
    bool IsLocal()
	{
         return mIsLocal;
	}
#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveTangentBase)
           & BOOST_SERIALIZATION_NVP(mNonlocalMatrices)
           & BOOST_SERIALIZATION_NVP(mIsLocal);
    }
#endif // ENABLE_SERIALIZATION

private:
    //! @brief ... vector of local tangents for nonlocal material model
    std::vector<ConstitutiveTangentLocal3x3> mNonlocalMatrices;
    //! @brief ... in case the nonlocal material model is only local (e.g. initial linear elastic part, unloading etc., only in mNonlocalMatrices[0] the current tangent is stored
    bool mIsLocal;
};

}

#endif // CONSTITUTIVETANGENTNONLOCAL3x3_H
