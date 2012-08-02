// $Id: $

#ifndef CONSTITUTIVETANGENTNONLOCAL__DEF
#define CONSTITUTIVETANGENTNONLOCAL__DEF

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

#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"
#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{
//! @brief ... tangent matrix for local constitutive formulations
//! @author Jörg F. Unger, BAM
//! @date July 2012
template <int TNumRows, int TNumColumns>
class ConstitutiveTangentNonlocal: public NuTo::ConstitutiveTangentBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class LinearElasticEngineeringStress;
    friend class ConstitutiveMisesPlasticity;
    friend class NonlocalDamagePlasticity;

public:
    //! @brief ... constructor
    ConstitutiveTangentNonlocal(int rNumNonlocalMatrices=0)
    {
    	mMatrices.resize(rNumNonlocalMatrices);
    }

    //! @brief ... destructor
    ~ConstitutiveTangentNonlocal()
    {}

    //! @brief ... sets the number of local matrices
    //! @return ... n
    void SetNumSubMatrices(int rNumSubMatrices) const
    {
    	mMatrices.resize(rNumSubMatrices);
    }


    //! @brief ... get the number of rows of the tangent matrix
    //! @return ... number of rows
    unsigned int GetNumberOfRows() const
    {
    	return TNumRows;
    }

    //! @brief ... get the number of columns of the tangent matrix
    //! @return ... number of columns
    unsigned int GetNumberOfColumns() const
    {
    	return TNumColumns*mMatrices.size();
    }

    //! @brief ... get the tangent matrix
    //! @brief ... pointer to the tangent matrix (column major storage)
    const double* GetData() const
    {
    	throw MechanicsException("[ConstitutiveTangentNonlocal::GetData] not implemented for nonlocal matrices, get the data for each submatrix separately.");
    }

    //! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
    NuTo::ConstitutiveTangentLocal<1,1>& GetSubMatrix_1x1(int rSubMatrix) override
    {
        if (rSubMatrix<(int)mMatrices.size())
        	return mMatrices[rSubMatrix].AsConstitutiveTangentLocal_1x1();
        else
        	throw MechanicsException("[ConstitutiveTangentNonlocal::GetSubMatrix_1x1] number of nonlocal matrices is less then requested submatrix.");
    }

    //! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
    NuTo::ConstitutiveTangentLocal<2,2>& GetSubMatrix_2x2(int rSubMatrix) override
	{
        if (rSubMatrix<(int)mMatrices.size())
        	return mMatrices[rSubMatrix].AsConstitutiveTangentLocal_2x2();
        else
        	throw MechanicsException("[ConstitutiveTangentNonlocal::GetSubMatrix_2x2] number of nonlocal matrices is less then requested submatrix.");
	}

    //! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
    NuTo::ConstitutiveTangentLocal<3,3>& GetSubMatrix_3x3(int rSubMatrix) override
	{
        if (rSubMatrix<(int)mMatrices.size())
        	return mMatrices[rSubMatrix].AsConstitutiveTangentLocal_3x3();
        else
        	throw MechanicsException("[ConstitutiveTangentNonlocal::GetSubMatrix_3x3] number of nonlocal matrices is less then requested submatrix.");
	}

    //! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
    NuTo::ConstitutiveTangentLocal<6,6>& GetSubMatrix_6x1(int rSubMatrix) override
	{
		if (rSubMatrix<(int)mMatrices.size())
			return mMatrices[rSubMatrix].AsConstitutiveTangentLocal_6x1();
		else
			throw MechanicsException("[ConstitutiveTangentNonlocal::GetSubMatrix_6x1] number of nonlocal matrices is less then requested submatrix.");
	}


    //! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
    NuTo::ConstitutiveTangentLocal<6,1>& GetSubMatrix_6x6(int rSubMatrix) override
	{
        if (rSubMatrix<(int)mMatrices.size())
        	return mMatrices[rSubMatrix].AsConstitutiveTangentLocal_6x6();
        else
        	throw MechanicsException("[ConstitutiveTangentNonlocal::GetSubMatrix_6x6] number of nonlocal matrices is less then requested submatrix.");
	}

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
           & BOOST_SERIALIZATION_NVP(mMatrices);
    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialize ConstitutiveTangentLocal" << std::endl;
    #endif
    }

#endif // ENABLE_SERIALIZATION
private:
     //! @brief ... tangent matrices
    std::vector<ConstitutiveTangentLocal<TNumRows,TNumColumns> > mMatrices;
};

}

//#ifdef ENABLE_SERIALIZATION
//BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveTangentLocal)
//#endif // ENABLE_SERIALIZATION

#endif // CONSTITUTIVETANGENTNONLOCAL__DEF
