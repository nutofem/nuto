// $Id: ConstitutiveTangentBase.h 89 2009-11-06 16:19:14Z eckardt4 $

#ifndef CONSTITUTIVETANGENTBASE_H_
#define CONSTITUTIVETANGENTBASE_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

namespace NuTo
{
//! @brief ... base class storing the tangent of the constitutive relationship
//! @author Jörg F. Unger, ISM
//! @date November 2009
class ConstitutiveTangentLocalBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief ... constructor
    ConstitutiveTangentLocalBase();

    //! @brief ... get the number of rows of the tangent matrix
    //! @return ... number of rows
    virtual unsigned int GetNumberOfRows() const = 0;

    //! @brief ... get the number of columns of the tangent matrix
    //! @return ... number of columns
    virtual unsigned int GetNumberOfColumns() const = 0;

    //! @brief ... returns the shape of the tangent matrix (symmetric or nonsymmetric)
    //! @return ... <B>true</B> if the matrix is symmetric / <B>false</B> if the matrix is not symmetric
    //! @sa mSymmetry
    inline bool GetSymmetry() const
    {
        return this->mSymmetry;
    }

    //! @brief ... get the tangent matrix
    //! @brief ... pointer to the tangent matrix (column major storage)
    virtual const double* GetData() const = 0;

    //! @brief ... set shape of the tangent matrix (symmetric or nonsymmetric)
    //! @param rSymmetry ... symmetry flag: <B>true</B> if the tangent matrix has a symmetric shape and <B>false</B> if the tangent is nonsymmetric.
    //! @sa mSymmetry
    inline void SetSymmetry(bool rSymmetry)
    {
        this->mSymmetry = rSymmetry;
    }

    //! @brief ... set tangent matrix
    //! @param rNumberOfRows ... number of rows
    //! @param rNumberOfColumns ... number of columns
    //! @param rTangentMatrix ... tangent matrix (column major storage)
    virtual void SetData(const double* rTangentMatrix) = 0;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(this->mSymmetry);
    }
#endif // ENABLE_SERIALIZATION
private:
    //! @brief ... determines the shape of the tangent matrix (symmetric or nonsymmetric)
    bool mSymmetry;
};

}

#endif // CONSTITUTIVETANGENTBASE_H_ 
