// $Id $
#ifndef ELEMENT_OUTPUT_FULLMATRIX_DOUBLE_H_
#define ELEMENT_OUTPUT_FULLMATRIX_DOUBLE_H_

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/elements/ElementOutputBase.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 29, 2010
//! @brief ...
class ElementOutputFullMatrixDouble : public ElementOutputBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    ElementOutputFullMatrixDouble()
    {
    	mSymmetric = false;
    	mConstant  = false;
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    ElementOutputFullMatrixDouble* Clone() const
    {
    	return new ElementOutputFullMatrixDouble(*this);
    }

    FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& GetFullMatrixDouble() override
	{
        return mMatrix;
	}
    void SetSymmetry(bool rSymmetric)
	{
        mSymmetric = rSymmetric;
	}

    bool GetSymmetry()const
	{
        return mSymmetric;
	}

    void SetConstant(bool rConstant)
	{
    	mConstant = rConstant;
	}

    bool GetConstant()const
	{
        return mConstant;
	}

private:
    bool mConstant;
    bool mSymmetric;
    FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> mMatrix;
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ElementOutputFullMatrixDouble)
#endif // ENABLE_SERIALIZATION
#endif /* ELEMENT_OUTPUT_FULLMATRIX_DOUBLE_H_ */
