// $Id $
#pragma once

#include "mechanics/elements/ElementOutputBase.h"

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

    ElementOutputFullMatrixDouble* Clone() const override
    {
    	return new ElementOutputFullMatrixDouble(*this);
    }

    Eigen::MatrixXd& GetFullMatrixDouble() override
	{
        return mMatrix;
	}
    void SetSymmetry(bool rSymmetric) override
	{
        mSymmetric = rSymmetric;
	}

    bool GetSymmetry()const override
	{
        return mSymmetric;
	}

    void SetConstant(bool rConstant) override
	{
    	mConstant = rConstant;
	}

    bool GetConstant()const override
	{
        return mConstant;
	}

private:
    bool mConstant;
    bool mSymmetric;
    Eigen::MatrixXd mMatrix;
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ElementOutputFullMatrixDouble)
#endif // ENABLE_SERIALIZATION
