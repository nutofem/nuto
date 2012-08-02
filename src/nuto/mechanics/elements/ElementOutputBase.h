// $Id $
#ifndef ELEMENT_OUTPUTBASE_H_
#define ELEMENT_OUTPUTBASE_H_

#include <vector>
#include "nuto/mechanics/elements/IpDataEnum.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 29, 2010
//! @brief ...
class ElementOutputBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    ElementOutputBase(){};

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    virtual FullMatrix<double>& GetFullMatrixDouble()
	{
        throw MechanicsException("[ElementOutputBase::GetFullMatrixDouble] element output matrix is not of type FullMatrix<double>");
	}

    virtual FullMatrix<int>& GetFullMatrixInt()
	{
        throw MechanicsException("[ElementOutputBase::GetFullMatrixDouble] element output matrix is not of type FullMatrix<double>");
	}
    virtual std::vector<int>& GetVectorInt()
	{
        throw MechanicsException("[ElementOutputBase::GetFullMatrixDouble] element output matrix is not of type std::vector<int>");
	}

    virtual NuTo::IpData::eIpStaticDataType GetIpDataType()
	{
        throw MechanicsException("[ElementOutputBase::GetIpDataType] ipdata is not stored.");
	}

    virtual void SetSymmetry(bool rSymmetric)
	{
        throw MechanicsException("[ElementOutputBase::SetSymmetry] symmetry is not stored.");
	}

    virtual bool GetSymmetry()const
	{
        throw MechanicsException("[ElementOutputBase::SetSymmetry] symmetry is not stored.");
	}

    virtual void SetConstant(bool rConstant)
	{
        throw MechanicsException("[ElementOutputBase::SetConstant] constness is not stored.");
	}

    virtual bool GetConstant()const
	{
        throw MechanicsException("[ElementOutputBase::GetConstant] constness is not stored.");
	}
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ElementOutputBase)
#endif // ENABLE_SERIALIZATION
#endif /* ELEMENT_OUTPUTBASE_H_ */
