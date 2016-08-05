#pragma once

#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataBase.h"

//! @brief Stores the scalar history variable of a constitutive relationship
//! @author Philip Huschke, BAM
//! @date July 2016
namespace NuTo
{

class ConstitutiveStaticDataHistoryVariableScalar : public ConstitutiveStaticDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
	//! @brief constructor
    ConstitutiveStaticDataHistoryVariableScalar();

    //! @brief copy constructor
    ConstitutiveStaticDataHistoryVariableScalar(ConstitutiveStaticDataHistoryVariableScalar const& rOther) = default;

    //! @brief clones (copies) the data
    ConstitutiveStaticDataHistoryVariableScalar* Clone() const override
    {
        return new ConstitutiveStaticDataHistoryVariableScalar(*this);
    }

    //! @brief assignment operator
    ConstitutiveStaticDataHistoryVariableScalar& operator= (ConstitutiveStaticDataHistoryVariableScalar const& rOther) = default;

    //! @brief check, if the static data is compatible with a given element and a given constitutive model
    virtual bool CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const;

    //! \brief Alternative cast
    ConstitutiveStaticDataHistoryVariableScalar* AsHistoryVariableScalar() override
    {
        return this;
    }


    //! \brief Alternative cast
    const ConstitutiveStaticDataHistoryVariableScalar* AsHistoryVariableScalar()const override
    {
        return this;
    }

    void SetHistoryVariable(const double rHistoryVariable)
    {
        mHistoryVariable = rHistoryVariable;
    }

    double GetHistoryVariable() const
    {
        return mHistoryVariable;
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION



protected:
    //! @brief History variable
    double mHistoryVariable;

};

}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveStaticDataHistoryVariableScalar)
#endif // ENABLE_SERIALIZATION
