#pragma once

#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataBase.h"

#include <map>
#include <memory>

namespace NuTo
{

//! @brief ... static data for combined constitutive laws
//! @author Volker Hirthammer, BAM
//! @date April 2016
class ConstitutiveStaticDataMultipleConstitutiveLaws : public ConstitutiveStaticDataBase
{
public:
    //! @brief constructor
    //! @param rConstitutiveLaws ... vector with constitutive laws which might need static data
    //! @param rElement ... related element
    //! @param rDim ... dimension of the static data (1D, 2D or 3D)
    ConstitutiveStaticDataMultipleConstitutiveLaws(const std::vector<NuTo::ConstitutiveBase*>& rConstitutiveLaws,
                                                   const ElementBase* rElement,
                                                   unsigned int rDim);

    //! @brief ... copy constructor
    ConstitutiveStaticDataMultipleConstitutiveLaws(const ConstitutiveStaticDataMultipleConstitutiveLaws& rOther);

    //! @brief ... move constructor
    ConstitutiveStaticDataMultipleConstitutiveLaws(ConstitutiveStaticDataMultipleConstitutiveLaws&& rOther) = delete;

    //! @brief ... copy assignment operator
    ConstitutiveStaticDataMultipleConstitutiveLaws& operator=(const ConstitutiveStaticDataMultipleConstitutiveLaws& rOther) = delete;

    //! @brief ... move assignment operator
    ConstitutiveStaticDataMultipleConstitutiveLaws& operator=(ConstitutiveStaticDataMultipleConstitutiveLaws&& rOther) = delete;

    //! @brief destructor
    ~ConstitutiveStaticDataMultipleConstitutiveLaws();








    //! @brief Allocates static data (if neccessary) for all the constitutive laws in the given vector
    //! @param rConstitutiveLaws ... vector with constitutive laws which might need static data
    //! @param rElement ... related element
    void AllocateStaticData1D(const std::vector<NuTo::ConstitutiveBase*>& rConstitutiveLaws,
                              const ElementBase* rElement);

    //! @brief Allocates static data (if neccessary) for all the constitutive laws in the given vector
    //! @param rConstitutiveLaws ... vector with constitutive laws which might need static data
    //! @param rElement ... related element
    void AllocateStaticData2D(const std::vector<NuTo::ConstitutiveBase*>& rConstitutiveLaws,
                              const ElementBase* rElement);
    //! @brief Allocates static data (if neccessary) for all the constitutive laws in the given vector
    //! @param rConstitutiveLaws ... vector with constitutive laws which might need static data
    //! @param rElement ... related element
    void AllocateStaticData3D(const std::vector<NuTo::ConstitutiveBase*>& rConstitutiveLaws,
                              const ElementBase* rElement);

    //!@ brief reinterpret as moisture transport
    virtual ConstitutiveStaticDataMoistureTransport* AsMoistureTransport();

    //!@ brief reinterpret as moisture transport
    virtual const ConstitutiveStaticDataMoistureTransport* AsMoistureTransport()const;

    //!@ brief reinterpret as multi physics
    virtual ConstitutiveStaticDataMultipleConstitutiveLaws* AsMultipleConstitutiveLaws() override
    {
        return this;
    }

    //!@ brief reinterpret as multi physics
    virtual const ConstitutiveStaticDataMultipleConstitutiveLaws *AsMultipleConstitutiveLaws()const override
    {
        return this;
    }

    //! @brief check, if the static data is compatible with a given element and a given constitutive model
    virtual bool CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType,
                                                NuTo::Element::eElementType rElementType) const override;


    //! @brief clones (copies) the data
    ConstitutiveStaticDataMultipleConstitutiveLaws* Clone()const override;


private:
    std::map<Constitutive::eConstitutiveStaticDataType,ConstitutiveStaticDataBase*>mStaticData;
};

}
