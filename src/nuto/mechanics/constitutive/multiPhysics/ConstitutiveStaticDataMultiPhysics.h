#ifndef CONSTITUTIVESTATICDATAMULTIPHYSICS_H
#define CONSTITUTIVESTATICDATAMULTIPHYSICS_H

#include "nuto/mechanics/constitutive/moistureTransport/ConstitutiveStaticDataMoistureTransport.h"
#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataMoistureTransport.h"

#include <map>
//#include <memory>

namespace NuTo
{

class ConstitutiveStaticDataMultiPhysics : public ConstitutiveStaticDataBase
{
public:
    //! @brief ... base constructor
    ConstitutiveStaticDataMultiPhysics();

    //! @brief ... copy constructor
    ConstitutiveStaticDataMultiPhysics(const ConstitutiveStaticDataMultiPhysics& rOther);

    //! @brief ... move constructor
    ConstitutiveStaticDataMultiPhysics(ConstitutiveStaticDataMultiPhysics&& rOther) = delete;

    //! @brief ... copy assignment operator
    ConstitutiveStaticDataMultiPhysics& operator=(const ConstitutiveStaticDataMultiPhysics& rOther) = delete;

    //! @brief ... move assignment operator
    ConstitutiveStaticDataMultiPhysics& operator=(ConstitutiveStaticDataMultiPhysics&& rOther) = delete;

    //! @brief ... destructor
    ~ConstitutiveStaticDataMultiPhysics();

    //! @brief Adds new static data to the multi physics static data
    //! @param rStaticDataType Type of the static data that should be added
    void AddNewStaticData(Constitutive::eConstitutiveStaticDataType rStaticDataType);

    //!@ brief reinterpret as moisture transport
    virtual ConstitutiveStaticDataMoistureTransport* AsMoistureTransport();

    //!@ brief reinterpret as moisture transport
    virtual const ConstitutiveStaticDataMoistureTransport* AsMoistureTransport()const;

    //!@ brief reinterpret as multi physics
    virtual ConstitutiveStaticDataMultiPhysics* AsMultiPhysics() override;

    //!@ brief reinterpret as multi physics
    virtual const ConstitutiveStaticDataMultiPhysics *AsMultiPhysics()const override;


    //! @brief clones (copies) the data
    ConstitutiveStaticDataMultiPhysics*                     Clone                                   ()const override;

    //! @brief check, if the static data is compatible with a given element and a given constitutive model
    virtual bool                                            CheckConstitutiveCompatibility          (NuTo::Constitutive::eConstitutiveType rConstitutiveType,
                                                                                                     NuTo::Element::eElementType rElementType) const override;

private:

    std::map<Constitutive::eConstitutiveStaticDataType,ConstitutiveStaticDataBase*>mStaticData;
};

} // namespace NuTo

#endif // CONSTITUTIVESTATICDATAMULTIPHYSICS_H
