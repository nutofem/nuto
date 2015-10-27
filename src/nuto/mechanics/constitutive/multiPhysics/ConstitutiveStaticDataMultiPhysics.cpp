#include "ConstitutiveStaticDataMultiPhysics.h"

//! @brief ... base constructor
NuTo::ConstitutiveStaticDataMultiPhysics::ConstitutiveStaticDataMultiPhysics()
    : ConstitutiveStaticDataBase()
{}

//! @brief ... copy constructor
NuTo::ConstitutiveStaticDataMultiPhysics::ConstitutiveStaticDataMultiPhysics(const NuTo::ConstitutiveStaticDataMultiPhysics &rOther)
    : ConstitutiveStaticDataBase(rOther)
{
    for (auto StaticData_it : rOther.mStaticData)
    {
        ConstitutiveStaticDataBase* ConstStaticData = StaticData_it.second->Clone();
        mStaticData.insert(std::pair<Constitutive::eConstitutiveStaticDataType,ConstitutiveStaticDataBase*>(StaticData_it.first,ConstStaticData));
    }
}


//! @brief ... destructor
NuTo::ConstitutiveStaticDataMultiPhysics::~ConstitutiveStaticDataMultiPhysics()
{
    for (auto StaticData_it : mStaticData)
    {
        delete StaticData_it.second;
    }
}


//! @brief Adds new static data to the multi physics static data
//! @param rStaticDataType Type of the static data that should be added
void NuTo::ConstitutiveStaticDataMultiPhysics::AddNewStaticData(NuTo::Constitutive::eConstitutiveStaticDataType rStaticDataType)
{
    switch (rStaticDataType)
    {
    case Constitutive::eConstitutiveStaticDataType::MOISTURE_TRANSPORT:
    {
        mStaticData.insert(std::pair<Constitutive::eConstitutiveStaticDataType,ConstitutiveStaticDataBase*>(rStaticDataType,new ConstitutiveStaticDataMoistureTransport));
        break;
    }
    default:
    {
        throw MechanicsException("[NuTo::ConstitutiveStaticDataMultiPhysics::AddNewStaticData] Unknown static data type.");
    }
    }
}



//!@ brief reinterpret as moisture transport
NuTo::ConstitutiveStaticDataMoistureTransport *NuTo::ConstitutiveStaticDataMultiPhysics::AsMoistureTransport()
{
    auto ConstLawStaticData_it = mStaticData.find(Constitutive::eConstitutiveStaticDataType::MOISTURE_TRANSPORT);
    if (ConstLawStaticData_it == mStaticData.end())
    {
        throw MechanicsException("[NuTo::ConstitutiveStaticDataMultiPhysics::AsMoistureTransport] There is no constitutive static data for moisture transport in the multi physics static data.");
    }
    return ConstLawStaticData_it->second->AsMoistureTransport();
}

//!@ brief reinterpret as moisture transport
const NuTo::ConstitutiveStaticDataMoistureTransport *NuTo::ConstitutiveStaticDataMultiPhysics::AsMoistureTransport() const
{
    auto ConstLawStaticData_it = mStaticData.find(Constitutive::eConstitutiveStaticDataType::MOISTURE_TRANSPORT);
    if (ConstLawStaticData_it == mStaticData.end())
    {
        throw MechanicsException("[NuTo::ConstitutiveStaticDataMultiPhysics::AsMoistureTransport] There is no constitutive static data for moisture transport in the multi physics static data.");
    }
    return ConstLawStaticData_it->second->AsMoistureTransport();
}


//!@ brief reinterpret as multi physics
NuTo::ConstitutiveStaticDataMultiPhysics *NuTo::ConstitutiveStaticDataMultiPhysics::AsMultiPhysics()
{
    return this;
}


//!@ brief reinterpret as multi physics
const NuTo::ConstitutiveStaticDataMultiPhysics *NuTo::ConstitutiveStaticDataMultiPhysics::AsMultiPhysics() const
{
    return this;
}

//! @brief clones (copies) the data
NuTo::ConstitutiveStaticDataMultiPhysics *NuTo::ConstitutiveStaticDataMultiPhysics::Clone() const
{
    return new ConstitutiveStaticDataMultiPhysics(*this);
}

//! @brief check, if the static data is compatible with a given element and a given constitutive model
bool NuTo::ConstitutiveStaticDataMultiPhysics::CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType) const
{
    if (rConstitutiveType==NuTo::Constitutive::MULTI_PHYSICS)
    {
        if (rElementType==NuTo::Element::ELEMENT2D)
            return true;
        else
            return false;
    }
    else
        return false;
}

