#include "ConstitutiveStaticDataMultipleConstitutiveLaws.h"







NuTo::ConstitutiveStaticDataMultipleConstitutiveLaws::ConstitutiveStaticDataMultipleConstitutiveLaws(const std::vector<NuTo::ConstitutiveBase *> &rConstitutiveLaws,
                                                                                                     const NuTo::ElementBase *rElement,
                                                                                                     unsigned int rDim)
    : ConstitutiveStaticDataBase()
{
    switch(rDim)
    {
    case 1:
        AllocateStaticData1D(rConstitutiveLaws,rElement);
        break;

    case 2:
        AllocateStaticData2D(rConstitutiveLaws,rElement);
        break;

    case 3:
        AllocateStaticData3D(rConstitutiveLaws,rElement);
        break;

    default:
        throw MechanicsException(__PRETTY_FUNCTION__,"Invalid diemnsion ("+std::to_string(rDim)+")");
    }
}


NuTo::ConstitutiveStaticDataMultipleConstitutiveLaws::ConstitutiveStaticDataMultipleConstitutiveLaws(const NuTo::ConstitutiveStaticDataMultipleConstitutiveLaws &rOther)
{
    for (auto itStaticData : rOther.mStaticData)
    {
        ConstitutiveStaticDataBase* constStaticData = itStaticData.second->Clone();
        mStaticData.insert(std::pair<Constitutive::eConstitutiveStaticDataType,ConstitutiveStaticDataBase*>(itStaticData.first,constStaticData));
    }
}




NuTo::ConstitutiveStaticDataMultipleConstitutiveLaws::~ConstitutiveStaticDataMultipleConstitutiveLaws()
{
    for (auto itStaticData : mStaticData)
    {
        delete itStaticData.second;
    }
}




void NuTo::ConstitutiveStaticDataMultipleConstitutiveLaws::AllocateStaticData1D(const std::vector<NuTo::ConstitutiveBase *> &rConstitutiveLaws,
                                                                                const NuTo::ElementBase *rElement)
{
//VHIRTHAMTODO check for nullptr instead
    for(unsigned int i=0; i<rConstitutiveLaws.size(); ++i)
    {
        switch(rConstitutiveLaws[i]->GetType())
        {
        case Constitutive::MOISTURE_TRANSPORT:
            mStaticData.emplace(Constitutive::eConstitutiveStaticDataType::MOISTURE_TRANSPORT,rConstitutiveLaws[i]->AllocateStaticData1D(rElement));
            break;

        // Do nothing
        case Constitutive::LINEAR_ELASTIC_ENGINEERING_STRESS:
        case Constitutive::SHRINKAGE_CAPILLARY_STRESS_BASED:
            break;
        default:
            throw MechanicsException(__PRETTY_FUNCTION__,std::string("Behaviour for constitutive law ")+
                                     Constitutive::ConstitutiveTypeToString(rConstitutiveLaws[i]->GetType())+" not specified");
        }
    }
}

void NuTo::ConstitutiveStaticDataMultipleConstitutiveLaws::AllocateStaticData2D(const std::vector<NuTo::ConstitutiveBase *> &rConstitutiveLaws,
                                                                                const NuTo::ElementBase *rElement)
{
    for(unsigned int i=0; i<rConstitutiveLaws.size(); ++i)
    {
        switch(rConstitutiveLaws[i]->GetType())
        {
        case Constitutive::MOISTURE_TRANSPORT:
            mStaticData.emplace(Constitutive::eConstitutiveStaticDataType::MOISTURE_TRANSPORT,rConstitutiveLaws[i]->AllocateStaticData2D(rElement));
            break;

        // Do nothing
        case Constitutive::LINEAR_ELASTIC_ENGINEERING_STRESS:
        case Constitutive::SHRINKAGE_CAPILLARY_STRESS_BASED:
            break;

        default:
            throw MechanicsException(__PRETTY_FUNCTION__,std::string("Behaviour for constitutive law ")+
                                     Constitutive::ConstitutiveTypeToString(rConstitutiveLaws[i]->GetType())+" not specified");
        }
    }
}

void NuTo::ConstitutiveStaticDataMultipleConstitutiveLaws::AllocateStaticData3D(const std::vector<NuTo::ConstitutiveBase *> &rConstitutiveLaws,
                                                                                const NuTo::ElementBase *rElement)
{
    for(unsigned int i=0; i<rConstitutiveLaws.size(); ++i)
    {
        switch(rConstitutiveLaws[i]->GetType())
        {
        case Constitutive::MOISTURE_TRANSPORT:
            mStaticData.emplace(Constitutive::eConstitutiveStaticDataType::MOISTURE_TRANSPORT,rConstitutiveLaws[i]->AllocateStaticData3D(rElement));
            break;
        // Do nothing
        case Constitutive::LINEAR_ELASTIC_ENGINEERING_STRESS:
        case Constitutive::SHRINKAGE_CAPILLARY_STRESS_BASED:
            break;
        default:
            throw MechanicsException(__PRETTY_FUNCTION__,std::string("Behaviour for constitutive law ")+
                                     Constitutive::ConstitutiveTypeToString(rConstitutiveLaws[i]->GetType())+" not specified");
        }
    }
}


NuTo::ConstitutiveStaticDataMoistureTransport *NuTo::ConstitutiveStaticDataMultipleConstitutiveLaws::AsMoistureTransport()
{
    auto itConstLawStaticData = mStaticData.find(Constitutive::eConstitutiveStaticDataType::MOISTURE_TRANSPORT);
    if (itConstLawStaticData == mStaticData.end())
    {
        throw MechanicsException(__PRETTY_FUNCTION__,"No constitutive static data for moisture transport found!");
    }
    return itConstLawStaticData->second->AsMoistureTransport();
}


const NuTo::ConstitutiveStaticDataMoistureTransport *NuTo::ConstitutiveStaticDataMultipleConstitutiveLaws::AsMoistureTransport() const
{
    auto itConstLawStaticData = mStaticData.find(Constitutive::eConstitutiveStaticDataType::MOISTURE_TRANSPORT);
    if (itConstLawStaticData == mStaticData.end())
    {
        throw MechanicsException(__PRETTY_FUNCTION__,"No constitutive static data for moisture transport found!");
    }
    return itConstLawStaticData->second->AsMoistureTransport();
}


bool NuTo::ConstitutiveStaticDataMultipleConstitutiveLaws::CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType,
                                                                                          NuTo::Element::eElementType rElementType) const
{
    switch(rConstitutiveType)
    {
    case Constitutive::CONSTITUTIVE_LAWS_ADDITIVE_OUTPUT:
        break;
    default:
        return false;
    }

    switch(rElementType)
    {
    case Element::CONTINUUMELEMENT:
        break;
    default:
        return false;
    }

    return true;
}

NuTo::ConstitutiveStaticDataMultipleConstitutiveLaws *NuTo::ConstitutiveStaticDataMultipleConstitutiveLaws::Clone() const
{
    return new ConstitutiveStaticDataMultipleConstitutiveLaws(*this);
}
