#include "ConstitutiveStaticDataMultipleConstitutiveLaws.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/elements/ElementEnum.h"






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
        mStaticData.insert(std::pair<Constitutive::eConstitutiveType,ConstitutiveStaticDataBase*>(itStaticData.first,constStaticData));
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
    for(unsigned int i=0; i<rConstitutiveLaws.size(); ++i)
    {
        ConstitutiveStaticDataBase* staticData = rConstitutiveLaws[i]->AllocateStaticData1D(rElement);
        if(staticData != nullptr)
        {
            mStaticData.emplace(rConstitutiveLaws[i]->GetType(), staticData);
        }
    }
}

void NuTo::ConstitutiveStaticDataMultipleConstitutiveLaws::AllocateStaticData2D(const std::vector<NuTo::ConstitutiveBase *> &rConstitutiveLaws,
                                                                                const NuTo::ElementBase *rElement)
{
    for(unsigned int i=0; i<rConstitutiveLaws.size(); ++i)
    {
        ConstitutiveStaticDataBase* staticData = rConstitutiveLaws[i]->AllocateStaticData2D(rElement);
        if(staticData != nullptr)
        {
            mStaticData.emplace(rConstitutiveLaws[i]->GetType(), staticData);
        }
    }
}

void NuTo::ConstitutiveStaticDataMultipleConstitutiveLaws::AllocateStaticData3D(const std::vector<NuTo::ConstitutiveBase *> &rConstitutiveLaws,
                                                                                const NuTo::ElementBase *rElement)
{
    for(unsigned int i=0; i<rConstitutiveLaws.size(); ++i)
    {
        ConstitutiveStaticDataBase* staticData = rConstitutiveLaws[i]->AllocateStaticData3D(rElement);
        if(staticData != nullptr)
        {
            mStaticData.emplace(rConstitutiveLaws[i]->GetType(), staticData);
        }
    }
}


NuTo::ConstitutiveStaticDataMoistureTransport *NuTo::ConstitutiveStaticDataMultipleConstitutiveLaws::AsMoistureTransport()
{
    auto itConstLawStaticData = mStaticData.find(Constitutive::eConstitutiveType::MOISTURE_TRANSPORT);
    if (itConstLawStaticData == mStaticData.end())
    {
        throw MechanicsException(__PRETTY_FUNCTION__,"No constitutive static data for moisture transport found!");
    }
    return itConstLawStaticData->second->AsMoistureTransport();
}


const NuTo::ConstitutiveStaticDataMoistureTransport *NuTo::ConstitutiveStaticDataMultipleConstitutiveLaws::AsMoistureTransport() const
{
    auto itConstLawStaticData = mStaticData.find(Constitutive::eConstitutiveType::MOISTURE_TRANSPORT);
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
    case Constitutive::eConstitutiveType::ADDITIVE_OUTPUT:
        break;
    default:
        return false;
    }

    switch(rElementType)
    {
    case Element::eElementType::CONTINUUMELEMENT:
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
