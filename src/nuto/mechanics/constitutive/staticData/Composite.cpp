#include "nuto/mechanics/constitutive/staticData/Composite.h"
#include "nuto/base/serializeStream/SerializeStreamIn.h"
#include "nuto/base/serializeStream/SerializeStreamOut.h"

using namespace NuTo::Constitutive::StaticData;

Composite* Composite::Create()
{
    return new Composite();
}


Composite* Composite::Clone() const
{
    auto newComposite = new Composite;
    for (const auto& it : this->mComponents)
    {
        newComposite->AddComponent(it.Clone());
    }
    return newComposite;
}


void Composite::AddComponent(Component* newComponent)
{
    mComponents.push_back(newComponent);
}


Component& Composite::GetComponent(unsigned int index)
{
    return mComponents.at(index);
}


unsigned int Composite::GetNumComponents() const
{
    return mComponents.size();
}


void Composite::ShiftToPast()
{
    for (auto& it : mComponents)
    {
        it.ShiftToPast();
    }
}


void Composite::ShiftToFuture()
{
    for (auto& it : mComponents)
    {
        it.ShiftToFuture();
    }
}


void Composite::AllocateAdditionalData(int numAdditionalData)
{
    for (auto& it : mComponents)
    {
        it.AllocateAdditionalData(numAdditionalData);
    }
}

void Composite::NuToSerializeSave(SerializeStreamOut &rStream)
{
    Component::NuToSerializeSave(rStream);
    SerializeComposite(rStream);
}

void Composite::NuToSerializeLoad(SerializeStreamIn &rStream)
{
    Component::NuToSerializeLoad(rStream);
    SerializeComposite(rStream);
}

template <typename TStream>
void Composite::SerializeComposite(TStream &rStream)
{
    for (Component& component : mComponents)
    {
        rStream.Serialize(component);
    }
}

template void Composite::SerializeComposite<NuTo::SerializeStreamIn>(SerializeStreamIn& rStream);
template void Composite::SerializeComposite<NuTo::SerializeStreamOut>(SerializeStreamOut& rStream);
