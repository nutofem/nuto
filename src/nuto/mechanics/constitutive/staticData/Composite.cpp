#include "nuto/mechanics/constitutive/staticData/Composite.h"
#include "nuto/base/serializeStream/SerializeStreamIn_Def.h"
#include "nuto/base/serializeStream/SerializeStreamOut_Def.h"

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

void Composite::WriteComponent(SerializeStreamOut& rStream)
{
    SerializeComponent(rStream);
}

void Composite::ReadComponent(SerializeStreamIn& rStream)
{
    SerializeComponent(rStream);
}

template <typename TStream>
void Composite::SerializeComponent(TStream& rStream)
{

}

template void Composite::SerializeComponent<NuTo::SerializeStreamIn>(SerializeStreamIn& rStream);
template void Composite::SerializeComponent<NuTo::SerializeStreamOut>(SerializeStreamOut& rStream);
