#include "nuto/mechanics/constitutive/staticData/Composite.h"

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
