#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"

template<typename IOEnum>
void NuTo::ConstitutiveIOMap<IOEnum>::Join(const ConstitutiveIOMap<IOEnum>& other)
{
    for(auto& it : other)
    {
        // if the map has no object associated with that key, simply add the key
        if (it.second == nullptr)
        {
            (*this)[it.first];
        }
        // copy construction from object in `other` map 
        else
        {
            this->insert(std::make_pair(it.first, it.second->clone()));
        }
    }
}

template class NuTo::ConstitutiveIOMap<NuTo::Constitutive::Input::eInput>;
