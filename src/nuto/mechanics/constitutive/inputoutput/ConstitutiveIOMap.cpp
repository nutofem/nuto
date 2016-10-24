#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"

template<typename IOEnum>
NuTo::ConstitutiveIOMap<IOEnum>::ConstitutiveIOMap(const ConstitutiveIOMap<IOEnum>& other)
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

template<typename IOEnum>
NuTo::ConstitutiveIOMap<IOEnum>& NuTo::ConstitutiveIOMap<IOEnum>::Merge(const ConstitutiveIOMap<IOEnum>& other)
{
    for(auto& it : other)
    {
        // if the key exists, and there is an object associated with it, throw
        // BEWARE: this uses short-circuit evaluation,
        // i.e. you can't change the order of the two conditionals
        if (this->count(it.first) != 0 and this->at(it.first) != nullptr)
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "Duplicate key in constitutive maps.");
        }
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
    return *this;
}

template<typename IOEnum>
bool NuTo::ConstitutiveIOMap<IOEnum>::Contains(IOEnum rEnum) const
{
    return this->find(rEnum) != this->end();
}


template class NuTo::ConstitutiveIOMap<NuTo::Constitutive::eInput>;
template class NuTo::ConstitutiveIOMap<NuTo::Constitutive::eOutput>;
