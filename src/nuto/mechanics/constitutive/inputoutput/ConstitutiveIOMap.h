#pragma once

#include <map>
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"

namespace NuTo
{

template<typename IOEnum>
class ConstitutiveIOMap : public std::map<IOEnum, std::unique_ptr<ConstitutiveIOBase>>
{
public:
    //! Join this map with an `other` map. Copies entries into the map.
    void Join(const ConstitutiveIOMap& other);
};

using ConstitutiveInputMap = ConstitutiveIOMap<Constitutive::Input::eInput>;
using ConstitutiveOutputMap = std::map<Constitutive::Output::eOutput, ConstitutiveIOBase*>;
} // namespace NuTo
