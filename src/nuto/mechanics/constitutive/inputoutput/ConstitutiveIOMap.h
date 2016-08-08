#pragma once

#include <map>
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"

namespace NuTo
{

template<typename IOEnum>
class ConstitutiveIOMap : public std::map<IOEnum, std::unique_ptr<ConstitutiveIOBase>>
{
public:
    ConstitutiveIOMap() = default;
    ConstitutiveIOMap(const ConstitutiveIOMap& other);
    void Merge(const ConstitutiveIOMap& other);
};

using ConstitutiveInputMap = ConstitutiveIOMap<Constitutive::Input::eInput>;
using ConstitutiveOutputMap = ConstitutiveIOMap<Constitutive::Output::eOutput>;
} // namespace NuTo
