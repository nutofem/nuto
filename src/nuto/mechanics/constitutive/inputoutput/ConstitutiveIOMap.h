#pragma once

#include <map>
#include <memory>

namespace NuTo
{

// Forward declarations
class ConstitutiveIOBase;
namespace Constitutive
{
    enum class eInput;
    enum class eOutput;
}



template<typename IOEnum>
class ConstitutiveIOMap : public std::map<IOEnum, std::unique_ptr<ConstitutiveIOBase>>
{
public:
    ConstitutiveIOMap() = default;
    ConstitutiveIOMap(const ConstitutiveIOMap& other);
    NuTo::ConstitutiveIOMap<IOEnum>& Merge(const ConstitutiveIOMap& other);
};

using ConstitutiveInputMap = ConstitutiveIOMap<Constitutive::eInput>;
using ConstitutiveOutputMap = ConstitutiveIOMap<Constitutive::eOutput>;
} // namespace NuTo
