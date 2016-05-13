#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveScalar.h"

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

NuTo::ConstitutiveIOBase& NuTo::ConstitutiveIOBase::operator=(const ConstitutiveIOBase& rOther)
{
    assert(GetNumColumns() == rOther.GetNumColumns() && "NumColumns must be equal");
    assert(GetNumRows() == rOther.GetNumRows() && "NumRows must be equal");

    for (int iCol = 0; iCol < GetNumColumns(); ++iCol)
        for (int iRow = 0; iRow < GetNumRows(); ++iRow)
            (*this)(iRow, iCol) = rOther(iRow,iCol);
    return *this;
}

double& NuTo::ConstitutiveIOBase::operator ()(int rRow, int rCol)
{
    throw MechanicsException(std::string("[")+__PRETTY_FUNCTION__+"] not supported.");
}


double NuTo::ConstitutiveIOBase::operator ()(int rRow, int rCol) const
{
    throw MechanicsException(std::string("[")+__PRETTY_FUNCTION__+"] not supported.");
}


double& NuTo::ConstitutiveIOBase::operator [](int rRow)
{
    throw MechanicsException(std::string("[")+__PRETTY_FUNCTION__+"] not supported.");
}

double NuTo::ConstitutiveIOBase::operator [](int rRow) const
{
    throw MechanicsException(std::string("[")+__PRETTY_FUNCTION__+"] not supported.");
}


void NuTo::ConstitutiveIOBase::AssertIsScalar(Constitutive::Output::eOutput rOutputEnum, std::string rMethodName) const
{
#ifdef DEBUG
    bool isNotScalar = dynamic_cast<const ConstitutiveScalar*>(this) == nullptr;
    if (isNotScalar)
        throw MechanicsException(std::string("[") + rMethodName + "] \n + Constitutive output " +
                Constitutive::OutputToString(rOutputEnum) + " is not a ConstitutiveScalar.");
#endif
}




#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::ConstitutiveIOBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveIOBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveIOBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveIOBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveIOBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveIOBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template <class Archive>
void NuTo::ConstitutiveIOBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
std::cout << "start serialize ConstitutiveIOBase" << std::endl;
#endif

#ifdef DEBUG_SERIALIZATION
std::cout << "finish serialize ConstitutiveIOBase" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
