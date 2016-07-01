#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveMatrixXd.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

//template<int TDim>
//NuTo::ConstitutiveIOBase* NuTo::ConstitutiveIOBase::makeConstitutiveIO(
//        NuTo::Constitutive::Output::eOutput outputType)
//{
//
//}
//
//template NuTo::ConstitutiveIOBase* NuTo::ConstitutiveIOBase::makeConstitutiveIO<1>(
//        NuTo::Constitutive::Output::eOutput outputType);
//template NuTo::ConstitutiveIOBase* NuTo::ConstitutiveIOBase::makeConstitutiveIO<2>(
//        NuTo::Constitutive::Output::eOutput outputType);
//template NuTo::ConstitutiveIOBase* NuTo::ConstitutiveIOBase::makeConstitutiveIO<3>(
//        NuTo::Constitutive::Output::eOutput outputType);

template<int TDim>
std::unique_ptr<NuTo::ConstitutiveIOBase> NuTo::ConstitutiveIOBase::makeConstitutiveIO(
        NuTo::Constitutive::Input::eInput inputType)
{
    switch (inputType)
    {
    case Constitutive::Input::ENGINEERING_STRAIN:
        return std::make_unique<EngineeringStrain<TDim>>();
    case Constitutive::Input::NONLOCAL_EQ_STRAIN:
        return std::make_unique<ConstitutiveScalar>();
    case Constitutive::Input::RELATIVE_HUMIDITY:
        return std::make_unique<ConstitutiveScalar>();
    case Constitutive::Input::RELATIVE_HUMIDITY_DT1:
        return std::make_unique<ConstitutiveScalar>();
    case Constitutive::Input::RELATIVE_HUMIDITY_GRADIENT:
        return std::make_unique<ConstitutiveVector<TDim>>();
    case Constitutive::Input::TEMPERATURE:
        return std::make_unique<ConstitutiveScalar>();
    case Constitutive::Input::TEMPERATURE_GRADIENT:
        return std::make_unique<ConstitutiveVector<TDim>>();
    case Constitutive::Input::TEMPERATURE_CHANGE:
        return std::make_unique<ConstitutiveScalar>();
    case Constitutive::Input::WATER_VOLUME_FRACTION:
        return std::make_unique<ConstitutiveScalar>();
    case Constitutive::Input::WATER_VOLUME_FRACTION_DT1:
        return std::make_unique<ConstitutiveScalar>();
    case Constitutive::Input::WATER_VOLUME_FRACTION_GRADIENT:
        return std::make_unique<ConstitutiveVector<TDim>>();
    case Constitutive::Input::INTERFACE_SLIP:
        return std::make_unique<ConstitutiveMatrixXd>();
    case Constitutive::Input::CRACK_PHASE_FIELD:
        return std::make_unique<ConstitutiveScalar>();
    case Constitutive::Input::CALCULATE_STATIC_DATA:
        return std::make_unique<ConstitutiveCalculateStaticData>(NuTo::CalculateStaticData::EULER_BACKWARD);
    default:
        throw MechanicsException(__PRETTY_FUNCTION__,
                "Don't know how to create Constitutive input for this input type");
    }
}

template std::unique_ptr<NuTo::ConstitutiveIOBase> NuTo::ConstitutiveIOBase::makeConstitutiveIO<1>(
        NuTo::Constitutive::Input::eInput inputType);
template std::unique_ptr<NuTo::ConstitutiveIOBase> NuTo::ConstitutiveIOBase::makeConstitutiveIO<2>(
        NuTo::Constitutive::Input::eInput inputType);
template std::unique_ptr<NuTo::ConstitutiveIOBase> NuTo::ConstitutiveIOBase::makeConstitutiveIO<3>(
        NuTo::Constitutive::Input::eInput inputType);


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
