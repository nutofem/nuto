#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"
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

template<int TDim>
std::unique_ptr<NuTo::ConstitutiveIOBase> NuTo::ConstitutiveIOBase::makeConstitutiveIO(
        NuTo::Constitutive::Output::eOutput outputType)
{
    using namespace Constitutive::Output;
    constexpr int VoigtDim = ConstitutiveIOBase::GetVoigtDim(TDim);
    switch (outputType)
    {
        // scalars
        case DAMAGE:
        case D_HEAT_D_TEMPERATURE:
        case D_INTERNAL_GRADIENT_RH_D_RH_BOUNDARY_NN_H0:
        case D_INTERNAL_GRADIENT_WV_D_WV_BOUNDARY_NN_H0:
        case D_INTERNAL_GRADIENT_RH_D_RH_BB_H0:
        case D_INTERNAL_GRADIENT_RH_D_RH_NN_H0:
        case D_INTERNAL_GRADIENT_RH_D_WV_NN_H0:
        case D_INTERNAL_GRADIENT_WV_D_WV_BB_H0:
        case D_INTERNAL_GRADIENT_WV_D_WV_NN_H0:
        case D_INTERNAL_GRADIENT_WV_D_RH_NN_H0:
        case D_INTERNAL_GRADIENT_RH_D_RH_NN_H1:
        case D_INTERNAL_GRADIENT_RH_D_WV_NN_H1:
        case D_INTERNAL_GRADIENT_WV_D_WV_NN_H1:
        case EXTRAPOLATION_ERROR:
        case ELASTIC_ENERGY_DAMAGED_PART:
        case HEAT_CHANGE:
        case INTERNAL_GRADIENT_RELATIVE_HUMIDITY_N:
        case INTERNAL_GRADIENT_RELATIVE_HUMIDITY_BOUNDARY_N:
        case INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_N:
        case INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_BOUNDARY_N:
        case NONLOCAL_PARAMETER_XI:
        case LOCAL_EQ_STRAIN:
            return std::make_unique<ConstitutiveScalar>();
        // vectors dim
        case INTERNAL_GRADIENT_RELATIVE_HUMIDITY_B:
        case INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_B:
        case D_INTERNAL_GRADIENT_RH_D_WV_BN_H0:
        case D_INTERNAL_GRADIENT_WV_D_WV_BN_H0:
        case HEAT_FLUX:
            return std::make_unique<ConstitutiveVector<TDim>>();
        // vectors voigtdim
        case D_ELASTIC_ENERGY_DAMAGED_PART_D_ENGINEERING_STRAIN:
        case D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN:
        case D_ENGINEERING_STRESS_D_PHASE_FIELD:
        case D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY:
        case D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION:
        case D_LOCAL_EQ_STRAIN_XI_D_STRAIN:
        case D_ENGINEERING_STRESS_D_TEMPERATURE:
        case D_LOCAL_EQ_STRAIN_D_STRAIN:
        case D_STRAIN_D_TEMPERATURE:
            return std::make_unique<ConstitutiveVector<VoigtDim>>();
        // visualize
        case ENGINEERING_PLASTIC_STRAIN_VISUALIZE:
        case ENGINEERING_STRAIN_VISUALIZE:
        case SHRINKAGE_STRAIN_VISUALIZE:
        case THERMAL_STRAIN:
            return std::make_unique<EngineeringStrain<3>>();
        case ENGINEERING_STRESS_VISUALIZE:
            return std::make_unique<EngineeringStress<3>>();
        // other
        case ENGINEERING_STRESS:
            return std::make_unique<EngineeringStress<TDim>>();
        case ENGINEERING_STRAIN:
            return std::make_unique<EngineeringStrain<TDim>>();
        case D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
            return std::make_unique<ConstitutiveMatrix<VoigtDim, VoigtDim>>();
        case D_HEAT_FLUX_D_TEMPERATURE_GRADIENT:
            return std::make_unique<ConstitutiveMatrix<TDim, TDim>>();
        case BOND_STRESS:
            return std::make_unique<ConstitutiveMatrixXd>();
        case INTERFACE_CONSTITUTIVE_MATRIX:
            return std::make_unique<ConstitutiveMatrixXd>();
        case UPDATE_STATIC_DATA:
            return 0;
        default:
            throw MechanicsException(__PRETTY_FUNCTION__,
                    "Don't know how to create constitutive output for "
                    + Constitutive::OutputToString(outputType));
    }
}

template std::unique_ptr<NuTo::ConstitutiveIOBase> NuTo::ConstitutiveIOBase::makeConstitutiveIO<1>(
        NuTo::Constitutive::Output::eOutput outputType);
template std::unique_ptr<NuTo::ConstitutiveIOBase> NuTo::ConstitutiveIOBase::makeConstitutiveIO<2>(
        NuTo::Constitutive::Output::eOutput outputType);
template std::unique_ptr<NuTo::ConstitutiveIOBase> NuTo::ConstitutiveIOBase::makeConstitutiveIO<3>(
        NuTo::Constitutive::Output::eOutput outputType);

template<int TDim>
std::unique_ptr<NuTo::ConstitutiveIOBase> NuTo::ConstitutiveIOBase::makeConstitutiveIO(
        NuTo::Constitutive::Input::eInput inputType)
{
    using namespace Constitutive::Input;
    switch (inputType)
    {
        // scalars
        case CRACK_PHASE_FIELD:
        case NONLOCAL_EQ_STRAIN:
        case RELATIVE_HUMIDITY:
        case RELATIVE_HUMIDITY_DT1:
        case WATER_VOLUME_FRACTION:
        case WATER_VOLUME_FRACTION_DT1:
        case TEMPERATURE:
        case TEMPERATURE_CHANGE:
            return std::make_unique<ConstitutiveScalar>();
        // vectors
        case RELATIVE_HUMIDITY_GRADIENT:
        case WATER_VOLUME_FRACTION_GRADIENT:
        case TEMPERATURE_GRADIENT:
            return std::make_unique<ConstitutiveVector<TDim>>();
        // other
        case ENGINEERING_STRAIN:
            return std::make_unique<EngineeringStrain<TDim>>();
        case INTERFACE_SLIP:
            return std::make_unique<ConstitutiveMatrixXd>();
        case CALCULATE_STATIC_DATA:
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

    this->mIsCalculated = rOther.GetIsCalculated();
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
