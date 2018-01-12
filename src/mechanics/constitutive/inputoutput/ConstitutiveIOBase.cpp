#include "mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrainAxSy.h"
#include "mechanics/constitutive/inputoutput/EngineeringStressAxSy.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveMatrixXd.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "mechanics/constitutive/inputoutput/ConstitutivePlaneState.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"

template <int TDim>
std::unique_ptr<NuTo::ConstitutiveIOBase>
NuTo::ConstitutiveIOBase::makeConstitutiveIO(NuTo::Constitutive::eOutput outputType)
{
    using namespace Constitutive;
    constexpr int VoigtDim = ConstitutiveIOBase::GetVoigtDim(TDim);
    switch (outputType)
    {
    // scalars
    case eOutput::DAMAGE:
    case eOutput::D_HEAT_D_TEMPERATURE:
    case eOutput::D_INTERNAL_GRADIENT_RH_D_RH_BOUNDARY_NN_H0:
    case eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BOUNDARY_NN_H0:
    case eOutput::D_INTERNAL_GRADIENT_RH_D_RH_BB_H0:
    case eOutput::D_INTERNAL_GRADIENT_RH_D_RH_NN_H0:
    case eOutput::D_INTERNAL_GRADIENT_RH_D_WV_NN_H0:
    case eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BB_H0:
    case eOutput::D_INTERNAL_GRADIENT_WV_D_WV_NN_H0:
    case eOutput::D_INTERNAL_GRADIENT_WV_D_RH_NN_H0:
    case eOutput::D_INTERNAL_GRADIENT_RH_D_RH_NN_H1:
    case eOutput::D_INTERNAL_GRADIENT_RH_D_WV_NN_H1:
    case eOutput::D_INTERNAL_GRADIENT_WV_D_WV_NN_H1:
    case eOutput::EXTRAPOLATION_ERROR:
    case eOutput::ELASTIC_ENERGY_DAMAGED_PART:
    case eOutput::HEAT_CHANGE:
    case eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_N:
    case eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_BOUNDARY_N:
    case eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_N:
    case eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_BOUNDARY_N:
    case eOutput::LOCAL_EQ_STRAIN:
    case eOutput::NONLOCAL_RADIUS:
        return std::make_unique<ConstitutiveScalar>();
    // vectors dim
    case eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_B:
    case eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_B:
    case eOutput::D_INTERNAL_GRADIENT_RH_D_WV_BN_H0:
    case eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BN_H0:
    case eOutput::HEAT_FLUX:
    case eOutput::ELECTRIC_DISPLACEMENT:
    case eOutput::ELECTRIC_FIELD:
        return std::make_unique<ConstitutiveVector<TDim>>();
    // vectors voigtdim
    case eOutput::D_ELASTIC_ENERGY_DAMAGED_PART_D_ENGINEERING_STRAIN:
    case eOutput::D_ENGINEERING_STRAIN_D_RELATIVE_HUMIDITY:
    case eOutput::D_ENGINEERING_STRAIN_D_WATER_VOLUME_FRACTION:
    case eOutput::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN:
    case eOutput::D_ENGINEERING_STRESS_D_PHASE_FIELD:
    case eOutput::D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY:
    case eOutput::D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION:
    case eOutput::D_ENGINEERING_STRESS_D_TEMPERATURE:
    case eOutput::D_LOCAL_EQ_STRAIN_D_STRAIN:
    case eOutput::D_STRAIN_D_TEMPERATURE:
        return std::make_unique<ConstitutiveVector<VoigtDim>>();
    case eOutput::D_ELASTIC_ENERGY_DAMAGED_PART_D_ENGINEERING_STRAIN_AxSy:
    case eOutput::D_ENGINEERING_STRESS_AxSy_D_NONLOCAL_EQ_STRAIN:
    case eOutput::D_ENGINEERING_STRESS_AxSy_D_PHASE_FIELD:
    case eOutput::D_LOCAL_EQ_STRAIN_D_STRAIN_AxSy:
    	return std::make_unique<ConstitutiveVector<4>>();
    // visualize
    case eOutput::ENGINEERING_PLASTIC_STRAIN_VISUALIZE:
    case eOutput::ENGINEERING_STRAIN_VISUALIZE:
    case eOutput::SHRINKAGE_STRAIN_VISUALIZE:
    case eOutput::THERMAL_STRAIN:
        return std::make_unique<EngineeringStrain<3>>();
    case eOutput::ENGINEERING_STRESS_VISUALIZE:
        return std::make_unique<EngineeringStress<3>>();
    // other
    case eOutput::ENGINEERING_STRESS:
        return std::make_unique<EngineeringStress<TDim>>();
    case eOutput::ENGINEERING_STRESS_AxSy:
    	return std::make_unique<EngineeringStressAxSy>();
    case eOutput::ENGINEERING_STRAIN:
        return std::make_unique<EngineeringStrain<TDim>>();
    case eOutput::ENGINEERING_STRAIN_AxSy:
    	return std::make_unique<EngineeringStrainAxSy>();
    case eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
    case eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_DT1:
        return std::make_unique<ConstitutiveMatrix<VoigtDim, VoigtDim>>();
    case eOutput::D_ENGINEERING_STRESS_AxSy_D_ENGINEERING_STRAIN_AxSy:
    case eOutput::D_ENGINEERING_STRESS_AxSy_D_ENGINEERING_STRAIN_AxSy_DT1:
        return std::make_unique<ConstitutiveMatrix<4, 4>>();
    case eOutput::D_ELECTRIC_DISPLACEMENT_D_ELECTRIC_FIELD:
        return std::make_unique<ConstitutiveMatrix<TDim, TDim>>();
    case eOutput::D_ELECTRIC_DISPLACEMENT_D_ENGINEERING_STRAIN:
        return std::make_unique<ConstitutiveMatrix<TDim, VoigtDim>>();
    case eOutput::D_ENGINEERING_STRESS_D_ELECTRIC_FIELD:
        return std::make_unique<ConstitutiveMatrix<VoigtDim, TDim>>();
    case eOutput::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT:
        return std::make_unique<ConstitutiveMatrix<TDim, TDim>>();
    case eOutput::BOND_STRESS:
        return std::make_unique<ConstitutiveMatrixXd>();
    case eOutput::INTERFACE_CONSTITUTIVE_MATRIX:
        return std::make_unique<ConstitutiveMatrixXd>();
    case eOutput::UPDATE_STATIC_DATA:
    case eOutput::AXISYMMETRIC_STATE:
        return nullptr;
    default:
        throw Exception(__PRETTY_FUNCTION__, "Don't know how to create constitutive output for " +
                                                              Constitutive::OutputToString(outputType));
    }
}

template std::unique_ptr<NuTo::ConstitutiveIOBase>
NuTo::ConstitutiveIOBase::makeConstitutiveIO<1>(NuTo::Constitutive::eOutput outputType);
template std::unique_ptr<NuTo::ConstitutiveIOBase>
NuTo::ConstitutiveIOBase::makeConstitutiveIO<2>(NuTo::Constitutive::eOutput outputType);
template std::unique_ptr<NuTo::ConstitutiveIOBase>
NuTo::ConstitutiveIOBase::makeConstitutiveIO<3>(NuTo::Constitutive::eOutput outputType);

template <int TDim>
std::unique_ptr<NuTo::ConstitutiveIOBase>
NuTo::ConstitutiveIOBase::makeConstitutiveIO(NuTo::Constitutive::eInput inputType)
{
    using namespace Constitutive;
    switch (inputType)
    {
    // scalars
    case eInput::CRACK_PHASE_FIELD:
    case eInput::NONLOCAL_EQ_STRAIN:
    case eInput::RELATIVE_HUMIDITY:
    case eInput::RELATIVE_HUMIDITY_BOUNDARY:
    case eInput::RELATIVE_HUMIDITY_DT1:
    case eInput::WATER_VOLUME_FRACTION:
    case eInput::WATER_VOLUME_FRACTION_DT1:
    case eInput::TEMPERATURE:
    case eInput::TEMPERATURE_CHANGE:
        return std::make_unique<ConstitutiveScalar>();
    // vectors
    case eInput::RELATIVE_HUMIDITY_GRADIENT:
    case eInput::WATER_VOLUME_FRACTION_GRADIENT:
    case eInput::TEMPERATURE_GRADIENT:
    case eInput::ELECTRIC_FIELD:
    case eInput::COORDINATES:
        return std::make_unique<ConstitutiveVector<TDim>>();
    // other
    case eInput::ENGINEERING_STRAIN:
    case eInput::ENGINEERING_STRAIN_DT1:
        return std::make_unique<EngineeringStrain<TDim>>();
    case eInput::ENGINEERING_STRAIN_AxSy:
    case eInput::ENGINEERING_STRAIN_AxSy_DT1:
    	return std::make_unique<EngineeringStrainAxSy>();
    case eInput::INTERFACE_SLIP:
        return std::make_unique<ConstitutiveMatrixXd>();
    case eInput::CALCULATE_STATIC_DATA:
        return std::make_unique<ConstitutiveCalculateStaticData>(NuTo::eCalculateStaticData::EULER_BACKWARD);
    case eInput::PLANE_STATE:
        return std::make_unique<ConstitutivePlaneState>(NuTo::ePlaneState::PLANE_STRESS);
    default:
        throw Exception(__PRETTY_FUNCTION__,
                                 "Don't know how to create Constitutive input for this input type");
    }
}

template std::unique_ptr<NuTo::ConstitutiveIOBase>
NuTo::ConstitutiveIOBase::makeConstitutiveIO<1>(NuTo::Constitutive::eInput inputType);
template std::unique_ptr<NuTo::ConstitutiveIOBase>
NuTo::ConstitutiveIOBase::makeConstitutiveIO<2>(NuTo::Constitutive::eInput inputType);
template std::unique_ptr<NuTo::ConstitutiveIOBase>
NuTo::ConstitutiveIOBase::makeConstitutiveIO<3>(NuTo::Constitutive::eInput inputType);


NuTo::ConstitutiveIOBase& NuTo::ConstitutiveIOBase::operator=(const ConstitutiveIOBase& rOther)
{
    assert(GetNumColumns() == rOther.GetNumColumns() && "NumColumns must be equal");
    assert(GetNumRows() == rOther.GetNumRows() && "NumRows must be equal");

    for (int iCol = 0; iCol < GetNumColumns(); ++iCol)
        for (int iRow = 0; iRow < GetNumRows(); ++iRow)
            (*this)(iRow, iCol) = rOther(iRow, iCol);

    this->mIsCalculated = rOther.GetIsCalculated();
    return *this;
}

Eigen::MatrixXd NuTo::ConstitutiveIOBase::CopyToEigenMatrix() const
{
    int rows = GetNumRows();
    int cols = GetNumColumns();
    Eigen::MatrixXd m(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            m(i, j) = this->operator()(i, j);
    return m;
}


double& NuTo::ConstitutiveIOBase::operator()(int, int)
{
    throw Exception(std::string("[") + __PRETTY_FUNCTION__ + "] not supported.");
}


double NuTo::ConstitutiveIOBase::operator()(int, int) const
{
    throw Exception(std::string("[") + __PRETTY_FUNCTION__ + "] not supported.");
}


double& NuTo::ConstitutiveIOBase::operator[](int)
{
    throw Exception(std::string("[") + __PRETTY_FUNCTION__ + "] not supported.");
}

double NuTo::ConstitutiveIOBase::operator[](int) const
{
    throw Exception(std::string("[") + __PRETTY_FUNCTION__ + "] not supported.");
}


void NuTo::ConstitutiveIOBase::AssertIsScalar(Constitutive::eOutput rOutputEnum, std::string rMethodName) const
{
    bool isNotScalar = dynamic_cast<const ConstitutiveScalar*>(this) == nullptr;
    if (isNotScalar)
        throw Exception(rMethodName, "Constitutive output " + Constitutive::OutputToString(rOutputEnum) +
                                                      " is not a ConstitutiveScalar.");
}


namespace NuTo
{
template <int TDim>
EngineeringStrain<TDim>& ConstitutiveIOBase::AsEngineeringStrain()
{
    throw Exception(__PRETTY_FUNCTION__, "invalid dimension");
}

template <>
EngineeringStrain<1>& ConstitutiveIOBase::AsEngineeringStrain<1>()
{
    return AsEngineeringStrain1D();
}

template <>
EngineeringStrain<2>& ConstitutiveIOBase::AsEngineeringStrain<2>()
{
    return AsEngineeringStrain2D();
}

template <>
EngineeringStrain<3>& ConstitutiveIOBase::AsEngineeringStrain<3>()
{
    return AsEngineeringStrain3D();
}

//EngineeringStrainAxSy& ConstitutiveIOBase::AsEngineeringStrainAxSy()
//{
//    return AsEngineeringStrainAxSy();
//}

template <int TDim>
const EngineeringStrain<TDim>& ConstitutiveIOBase::AsEngineeringStrain() const
{
    throw Exception(__PRETTY_FUNCTION__, "invalid dimension");
}

template <>
const EngineeringStrain<1>& ConstitutiveIOBase::AsEngineeringStrain<1>() const
{
    return AsEngineeringStrain1D();
}

template <>
const EngineeringStrain<2>& ConstitutiveIOBase::AsEngineeringStrain<2>() const
{
    return AsEngineeringStrain2D();
}

template <>
const EngineeringStrain<3>& ConstitutiveIOBase::AsEngineeringStrain<3>() const
{
    return AsEngineeringStrain3D();
}

//const EngineeringStrainAxSy& ConstitutiveIOBase::AsEngineeringStrainAxSy() const
//{
//    return AsEngineeringStrainAxSy();
//}
}
