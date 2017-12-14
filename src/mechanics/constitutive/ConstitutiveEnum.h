#pragma once

#include <map>


namespace NuTo
{

namespace Constitutive
{
enum class eConstitutiveType
{
    ADDITIVE_INPUT_EXPLICIT, //!< container for multiple constitutive laws linked by a shared input
    ADDITIVE_INPUT_IMPLICIT, //!< container for multiple constitutive laws linked by a shared input
    ADDITIVE_OUTPUT, //!< container for multiple constitutive laws linked by addition of their outputs
    CREEP, //!< Kelvin chain based creep model
    DAMAGE_VISCO_PLASTICITY_ENGINEERING_STRESS, //!< viscoplastic damage model
    DAMAGE_VISCO_PLASTICITY_HARDENING_ENGINEERING_STRESS, //!< viscoplastic damage model with hardening
    FIBRE_MATRIX_BOND_STRESS_SLIP, //!< material model for the matrix-fibre interface
    GRADIENT_DAMAGE_ENGINEERING_STRESS, //!< gradient damage model
    GRADIENT_DAMAGE_FATIGUE_ENGINEERING_STRESS, //!< gradient damage model
    GRADIENT_DAMAGE_ENGINEERING_STRESS_FATIGUE, //!< gradient damage model for fatigued
    GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS, //!< gradient damage plasticity model
    HEAT_CONDUCTION, //!< Heat conduction
    LATTICE_CONCRETE, //!< material law for lattice model
    LINEAR_DAMPING_ENGINEERING_STRESS, //!< linear damping
    LINEAR_ELASTIC_ENGINEERING_STRESS, //!< linear-elastic behavior
    LINEAR_ELASTIC_INHOMOGENEOUS, //!<
    LINEAR_SPRING, //!< linear spring model
    LOCAL_DAMAGE_MODEL, //!< local damage model
    MISES_PLASTICITY_ENGINEERING_STRESS, //!< mises plasticity with isotropic and kinematic hardening
    MOISTURE_TRANSPORT, //!< moisture transport model
    MULTISCALE, //!< multiscale model, where the average stress is calculated from a full fine scale model
    NONLOCAL_DAMAGE_PLASTICITY_ENGINEERING_STRESS, //!< nonlocal damage model with plasticity in the effective stress
    //! space
    PHASE_FIELD, //!< phase field model
    SHRINKAGE_CAPILLARY_STRAIN_BASED, //!< strain based drying shrinkage - capillary term
    SHRINKAGE_CAPILLARY_STRESS_BASED, //!< stress based drying shrinkage - capillary term
    STRAIN_GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS, //!< strain gradient damage plasticity model (damage and
    //! plasticity are function of nonlocal total strain)
    THERMAL_STRAINS, //!< strain induced by temperature change
    LINEAR_ELASTIC_ANISOTROPIC, //!< linear elastic fully anisotropic material
    LINEAR_DIELECTRIC, //!< linear isotropic dielectric material (insulating but polarizable)
    LINEAR_PIEZOELECTRIC, //!< linear piezoelectric material (fully anisotropic)
    WAHWAHWAH
};

const std::map<eConstitutiveType, std::string> GetConstitutiveTypeMap();
std::string ConstitutiveTypeToString(eConstitutiveType rOutput);
eConstitutiveType ConstitutiveTypeToEnum(std::string rOutput);


enum class eConstitutiveStaticDataType
{
    MOISTURE_TRANSPORT,
    MULTI_PHYSICS
};


enum class eNonlocalDamageYieldSurface
{
    RANKINE_ROUNDED, //!< Rankine yield surface (rounded in tension
    ONLY_DP, //!< only Drucker-Prager yield surface
    COMBINED_ROUNDED, //!< Drucker-Prager and rounded Rankine
    COMBINED_SHARP, //!< Drucker-Prager and Rankine (no rounding)
    TWO_DP //!< two Drucker Prager cones
};

enum class eSolutionPhaseType
{
    HOMOGENIZED_LINEAR_ELASTIC, //!< linear (homogenized tangent stiffness is used to calculate stiffness and stress)
    NONLINEAR_NO_CRACK, //!< nonlinear with disp boundary conditions, but no crack can localize
    NONLINEAR_CRACKED //!< nonlinear with crack enrichment
};

enum class ePhaseFieldEnergyDecomposition
{
    ISOTROPIC, //!< isotropic degradation
    ANISOTROPIC_SPECTRAL_DECOMPOSITION, //!< spectral decompostion of the elastic strain tensor proposed by Miehe et al.
};

enum class eConstitutiveParameter
{
    ALPHA, //!<
    ARTIFICIAL_VISCOSITY, //!<
    BIAXIAL_COMPRESSIVE_STRENGTH, //!<
    BOUNDARY_DIFFUSION_COEFFICIENT_RH, //!<
    BOUNDARY_DIFFUSION_COEFFICIENT_WV, //!<
    COMPRESSIVE_STRENGTH, //!<
    DAMPING_COEFFICIENT, //!<
    DAMAGE_DISTRIBUTION, //!<
    DAMAGE_LAW, //!<
    DENSITY, //!<
    DENSITY_WATER, //!<
    DIELECTRIC_TENSOR, //!<
    DIFFUSION_COEFFICIENT_RH, //!<
    DIFFUSION_COEFFICIENT_WV, //!<
    DIFFUSION_EXPONENT_RH, //!<
    DIFFUSION_EXPONENT_WV, //!<
    ENABLE_MODIFIED_TANGENTIAL_STIFFNESS, //!<
    ENABLE_SORPTION_HYSTERESIS, //!<
    FATIGUE_EXTRAPOLATION, //!<
    FRACTURE_ENERGY, //!<
    GRADIENT_CORRECTION_ADSORPTION_DESORPTION, //!<
    GRADIENT_CORRECTION_DESORPTION_ADSORPTION, //!<
    HARDENING_EXPONENT, //!<
    HARDENING_VALUE, //!<
    HEAT_CAPACITY, //!< specific heat capacity \f$c_T\f$
    INITIAL_HARDENING_MODULUS, //!<
    INITIAL_YIELD_STRENGTH, //!<
    KELVIN_CHAIN_DAMPING, //!<
    KELVIN_CHAIN_RETARDATIONTIME, //!<
    KELVIN_CHAIN_STIFFNESS, //!<
    LENGTH_SCALE_PARAMETER, //!<
    MACROSCOPIC_BULK_MODULUS, //!<
    MASS_EXCHANGE_RATE, //!<
    MAX_BOND_STRESS, //!<
    NONLOCAL_RADIUS, //!<
    NORMAL_STIFFNESS, //!<
    PIEZOELECTRIC_TENSOR, //!<
    POISSONS_RATIO, //!<
    POLYNOMIAL_COEFFICIENTS_ADSORPTION, //!<
    POLYNOMIAL_COEFFICIENTS_DESORPTION, //!<
    PORE_VOLUME_FRACTION, //!<
    RESIDUAL_BOND_STRESS, //!<
    DENSITY_SATURATED_WATER_VAPOR, //!<
    SLIP_AT_MAX_BOND_STRESS, //!<
    SLIP_AT_RESIDUAL_BOND_STRESS, //!<
    SOLID_PHASE_BULK_MODULUS, //!<
    SPRING_STIFFNESS, //!<
    SPRING_DIRECTION, //!<
    STIFFNESS, //!<
    TENSILE_STRENGTH, //!<
    TEMPERATURE, //!<
    THERMAL_EXPANSION_COEFFICIENT, //!<
    THERMAL_CONDUCTIVITY, //!< \f$k \text{ in } \mathbf{q} = - k \nabla T \f$
    VISCOPLASTIC_YIELD_SURFACE_OFFSET, //!<
    VISCOSITY, //!<
    VISCOSITY_EXPONENT, //!<
    YOUNGS_MODULUS, //!<
    ENDURANCE_STRESS, //!<
    FATIGUE_PARAMETER, //!<

};

const std::map<eConstitutiveParameter, std::string> GetConstitutiveParameterMap();
std::string ConstitutiveParameterToString(eConstitutiveParameter rParameter);
eConstitutiveParameter ConstitutiveParameterToEnum(std::string rParameter);


enum class eInput
{
    NONE, //!< for constitutive law additive input explicit (AddConstitutiveLaw-function)
    ELECTRIC_FIELD, //!<
    ENGINEERING_STRAIN, //!<
    ENGINEERING_STRAIN_DT1, //!< first time derivative of the engineering strain
    TEMPERATURE, //!< Temperature \f$T\f$
    TEMPERATURE_GRADIENT, //!< Temperature gradient \f$\nabla T\f$
    TEMPERATURE_CHANGE, //!< First time derivative \f$\frac{\partial T}{\partial t}\f$
    NONLOCAL_EQ_PLASTIC_STRAIN, //!<
    NONLOCAL_EQ_STRAIN, //!<
    NONLOCAL_TOTAL_STRAIN_1D, //!<
    ENGINEERING_STRESS_1D, //!< usually the stress is an output, that 's why the additional input term is required
    DEFORMATION_GRADIENT_REAL_1D, //!<
    NONLOCAL_TOTAL_STRAIN_REAL_1D, //!<
    NONLOCAL_TOTAL_STRAIN_VIRT_1D, //!<
    RELATIVE_HUMIDITY, //!<
    RELATIVE_HUMIDITY_BOUNDARY, //!< control node relative humidity
    RELATIVE_HUMIDITY_DT1, //!< first time derivative
    RELATIVE_HUMIDITY_GRADIENT, //!<
    WATER_VOLUME_FRACTION, //!<
    WATER_VOLUME_FRACTION_DT1, //!< first time derivative
    WATER_VOLUME_FRACTION_GRADIENT, //!<
    INTERFACE_SLIP, //!<
    CRACK_PHASE_FIELD, //!<
    ELASTIC_ENERGY_DENSITY, //!<
    CALCULATE_STATIC_DATA,
    PLANE_STATE, //!< for telling the law whether there is plane stress or plane strain
    TIME,
    TIME_STEP,
    CALCULATE_INITIALIZE_VALUE_RATES, //!
    COORDINATES
};


std::string InputToString(const eInput& e);


enum class eOutput
{
    ELECTRIC_DISPLACEMENT,
    ELECTRIC_POTENTIAL,
    ELECTRIC_FIELD,
    D_ELECTRIC_DISPLACEMENT_D_ELECTRIC_FIELD,
    D_ELECTRIC_DISPLACEMENT_D_ENGINEERING_STRAIN,
    PIEZOELECTRIC_STRESS,
    ENGINEERING_STRESS,
    ENGINEERING_STRESS_VISUALIZE,
    D_ENGINEERING_STRESS_D_PHASE_FIELD,
    ENGINEERING_STRAIN,
    D_ENGINEERING_STRAIN_D_RELATIVE_HUMIDITY,
    D_ENGINEERING_STRAIN_D_WATER_VOLUME_FRACTION,
    ENGINEERING_STRAIN_VISUALIZE,
    SHRINKAGE_STRAIN_VISUALIZE,
    THERMAL_STRAIN, //!< \f$\varepsilon_{th} = \alpha \Delta T I\f$
    D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN,
    D_ENGINEERING_STRESS_D_ELECTRIC_FIELD,
    D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_DT1,
    D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY,
    D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION,
    ENGINEERING_PLASTIC_STRAIN_VISUALIZE,

    D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN,
    D_LOCAL_EQ_STRAIN_D_STRAIN,
    ENGINEERING_VISCOPLASTIC_STRAIN_3D, //!<
    ENGINEERING_TOTAL_INELASTIC_STRAIN_3D, //!<
    ENGINEERING_STRESS_ELASTIC_1D, //!<
    ENGINEERING_STRESS_ELASTIC_2D, //!<
    ENGINEERING_STRESS_ELASTIC_3D, //!<
    D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_ELASTIC_1D,
    D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_ELASTIC_2D,
    D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_ELASTIC_3D,
    D_ENGINEERING_STRESS_D_THERMAL_STRAIN,
    D_ENGINEERING_STRESS_D_NONLOCAL_EQ_PLASTIC_STRAIN_1D,
    D_ENGINEERING_STRESS_D_NONLOCAL_TOTAL_STRAIN_1D,
    D_ENGINEERING_STRESS_D_TEMPERATURE,
    D_STRAIN_D_TEMPERATURE,
    D_ENGINEERING_STRESS_D_TEMPERATURE_1D,
    D_ENGINEERING_STRESS_D_TEMPERATURE_2D,
    D_ENGINEERING_STRESS_D_TEMPERATURE_3D,
    HEAT_FLUX,
    HEAT_CHANGE, //!< First time derivative of heat, i.e.
    //!< \f$\frac{\partial Q}{\partial t} = \rho c_T \frac{\partial T}{\partial t}\f$
    D_HEAT_FLUX_D_TEMPERATURE_GRADIENT, //!< conductivity matrix
    D_HEAT_D_TEMPERATURE, //!< heat capacity
    DAMAGE,
    EXTRAPOLATION_ERROR,
    UPDATE_STATIC_DATA,
    UPDATE_TMP_STATIC_DATA,
    LOCAL_EQ_PLASTIC_STRAIN,
    LOCAL_EQ_TOTAL_INELASTIC_STRAIN,
    LOCAL_EQ_STRAIN,
    D_LOCAL_EQ_PLASTIC_STRAIN_D_STRAIN_1D,
    ENGINEERING_STRESS_REAL_1D,
    ENGINEERING_STRAIN_VIRT_1D,
    D_ENGINEERING_STRESS_REAL_D_ENGINEERING_STRAIN_REAL_1D,
    D_ENGINEERING_STRESS_REAL_D_NONLOCAL_TOTAL_STRAIN_REAL_1D,
    D_ENGINEERING_STRAIN_VIRT_D_STRESS_REAL_1D,
    D_ENGINEERING_STRAIN_VIRT_D_NONLOCAL_TOTAL_STRAIN_VIRT_1D,
    D_INTERNAL_GRADIENT_RH_D_RH_BB_H0,
    D_INTERNAL_GRADIENT_RH_D_RH_NN_H0,
    D_INTERNAL_GRADIENT_RH_D_WV_BN_H0,
    D_INTERNAL_GRADIENT_RH_D_WV_NN_H0,
    D_INTERNAL_GRADIENT_WV_D_RH_NN_H0,
    D_INTERNAL_GRADIENT_WV_D_WV_BB_H0,
    D_INTERNAL_GRADIENT_WV_D_WV_BN_H0,
    D_INTERNAL_GRADIENT_WV_D_WV_NN_H0,
    D_INTERNAL_GRADIENT_RH_D_RH_NN_H1,
    D_INTERNAL_GRADIENT_RH_D_WV_NN_H1,
    D_INTERNAL_GRADIENT_WV_D_WV_NN_H1,
    INTERNAL_GRADIENT_RELATIVE_HUMIDITY_B,
    INTERNAL_GRADIENT_RELATIVE_HUMIDITY_N,
    INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_B,
    INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_N,
    D_INTERNAL_GRADIENT_RH_D_RH_BOUNDARY_NN_H0,
    D_INTERNAL_GRADIENT_WV_D_WV_BOUNDARY_NN_H0,
    INTERNAL_GRADIENT_RELATIVE_HUMIDITY_BOUNDARY_N,
    INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_BOUNDARY_N,
    FATIGUE_SAVE_STATIC_DATA,
    FATIGUE_RESTORE_STATIC_DATA,
    FATIGUE_EXTRAPOLATE_STATIC_DATA,
    INTERFACE_CONSTITUTIVE_MATRIX,
    BOND_STRESS,
    SLIP,
    ELASTIC_ENERGY_DAMAGED_PART,
    D_ELASTIC_ENERGY_DAMAGED_PART_D_ENGINEERING_STRAIN,
    NONLOCAL_RADIUS
};


std::string OutputToString(const eOutput e);
} // Constitutive


} // NuTo
