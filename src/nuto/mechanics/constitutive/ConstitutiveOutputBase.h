// $Id$

#ifndef ConstitutiveOutputBase_H_
#define ConstitutiveOutputBase_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include <string>

#include "nuto/base/ErrorEnum.h"
#include "nuto/math/FullMatrix_Def.h"
#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"

namespace NuTo
{
// forward declarations
class ConstitutiveHeatFluxTemperatureGradient;
class ConstitutiveEngineeringStressStrain;
class ConstitutiveLatticeStressStrain;
class ConstitutiveTangentLocal1x1;
class ConstitutiveTangentLocal2x2;
class ConstitutiveTangentLocal3x3;
class ConstitutiveTangentLocal6x6;
class Damage;
class DeformationGradient1D;
class DeformationGradient2D;
class DeformationGradient3D;
class ElementBase;
class EngineeringStrain1D;
class EngineeringStrain2D;
class EngineeringStrain3D;
class EngineeringStress1D;
class EngineeringStress2D;
class EngineeringStress3D;
class HeatFlux3D;
class LocalEqPlasticStrain;
class LocalEqStrain;
class LocalEqTotalInelasticStrain;
class Logger;
class SecondPiolaKirchhoffStress3D;
class StructureBase;
template <int TNumRows,int TNumColumns>
class ConstitutiveTangentLocal;

//! @brief ... base class for the constitutive relationship, e.g. material laws
//! @author Jörg F. Unger, BAM
//! @date July 2012
class ConstitutiveOutputBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief ... constructor
    ConstitutiveOutputBase();

    //! @brief ... constructor
    virtual ~ConstitutiveOutputBase()
    {}

    virtual EngineeringStrain1D& GetEngineeringStrain1D();
    virtual EngineeringStrain2D& GetEngineeringStrain2D();
    virtual EngineeringStrain3D& GetEngineeringStrain3D();
    virtual EngineeringStress1D& GetEngineeringStress1D();
    virtual EngineeringStress2D& GetEngineeringStress2D();
    virtual EngineeringStress3D& GetEngineeringStress3D();
    virtual ConstitutiveTangentLocal3x3& GetConstitutiveTangentLocal3x3();
    virtual ConstitutiveTangentLocal6x6& GetConstitutiveTangentLocal6x6();
    virtual LocalEqPlasticStrain& GetLocalEqPlasticStrain();
    virtual LocalEqStrain& GetLocalEqStrain();
    virtual LocalEqTotalInelasticStrain& GetLocalEqTotalInelasticStrain();
    virtual HeatFlux3D& GetHeatFlux3D();
    virtual Damage& GetDamage();

    //! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
    virtual NuTo::ConstitutiveTangentLocal<1,1>& AsConstitutiveTangentLocal_1x1();

    //! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
    virtual NuTo::ConstitutiveTangentLocal<1,2>& AsConstitutiveTangentLocal_1x2();

    //! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
    virtual NuTo::ConstitutiveTangentLocal<2,1>& AsConstitutiveTangentLocal_2x1();

    //! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
    virtual NuTo::ConstitutiveTangentLocal<2,2>& AsConstitutiveTangentLocal_2x2();

    //! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
    virtual NuTo::ConstitutiveTangentLocal<3,1>& AsConstitutiveTangentLocal_3x1();

    //! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
    virtual NuTo::ConstitutiveTangentLocal<3,3>& AsConstitutiveTangentLocal_3x3();

    //! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
    virtual NuTo::ConstitutiveTangentLocal<6,1>& AsConstitutiveTangentLocal_6x1();

    //! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
    virtual NuTo::ConstitutiveTangentLocal<6,6>& AsConstitutiveTangentLocal_6x6();

    //! @brief return part of the nonlocal matrix
    virtual NuTo::ConstitutiveTangentLocal<1,1>& GetSubMatrix_1x1(int rSubMatrix);

    //! @brief return part of the nonlocal matrix
    virtual NuTo::ConstitutiveTangentLocal<2,2>& GetSubMatrix_2x2(int rSubMatrix);

    //! @brief return part of the nonlocal matrix
    virtual NuTo::ConstitutiveTangentLocal<3,3>& GetSubMatrix_3x3(int rSubMatrix);

    //! @brief return part of the nonlocal matrix
    virtual NuTo::ConstitutiveTangentLocal<6,1>& GetSubMatrix_6x1(int rSubMatrix);

    //! @brief return part of the nonlocal matrix
    virtual NuTo::ConstitutiveTangentLocal<6,6>& GetSubMatrix_6x6(int rSubMatrix);

    //! @brief return number of nonlocal matrices (one for each nonlocal integration point)
    virtual int GetNumSubMatrices()const;

    //! @brief set if a nonlocal matrix is actually local (because of the loading regime)
    virtual void SetLocalSolution(bool rLocalSolution);

    //! @brief return if a nonlocal matrix is actually local (because of the loading regime)
    virtual bool GetLocalSolution()const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
};
}
#endif // ConstitutiveOutputBase_H_
