// $Id$

#ifndef CONSTITUTIVESTATICDATABASE_H_
#define CONSTITUTIVESTATICDATABASE_H_

#include <string>
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/elements/ElementEnum.h"

#ifdef ENABLE_VISUALIZE
#include <boost/ptr_container/ptr_list.hpp>
#endif // ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

//! @brief ... base class, storing the static data (history variables) of a constitutive relationship
//! @author Stefan Eckardt, ISM
//! @date November 2009
namespace NuTo
{
class ConstitutiveStaticDataLatticeConcrete2D;
class ConstitutiveStaticDataGradientDamagePlasticity1D;
class ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain;
class ConstitutiveStaticDataMisesPlasticity3D;
class ConstitutiveStaticDataMultiscale2DPlaneStrain;
class ElementBase;
class VisualizeUnstructuredGrid;
class VisualizeComponentBase;
template<class T> class  FullMatrix;

class ConstitutiveStaticDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief constructor
    ConstitutiveStaticDataBase()
    {}

    //! @brief destructor (virtual, in order to make the class a polymorphic type)
    virtual ~ConstitutiveStaticDataBase()
    {};

    //! @brief constructor
    virtual ConstitutiveStaticDataBase* Clone()const=0;

    NuTo::ConstitutiveStaticDataBase& operator= (NuTo::ConstitutiveStaticDataBase const& rOther);

    //! @brief check, if the static data is compatible with a given element and a given constitutive model
    virtual bool CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType constitutiveType, NuTo::Element::eElementType elementType)const=0;

    //! @brief sets the fine scale model (deserialization from a binary file)
    virtual void SetFineScaleModel(std::string rFileName, double rMacroLength, double rCenter[2], std::string rIPName);

    //! @brief sets the fine scale model (deserialization from a binary file)
    virtual void SetFineScaleParameter(const std::string& rName, double rParameter);

    //! @brief sets the fine scale model (deserialization from a binary file)
    virtual void SetFineScaleParameter(const std::string& rName, std::string rParameter);

    //!@ brief reinterpret as nonlocal damage2d static data
    virtual ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain* AsNonlocalDamagePlasticity2DPlaneStrain();

    //!@ brief reinterpret as nonlocal damage2d static data
    virtual const ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain* AsNonlocalDamagePlasticity2DPlaneStrain()const;

    //!@ brief reinterpret as gradient damage 1d static data
    virtual ConstitutiveStaticDataGradientDamagePlasticity1D* AsGradientDamagePlasticity1D();

    //!@ brief reinterpret as gradient damage1d static data
    virtual const ConstitutiveStaticDataGradientDamagePlasticity1D* AsGradientDamagePlasticity1D()const;

    //!@ brief reinterpret as lattice concrete 2D static data
    virtual NuTo::ConstitutiveStaticDataLatticeConcrete2D* AsConstitutiveStaticDataLatticeConcrete2D();

    //!@ brief reinterpret as lattice concrete 2D static data
    virtual const NuTo::ConstitutiveStaticDataLatticeConcrete2D* AsConstitutiveStaticDataLatticeConcrete2D()const;

    //!@ brief reinterpret as nonlocal damage2d static data
    virtual ConstitutiveStaticDataMultiscale2DPlaneStrain* AsMultiscale2DPlaneStrain();

    //!@ brief reinterpret as nonlocal damage2d static data
    virtual const ConstitutiveStaticDataMultiscale2DPlaneStrain* AsMultiscale2DPlaneStrain()const;

    //!@ brief reinterpret
    virtual ConstitutiveStaticDataMisesPlasticity3D* AsConstitutiveStaticDataMisesPlasticity3D();

    //!@ brief reinterpret
    virtual const ConstitutiveStaticDataMisesPlasticity3D* AsConstitutiveStaticDataMisesPlasticity3D() const;

#ifdef ENABLE_VISUALIZE
    //Visualize for all integration points the fine scale structure
    virtual void VisualizeIpMultiscale(VisualizeUnstructuredGrid& rVisualize,
    		const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat, bool rVisualizeDamage)const;
#endif


#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
};

}

#endif // CONSTITUTIVESTATICDATABASE_H_ 
