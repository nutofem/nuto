// $Id$

#ifndef VISUALIZEBASE_H_
#define VISUALIZEBASE_H_

#include <string>

namespace NuTo
{

//! @brief ... base class for visualize objects
//! @author Stefan Eckardt, ISM
//! @date November 2009
class VisualizeBase
{
public:
    enum eVisualizeWhat
    {
        DISPLACEMENTS,					//!< visualize displacements
        ENGINEERING_STRESS,				//!< visualize engineering stress tensor
        ENGINEERING_STRAIN,				//!< visualize engineering strain tensor
        NONLOCAL_WEIGHT,                //!< visualize nonlocal weights
        NONLOCAL_EQ_STRAIN,             //!< visualize nonlocal equivalent strains
        DAMAGE,							//!< visualize damage
        ENGINEERING_PLASTIC_STRAIN,		//!< visualize engineering plastic strain
        CONSTITUTIVE,					//!< visualize constitutive id
        SECTION,						//!< visualize section id
        ELEMENT,						//!< visualize element id
        CRACK,							//!< visualize crack id
        PRINCIPAL_ENGINEERING_STRESS,	//!< visualize principal stresses
        PARTICLE_RADIUS,                //!< visualize radius of particles/nodes
        LATTICE_STRAIN,                 //!< visualize strain of lattice models
        LATTICE_STRESS,                 //!< visualize stress of lattice models
        LATTICE_PLASTIC_STRAIN,         //!< visualize plastic strain of lattice models
        ROTATION,                       //!< visualize rotations
        VELOCITY,                       //!< visualize velocity
        ANGULAR_VELOCITY,               //!< visualize angular velocity
        ACCELERATION,                   //!< visualize acceleration
        ANGULAR_ACCELERATION,           //!< visualize angular acceleration
        TEMPERATURE,                    //!< visualize temperature
        HEAT_FLUX                       //!< visualize heat flux
    };

    //! @brief ... export to Vtk datafile
    //! @param rFilename ... filename
    void ExportVtkDataFile(const std::string& rFilename) const;
};

}

#endif // VISUALIZEBASE_H_ 
