// $Id$

#pragma once

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
        BOND_STRESS,                    //!< visualize bond stress
        SLIP,                           //!< visualize slip (relative displacement)
        DISPLACEMENTS,                  //!< visualize displacements
        ENGINEERING_STRESS,             //!< visualize engineering stress tensor
        ENGINEERING_STRAIN,             //!< visualize engineering strain tensor
        SHRINKAGE_STRAIN,               //!< visualize shrinkage strain tensor
        THERMAL_STRAIN,                 //!< visualize thermal strain tensor
        NONLOCAL_WEIGHT,                //!< visualize nonlocal weights
        NONLOCAL_EQ_STRAIN,             //!< visualize nonlocal equivalent strains
        LOCAL_EQ_STRAIN,                //!< visualize local equivalent strains
		TOTAL_INELASTIC_EQ_STRAIN,		//!< visualize inelastic equivalent strain
        DAMAGE,                         //!< visualize damage
        ENGINEERING_PLASTIC_STRAIN,     //!< visualize engineering plastic strain
        CONSTITUTIVE,                   //!< visualize constitutive id
        SECTION,                        //!< visualize section id
        ELEMENT,                        //!< visualize element id
        CRACK,                          //!< visualize crack id
        PRINCIPAL_ENGINEERING_STRESS,   //!< visualize principal stresses
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
        HEAT_FLUX,                      //!< visualize heat flux
        RELATIVE_HUMIDITY,              //!< visualize relative humidity
        WATER_VOLUME_FRACTION           //!< visualize water volume fraction
    };

    enum eVisualizationType
    {
        VORONOI_CELL,                   //!< Decomposes the element into smaller cells
        EXTRAPOLATION_TO_NODES,         //!< Extrapolates integration point data to nodes
        POINTS,                         //!< Visualize integration point data as vertex elements
    };

    //! @brief ... export to Vtk datafile
    //! @param rFilename ... filename
    void ExportVtkDataFile(const std::string& rFilename) const;
};

}

