#pragma once

namespace NuTo
{
namespace Material
{
//! Common material parameters for softening materials
//! @remark For simplicity, we use the common short symbols here. For readability, we deviate from the naming convention
//! (m... ) for members.
struct Softening
{
    //! Young's modulus [pressure]
    double E;

    //! Poisson's ratio [-]
    double nu;

    //! tensile strength [pressure]
    double ft;

    //! compressive strength [pressure]
    double fc;

    //! fracture energy parameter [energy per area / pressure times length ]
    double gf;

    //! nonlocal radius c = l^2 [length]
    double c;

    //! residual strength [pressure]
    double fMin;
};

//! Sets the softening material parameters to approximate normal strength concrete
//! @remark units: length = mm, pressure = MPa
inline Softening DefaultConcrete()
{
    Softening m;
    m.E = 30000;
    m.nu = 0.2;
    m.ft = 4;
    m.fc = 40;
    m.gf = 0.1;
    m.c = 1;
    m.fMin = 0.01 * m.ft;
    return m;
}

} /* NuTo */
} /* Material */
