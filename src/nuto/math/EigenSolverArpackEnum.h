#pragma once

namespace NuTo
{
namespace EIGEN_SOLVER_ARPACK
{
enum class eDriver
{
	DSDRV1=0,            //!< regular eigenvalue problem (see documentation of arpack (section 3.5 page 31)
	DSDRV2,
	DSDRV3,
	DSDRV4,
	DSDRV5,
	DSDRV6,
	DNDRV1,              //!< regular eigenvalue problem (see documentation of arpack (section 3.5 page 31)
	DNDRV2,
	DNDRV3,
	DNDRV4,
	DNDRV5,
	DNDRV6
};

enum class eWhich
{
    LA=0,            //!< regular eigenvalue problem (see documentation of arpack (section 3.5 page 31)
    SA,
    LM,
    SM,
    BE,
    LR,
    SR,
    LI,
    SI
};


}//EigenSolverArpack
}//Nuto
