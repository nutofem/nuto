#include "BoostUnitTest.h"

#include "mechanics/constitutive/laws/LinearElasticAnisotropic.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveMatrix.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"

using namespace NuTo;
BOOST_AUTO_TEST_CASE(tangent_matrix)
{
    // set up constitutive law
    LinearElasticAnisotropic anisotropicLaw;

    Eigen::MatrixXd stiffnessTensor(6,6);
    stiffnessTensor << 1 , 1 , 1 , 0   , 0   , 0 ,
                       1 , 1 , 1 , 0   , 0   , 0 ,
                       1 , 1 , 1 , 0   , 0   , 0 ,
                       0 , 0 , 0 , 1   , 0   , 0 ,
                       0 , 0 , 0 , 0   , 1   , 0 ,
                       0 , 0 , 0 , 0   , 0   , 1 ;

    anisotropicLaw.SetParameterMatrixDouble(Constitutive::eConstitutiveParameter::STIFFNESS,stiffnessTensor);

    ConstitutiveInputMap input_map;
    ConstitutiveOutputMap output_map;
    output_map[Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN] =
            ConstitutiveIOBase::makeConstitutiveIO<3>(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN);

    // evaluate tangent
    anisotropicLaw.Evaluate<3>(input_map, output_map);

    Eigen::Matrix<double, 6, 6> calculated_tangent =
            *static_cast<ConstitutiveMatrix<6,6>*>(output_map.at(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN).get());
    BOOST_CHECK_EQUAL(calculated_tangent, stiffnessTensor);

}
