#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "nuto/mechanics/constitutive/laws/HeatConduction.h"

void TestTangentMatrix(double conductivity)
{
    // create constitutive law
    NuTo::HeatConduction heat_conduction;
    heat_conduction.SetParameterDouble(
            NuTo::Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY, conductivity);
    NuTo::ConstitutiveInputMap input_map;
    NuTo::ConstitutiveOutputMap output_map;
    NuTo::ConstitutiveMatrix<3, 3> tangent;
    output_map[NuTo::Constitutive::Output::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT]
            = &tangent;

    // evaluate tangent
    NuTo::ElementBase* some_element = nullptr;
    heat_conduction.Evaluate<3>(some_element, 0, input_map, output_map);

    // compare to expected output
    Eigen::Matrix<double, 3, 3> matrix;
    matrix << conductivity, 0, 0,
              0, conductivity, 0,
              0, 0, conductivity;
    if(matrix != tangent)
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Tangent matrix does not have expected value!");
}

int main()
{
    try
    {
        TestTangentMatrix(5.0);
        try
        {
            TestTangentMatrix(-1.0);
        }
        catch (NuTo::MechanicsException& e)
        {
            std::string message = e.what();
            // I'm expecting it to fail, throwing "THERMAL_CONDUCTIVITY must be > 0"
            if(message.find("THERMAL_CONDUCTIVITY must be > 0") == std::string::npos)
                throw e;
        }
    }
    catch (NuTo::Exception& e)
    {
        std::cout << e.what();
        return 1;
    }

    return 0;
}
