
#include <vector>
#include <functional>

#if defined HAVE_PARDISO
    #include <nuto/math/SparseDirectSolverPardiso.h>
#elif defined HAVE_MUMPS
    #include <nuto/math/SparseDirectSolverMUMPS.h>
#else
    std::cout << "Solver not available - can't solve system of equations " << std::endl;
#endif

#include <boost-1_55/boost/progress.hpp>

#include <nuto/mechanics/structures/unstructured/Structure.h>
#include <nuto/mechanics/constitutive/moistureTransport/ConstitutiveStaticDataMoistureTransport.h>

#include <nuto/metamodel/PolynomialLeastSquaresFitting.h>
#include "nuto/mechanics/timeIntegration/CrankNicolsonEvaluate.h"


class MeshGenerator
{
    static std::vector<int> GetNodeCoordinatesCuboid(NuTo::Structure &rStructure,
                                                     NuTo::FullVector<int,3> rNumElements,
                                                     NuTo::FullVector<double,3> rLength);

    static std::vector<int> GetNodeCoordinatesCuboidMapped(NuTo::Structure &rStructure,
                                                           NuTo::FullVector<int,3> rNumElements,
                                                           NuTo::FullVector<double,3> rLength,
                                                           std::function<NuTo::FullVector<double,Eigen::Dynamic>(double, double, double)> rMappingFunction);

public:
    static void MeshCubuid(NuTo::Structure &rStructure,
                           int rSection,
                           int rConstitutiveLaw,
                           int rInterpolationtype,
                           NuTo::FullVector<int,3> rNumElements,
                           NuTo::FullVector<double,3> rLength);

    static void MeshCylinder(NuTo::Structure &rStructure,
                             int rSection,
                             int rConstitutiveLaw,
                             int rInterpolationtype,
                             NuTo::FullVector<int,3> rNumElements,
                             NuTo::FullVector<double,3> rLength);
};


std::vector<int> MeshGenerator::GetNodeCoordinatesCuboid(NuTo::Structure &rStructure,
                                                         NuTo::FullVector<int, 3> rNumElements,
                                                         NuTo::FullVector<double, 3> rLength)
{
    assert(rNumElements(0)>0);
    assert(rNumElements(1)>0);
    assert(rNumElements(2)>0);
    assert(rLength(0)>0.0);
    assert(rLength(1)>0.0);
    assert(rLength(2)>0.0);

    std::vector<int> NodeIDs;


    double delta_x = rLength(0)/ static_cast<double>(rNumElements(0));
    double delta_y = rLength(1)/ static_cast<double>(rNumElements(1));
    double delta_z = rLength(2)/ static_cast<double>(rNumElements(2));

    int NumNodesX = rNumElements(0)+1;
    int NumNodesY = rNumElements(1)+1;
    int NumNodesZ = rNumElements(2)+1;

    NodeIDs.resize(NumNodesX * NumNodesY * NumNodesZ);

    int NodeNum = 0;

    NuTo::FullVector<double, Eigen::Dynamic> Coordinates(3);

    for(int z_count = 0; z_count < NumNodesZ; z_count++)
    {
        for(int y_count = 0; y_count < NumNodesY; y_count++)
        {
            for(int x_count = 0; x_count < NumNodesX; x_count++)
            {
                NodeNum = x_count + y_count* NumNodesX + z_count * NumNodesX * NumNodesY;

                Coordinates = {{x_count * delta_x,
                                y_count * delta_y,
                                z_count * delta_z}};

                NodeIDs[NodeNum] = rStructure.NodeCreate(Coordinates);
            }
        }
    }

    return NodeIDs;
}

std::vector<int> MeshGenerator::GetNodeCoordinatesCuboidMapped(NuTo::Structure &rStructure,
                                                               NuTo::FullVector<int, 3> rNumElements,
                                                               NuTo::FullVector<double, 3> rLength,
                                                               std::function<NuTo::FullVector<double,Eigen::Dynamic>(double, double, double)> rMappingFunction)
{
    assert(rNumElements(0)>0);
    assert(rNumElements(1)>0);
    assert(rNumElements(2)>0);
    assert(rLength(0)>0.0);
    assert(rLength(1)>0.0);
    assert(rLength(2)>0.0);

    std::vector<int> NodeIDs;


    double delta_x = 2.0/ static_cast<double>(rNumElements(0));
    double delta_y = 2.0/ static_cast<double>(rNumElements(1));
    double delta_z = 2.0/ static_cast<double>(rNumElements(2));

    int NumNodesX = rNumElements(0)+1;
    int NumNodesY = rNumElements(1)+1;
    int NumNodesZ = rNumElements(2)+1;

    NodeIDs.resize(NumNodesX * NumNodesY * NumNodesZ);

    std::cout << std::endl << "Creating Nodes of Mesh" << std::endl << std::endl;

    boost::progress_display show_progress( NodeIDs.size() );


    int NodeNum = 0;

    NuTo::FullVector<double, Eigen::Dynamic> Coordinates(3);

    for(int z_count = 0; z_count < NumNodesZ; z_count++)
    {
        for(int y_count = 0; y_count < NumNodesY; y_count++)
        {
            for(int x_count = 0; x_count < NumNodesX; x_count++)
            {
                NodeNum = x_count + y_count* NumNodesX + z_count * NumNodesX * NumNodesY;

                Coordinates = rMappingFunction(static_cast<double>(x_count) * delta_x - 1.0,
                                               static_cast<double>(y_count) * delta_y - 1.0,
                                               static_cast<double>(z_count) * delta_z - 1.0);

                NodeIDs[NodeNum] = rStructure.NodeCreate(Coordinates);

                ++show_progress;
            }
        }
    }

    return NodeIDs;
}

void MeshGenerator::MeshCubuid(NuTo::Structure& rStructure,
                               int rSection,
                               int rConstitutiveLaw,
                               int rInterpolationtype,
                               NuTo::FullVector<int, 3> rNumElements,
                               NuTo::FullVector<double, 3> rLength)
{


    std::vector<int> NodeIDs = MeshGenerator::GetNodeCoordinatesCuboid(rStructure,
                                                                       rNumElements,
                                                                       rLength);

    NuTo::FullVector<int,Eigen::Dynamic> Nodes(8);

    int NumNodesX = rNumElements(0)+1;
    int NumNodesY = rNumElements(1)+1;


    std::cout << std::endl << "Creating Elements of Mesh" << std::endl << std::endl;

    boost::progress_display show_progress( rNumElements(0) * rNumElements(1) * rNumElements(2) );


    for(int z_count = 0; z_count < rNumElements(2); z_count++)
    {
        for(int y_count = 0; y_count < rNumElements(1); y_count++)
        {
            for(int x_count = 0; x_count < rNumElements(0); x_count++)
            {
                Nodes(0) = NodeIDs[x_count       +  y_count       * NumNodesX   +  z_count       * NumNodesX * NumNodesY];
                Nodes(1) = NodeIDs[x_count + 1   +  y_count       * NumNodesX   +  z_count       * NumNodesX * NumNodesY];
                Nodes(2) = NodeIDs[x_count + 1   + (y_count + 1)  * NumNodesX   +  z_count       * NumNodesX * NumNodesY];
                Nodes(3) = NodeIDs[x_count       + (y_count + 1)  * NumNodesX   +  z_count       * NumNodesX * NumNodesY];
                Nodes(4) = NodeIDs[x_count       + y_count        * NumNodesX   + (z_count + 1)  * NumNodesX * NumNodesY];
                Nodes(5) = NodeIDs[x_count + 1   + y_count        * NumNodesX   + (z_count + 1)  * NumNodesX * NumNodesY];
                Nodes(6) = NodeIDs[x_count + 1   + (y_count + 1)  * NumNodesX   + (z_count + 1)  * NumNodesX * NumNodesY];
                Nodes(7) = NodeIDs[x_count       + (y_count + 1)  * NumNodesX   + (z_count + 1)  * NumNodesX * NumNodesY];

                int elementNumber = rStructure.ElementCreate(rInterpolationtype, Nodes, "CONSTITUTIVELAWIP","STATICDATA");

                // set element section
                rStructure.ElementSetSection(elementNumber,rSection);

                rStructure.ElementSetConstitutiveLaw(elementNumber,rConstitutiveLaw);
                ++show_progress;
            }
        }
    }


}

void MeshGenerator::MeshCylinder(NuTo::Structure &rStructure,
                                 int rSection,
                                 int rConstitutiveLaw,
                                 int rInterpolationtype,
                                 NuTo::FullVector<int, 3> rNumElements,
                                 NuTo::FullVector<double, 3> rLength)
{


    auto MappingFunction = [&rLength](double rX, double rY, double rZ) -> NuTo::FullVector<double, Eigen::Dynamic>
                            {
                                NuTo::FullVector<double,Eigen::Dynamic> CoordVec(3);
                                CoordVec = {{rX * sqrt(1 - (rY * rY) / 2.0 ) * rLength(0) / 2.0,
                                             rY * sqrt(1 - (rX * rX) / 2.0 ) * rLength(1) / 2.0,
                                             rZ * rLength(2) / 2.0}};
                                return CoordVec;
                            };

    std::vector<int> NodeIDs = MeshGenerator::GetNodeCoordinatesCuboidMapped(rStructure,
                                                                             rNumElements,
                                                                             rLength,
                                                                             MappingFunction);

    NuTo::FullVector<int,Eigen::Dynamic> Nodes(8);

    int NumNodesX = rNumElements(0)+1;
    int NumNodesY = rNumElements(1)+1;


    std::cout << std::endl << "Creating Elements of Mesh" << std::endl << std::endl;

    boost::progress_display show_progress( rNumElements(0) * rNumElements(1) * rNumElements(2) );

    for(int z_count = 0; z_count < rNumElements(2); z_count++)
    {
        for(int y_count = 0; y_count < rNumElements(1); y_count++)
        {
            for(int x_count = 0; x_count < rNumElements(0); x_count++)
            {
                Nodes(0) = NodeIDs[x_count       +  y_count       * NumNodesX   +  z_count       * NumNodesX * NumNodesY];
                Nodes(1) = NodeIDs[x_count + 1   +  y_count       * NumNodesX   +  z_count       * NumNodesX * NumNodesY];
                Nodes(2) = NodeIDs[x_count + 1   + (y_count + 1)  * NumNodesX   +  z_count       * NumNodesX * NumNodesY];
                Nodes(3) = NodeIDs[x_count       + (y_count + 1)  * NumNodesX   +  z_count       * NumNodesX * NumNodesY];
                Nodes(4) = NodeIDs[x_count       + y_count        * NumNodesX   + (z_count + 1)  * NumNodesX * NumNodesY];
                Nodes(5) = NodeIDs[x_count + 1   + y_count        * NumNodesX   + (z_count + 1)  * NumNodesX * NumNodesY];
                Nodes(6) = NodeIDs[x_count + 1   + (y_count + 1)  * NumNodesX   + (z_count + 1)  * NumNodesX * NumNodesY];
                Nodes(7) = NodeIDs[x_count       + (y_count + 1)  * NumNodesX   + (z_count + 1)  * NumNodesX * NumNodesY];

                int elementNumber = rStructure.ElementCreate(rInterpolationtype, Nodes, "CONSTITUTIVELAWIP","STATICDATA");

                // set element section
                rStructure.ElementSetSection(elementNumber,rSection);

                rStructure.ElementSetConstitutiveLaw(elementNumber,rConstitutiveLaw);

                ++show_progress;

            }
        }
    }

}




int main()
{
    try
    {

        /*------------------------------------*\
        |                Setup                 |
        \*------------------------------------*/


        NuTo::FullVector<int,3>     NumElements = {{    2,
                                                        2,
                                                        6}};


        NuTo::FullVector<double,3>  Length      = {{    0.1,
                                                        0.1,
                                                        0.3 }};

        bool            EnableSorptionHysteresis            = false;
        bool            EnableModiefiedTangentialStiffness  = false;

        double          delta_t                         = 1.0/1.0 *     1.0 * 24.0 * 60.0 * 60.0;
        double          t_final                         = 1.0/1.0 *     3.0 * 24.0 * 60.0 * 60.0;
        double          BC_TransitionTime               =                     24.0 * 60.0 * 60.0;


        // initial node values
        double          InitialRelativeHumidity         =    0.95;
        double          InitialWaterVolumeFraction      =    0.03;


        // constitutive law values
        double          MassExchangeRate                =    3.42e-7    ;
        double          Porosity                        =    0.25      ;
        double          VaporPhaseDiffusionCoefficient  =    3.9e-10     ;
        double          VaporPhaseDiffusionExponent     =    1.0        ;
        double          VaporPhaseSaturationDensity     =    0.0173     ;
        double          WaterPhaseDensity               =  999.97       ;
        double          WaterPhaseDiffusionCoefficient  =    1.17e-7    ;
        double          WaterPhaseDiffusionExponent     =    2.0        ;

        // Boundary Condition Values
        double          BC_RelativeHumidity             =    0.45;
        double          BC_Surface_Moisture_Transfer_RH =    1.0e-10 * delta_t;
        double          BC_Surface_Moisture_Transfer_WVF=    1.0e-7 * delta_t;
        bool            SorptionHistoryDesorption       =    true;

        // sorption hysteresis
        double          Ka                              =    0.26       ;
        double          Kd                              =    0.56       ;
        NuTo::PolynomialLeastSquaresFitting AdsorptionFit;
        NuTo::PolynomialLeastSquaresFitting DesorptionFit;


        double          MaxResidual                     =    1.0e-18;



        /*------------------------------------*\
        |         Fit Sorption Curves          |
        \*------------------------------------*/


        NuTo::FullVector <double,Eigen::Dynamic> x_Values_Ad({0.1, 0.2 ,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9});
        NuTo::FullVector <double,Eigen::Dynamic> y_Values_Ad({0.017, 0.03, 0.04, 0.048, 0.056, 0.066, 0.077, 0.092,0.114});

        NuTo::FullVector <double,Eigen::Dynamic> x_Values_De({0.1, 0.2 ,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9});
        NuTo::FullVector <double,Eigen::Dynamic> y_Values_De({0.022, 0.039, 0.052, 0.062, 0.072, 0.083, 0.097, 0.118,0.145});

        AdsorptionFit.SetSupportPoints(1,1,x_Values_Ad.Trans(),y_Values_Ad.Trans());
        AdsorptionFit.SetDegree(3);
        AdsorptionFit.AddBoundaryCondition(0.0,0.0);
        AdsorptionFit.AddBoundaryCondition(1.0,0.141);
        AdsorptionFit.BuildDerived();

        DesorptionFit.SetSupportPoints(1,1,x_Values_De.Trans(),y_Values_De.Trans());
        DesorptionFit.SetDegree(3);
        DesorptionFit.AddBoundaryCondition(0.0,0.0);
        DesorptionFit.AddBoundaryCondition(1.0,0.182);
        DesorptionFit.BuildDerived();




        /*------------------------------------*\
        |           Create Structure           |
        \*------------------------------------*/


        NuTo::Structure myStructure(3);
        myStructure.SetNumTimeDerivatives(2);

        // disable output of calculation times
        myStructure.SetShowTime(false);
        myStructure.UseMaximumIndependentSets(true);





        /*------------------------------------*\
        |            Create Section            |
        \*------------------------------------*/

        int mySection = myStructure.SectionCreate("Volume");
        myStructure.ElementTotalSetSection(mySection);





        /*------------------------------------*\
        |       Set Interpolation Type         |
        \*------------------------------------*/

        int myInterpolationType = myStructure.InterpolationTypeCreate("Brick3d");
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::RELATIVEHUMIDITY, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::WATERVOLUMEFRACTION, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);




        /*------------------------------------*\
        |       Create Constitutive Law        |
        \*------------------------------------*/

        int myConstitutiveLaw = myStructure.ConstitutiveLawCreate("MoistureTransport");

        // set variables
        myStructure.ConstitutiveLawSetParameterBool    (myConstitutiveLaw,"ENABLE_MODIFIED_TANGENTIAL_STIFFNESS",EnableModiefiedTangentialStiffness);      // sets whether modified tangential stiffness should be used or not
        myStructure.ConstitutiveLawSetParameterBool    (myConstitutiveLaw,"enable_sorption_hysteresis",EnableSorptionHysteresis);                          // sets whether sorption hysteresis should be used or not


        myStructure.ConstitutiveLawSetParameterDouble  (myConstitutiveLaw,"boundary_TRANSPORT_CONSTANT_GAS_PHASE",BC_Surface_Moisture_Transfer_RH);        // set water phase density
        myStructure.ConstitutiveLawSetParameterDouble  (myConstitutiveLaw,"BOUNDARY_TRANSPORT_CONSTANT_WATER_PHASE",BC_Surface_Moisture_Transfer_WVF);     // set water phase density
        myStructure.ConstitutiveLawSetParameterDouble  (myConstitutiveLaw,"DENSITY_WATER_PHASE",WaterPhaseDensity);                                        // set water phase density
        myStructure.ConstitutiveLawSetParameterDouble  (myConstitutiveLaw,"DIFFUSION_CONSTANT_GAS_PHASE",VaporPhaseDiffusionCoefficient);                  // set vapor phase diffusion coefficient
        myStructure.ConstitutiveLawSetParameterDouble  (myConstitutiveLaw,"DIFFUSION_CONSTANT_WATER_PHASE",WaterPhaseDiffusionCoefficient);                // set water phase diffusion coefficient
        myStructure.ConstitutiveLawSetParameterDouble  (myConstitutiveLaw,"DIFFUSION_EXPONENT_GAS_PHASE",VaporPhaseDiffusionExponent);                     // set vapor phase diffusion exponent
        myStructure.ConstitutiveLawSetParameterDouble  (myConstitutiveLaw,"DIFFUSION_EXPONENT_WATER_PHASE",WaterPhaseDiffusionExponent);                   // set water phase diffusion exponent
        myStructure.ConstitutiveLawSetParameterDouble  (myConstitutiveLaw,"GRADIENT_CORRECTION_ADSORPTION_DESORPTION",Kd);                                 // set gradient correction when changing from adsorption to desorption
        myStructure.ConstitutiveLawSetParameterDouble  (myConstitutiveLaw,"GRADIENT_CORRECTION_DESORPTION_ADSORPTION",Ka);                                 // set gradient correction when changing from desorption to adsorption
        myStructure.ConstitutiveLawSetParameterDouble  (myConstitutiveLaw,"MASS_EXCHANGE_RATE",MassExchangeRate);                                          // set mass exchange rate
        myStructure.ConstitutiveLawSetParameterDouble  (myConstitutiveLaw,"POROSITY",Porosity);                                                            // set porosity
        myStructure.ConstitutiveLawSetParameterDouble  (myConstitutiveLaw,"SATURATION_DENSITY_GAS_PHASE",VaporPhaseSaturationDensity);                     // set vapor phase saturation density

        myStructure.ConstitutiveLawSetParameterFullVectorDouble    (myConstitutiveLaw,"polynomial_COEFFICIENTS_ADSORPTION",AdsorptionFit.GetPolynomialCoefficients());               // set adsorption coefficients
        myStructure.ConstitutiveLawSetParameterFullVectorDouble    (myConstitutiveLaw,"POLYNOMIAL_COEFFICIENTS_DESORPTION",DesorptionFit.GetPolynomialCoefficients());               // set desorption coefficients


        // Calculate equilibrium water volume fraction
        InitialWaterVolumeFraction   = myStructure.ConstitutiveLawGetEquilibriumWaterVolumeFraction(myConstitutiveLaw,InitialRelativeHumidity,myStructure.ConstitutiveLawGetParameterFullVectorDouble(myConstitutiveLaw,NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));





        /*------------------------------------*\
        |             Create Mesh              |
        \*------------------------------------*/

        MeshGenerator::MeshCylinder(myStructure,mySection,myConstitutiveLaw,myInterpolationType,NumElements,Length);
        //MeshGenerator::MeshCubuid(myStructure,mySection,myConstitutiveLaw,myInterpolationType,NumElements,Length);


        myStructure.ElementTotalConvertToInterpolationType(1.0e-12,0.001);


        //return 0;

        /*------------------------------------*\
        |       Set nodal start values         |
        \*------------------------------------*/


        unsigned int NNodes = myStructure.GetNumNodes();

        for (unsigned int i=0; i<NNodes; i++)
        {
/*
            double NodeX = myStructure.NodeGetNodePtr(i)->GetCoordinates()[0] / Length(0);
            double NodeY = myStructure.NodeGetNodePtr(i)->GetCoordinates()[1] / Length(1);
            double XVal = 0.5 * sin(NodeX * 3.14) + 0.5 * NodeX * NodeX;
            double YVal = (cos((NodeY + 0.5)* 3.14) + 1.0) * (0.5 +0.5*NodeY);
            double NodeMultiplier = 0.25 +  XVal * 0.75 * YVal;
*/
            double NodeMultiplier = 1.0;
            if(myStructure.NodeGetNodePtr(i)->GetNumRelativeHumidity() != 0)
            {
                myStructure.NodeGetNodePtr(i)->SetRelativeHumidity(0,InitialRelativeHumidity*NodeMultiplier);
            }
            if(myStructure.NodeGetNodePtr(i)->GetNumWaterVolumeFraction() != 0)
            {
                myStructure.NodeGetNodePtr(i)->SetWaterVolumeFraction(0,InitialWaterVolumeFraction*NodeMultiplier);
            }
        }



        /*------------------------------------*\
        |      Create boundary elements        |
        \*------------------------------------*/

        int nodeGroupBoundary = myStructure.GroupCreate("NODES");


        auto GetBoundaryNodesLambda = [Length](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNumCoordinates()>0)
                                    {
                                        double x = rNodePtr->GetCoordinate(0);
                                        double y = rNodePtr->GetCoordinate(1);
                                        double z = rNodePtr->GetCoordinate(2);

                                        double r = sqrt(x*x + y*y);

                                        if ((z >= -Length(2) / 2.0 -Tol   && z <= -Length(2) / 2.0 + Tol) ||
                                            (z >=  Length(2) / 2.0 -Tol   && z <=  Length(2) / 2.0 + Tol) ||
                                            (r >=  Length(0) / 2.0 -Tol   && r <=  Length(0) / 2.0 + Tol) )
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };

        /*
        auto GetBoundaryNodesLambda = [Length](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNumCoordinates()>0)
                                    {
                                        double z = rNodePtr->GetCoordinate(2);

                                        if ((z >=             -Tol   && z <=            + Tol) ||
                                            (z >=  Length(2)  -Tol   && z <=  Length(2) + Tol) )
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };
        */

        myStructure.GroupAddNodeFunction(nodeGroupBoundary,GetBoundaryNodesLambda);

        // Group all elements with boundary nodes
        int elemGroupBoundary = myStructure.GroupCreate("ELEMENTS");
        myStructure.GroupAddElementsFromNodes(elemGroupBoundary, nodeGroupBoundary, false);

        // additional node that stores the environmental relative humidity
        int BoundaryNodeID = myStructure.NodeCreateDOFs("RELATIVEHUMIDITY");

        myStructure.BoundaryElementsCreate(elemGroupBoundary,nodeGroupBoundary,myStructure.NodeGetNodePtr(BoundaryNodeID));



        /*------------------------------------*\
        |           Set static data            |
        \*------------------------------------*/


        // Loop over all integration points
        for (int i=0; i<myStructure.GetNumElements(); i++)
        {
            for (int theIP=0; theIP< myStructure.ElementGetElementPtr(i)->GetNumIntegrationPoints(); theIP++)
            {
                NuTo::ConstitutiveStaticDataMoistureTransport *StaticData = myStructure.ElementGetElementPtr(i)->GetStaticData(theIP)->AsMoistureTransport();
                StaticData->SetLastSorptionCoeff(myStructure.ConstitutiveLawGetParameterFullVectorDouble(myConstitutiveLaw,NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
                StaticData->SetActualSorptionCoeff(myStructure.ConstitutiveLawGetParameterFullVectorDouble(myConstitutiveLaw,NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
                StaticData->SetLastRelHumValue(InitialRelativeHumidity);
                StaticData ->SetSorptionHistoryDesorption(SorptionHistoryDesorption);
            }
        }



        /*------------------------------------*\
        |           Set constraints            |
        \*------------------------------------*/

        auto RelativeHumidityFunc = [BC_RelativeHumidity,
                                     BC_TransitionTime,
                                     InitialRelativeHumidity,
                                     t_final]
                                     (double rTime)->double
                                  {

                                        if(rTime == 0.0)
                                        {
                                            return InitialRelativeHumidity;
                                        }
                                        else
                                        {
                                            if(rTime< BC_TransitionTime)
                                            {
                                                return InitialRelativeHumidity - sin(rTime / BC_TransitionTime * 3.14 /2.0) * (InitialRelativeHumidity-BC_RelativeHumidity);
                                            }
                                            {
                                                return BC_RelativeHumidity;
                                            }
                                        }
                                  };

        int BoundaryConstraint = myStructure.ConstraintLinearSetRelativeHumidityNode(BoundaryNodeID,1.0);



        /*------------------------------------*\
        |        Multi processor setup         |
        \*------------------------------------*/

#ifdef _OPENMP
        myStructure.SetNumProcessors(4);
        std::cout << "OpenMP enabled" << std::endl;
#else
        myStructure.SetNumProcessors(1);
#endif
        myStructure.CalculateMaximumIndependentSets();
        myStructure.UseMaximumIndependentSets(true);
        myStructure.NodeBuildGlobalDofs();




        /*------------------------------------*\
        |              Set solver              |
        \*------------------------------------*/

#if defined HAVE_PARDISO
        NuTo::SparseDirectSolverPardiso Solver(4);
#elif defined HAVE_MUMPS
        NuTo::SparseDirectSolverMUMPS Solver;
#else
        std::cout << "Solver not available - can't solve system of equations " << std::endl;
#endif
        Solver.SetShowTime(false);





        /*------------------------------------*\
        |           Time integration           |
        \*------------------------------------*/

        myStructure.AddVisualizationComponentRelativeHumidity();
        myStructure.AddVisualizationComponentWaterVolumeFraction();

        NuTo::CrankNicolsonEvaluate myTimeIntegrationScheme(&myStructure);


        myTimeIntegrationScheme.SetPerformLineSearch(false);
        myTimeIntegrationScheme.SetCheckEquilibriumOnStart(false);
        myTimeIntegrationScheme.SetToleranceForce(MaxResidual);
        myTimeIntegrationScheme.SetMaxNumIterations(40);

        myTimeIntegrationScheme.AddTimeDependentConstraint(BoundaryConstraint,  RelativeHumidityFunc);


        //set result directory
        bool deleteResultDirectoryFirst(true);
        myTimeIntegrationScheme.SetResultDirectory("./ResultsMoistureTransport3D",deleteResultDirectoryFirst);


        myTimeIntegrationScheme.SetTimeStep(delta_t);
        myTimeIntegrationScheme.SetMinTimeStepPlot(t_final);
        myTimeIntegrationScheme.Solve(t_final);



    }
    catch(NuTo::Exception e)
    {
        std::cout << e.ErrorMessage() << std::endl;
        return 1;
    }





}



