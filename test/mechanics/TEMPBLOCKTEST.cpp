
//#include <nuto/mechanics/timeIntegration/ImplicitStandard.h>
#include <nuto/mechanics/constitutive/multiPhysics/ConstitutiveStaticDataMultiPhysics.h>
#include <nuto/mechanics/constitutive/moistureTransport/ConstitutiveStaticDataMoistureTransport.h>
#include "nuto/geometryConcrete/GeometryConcrete.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
//#include "nuto/mechanics/timeIntegration/CrankNicolson.h"
#include "nuto/mechanics/timeIntegration/CrankNicolsonEvaluate.h"
#include <nuto/metamodel/PolynomialLeastSquaresFitting.h>
#include <nuto/mechanics/nodes/NodeBase.h>
#include <boost/filesystem.hpp>




void CreateNodes(NuTo::Structure& rStructure,
                    double rLX,
                    double rLY,
                    unsigned int rNumElementsX,
                    unsigned int rNumElementsY)
{

    double delta_x = rLX/ static_cast<double>(rNumElementsX);
    double delta_y = rLY/ static_cast<double>(rNumElementsY);

    int NumNodesX = rNumElementsX+1;
    int NumNodesY = rNumElementsY+1;

    //NodeIDs.resize(NumNodesX * NumNodesY);

    //int NodeNum = 0;

    NuTo::FullVector<double, Eigen::Dynamic> Coordinates(2);

    for(int y_count = 0; y_count < NumNodesY; y_count++)
    {
        for(int x_count = 0; x_count < NumNodesX; x_count++)
        {
            //NodeNum = x_count + y_count* NumNodesX;

            Coordinates = {{x_count * delta_x,
                            y_count * delta_y}};

            rStructure.NodeCreate(Coordinates);
        }
    }
}

int GetNearestNode(double rTargetX, double rTargetY, NuTo::Structure& rStructure, NuTo::FullVector<int,Eigen::Dynamic> rGroupIdsAll)
{

    NuTo::NodeBase* TmpNode = rStructure.NodeGetNodePtr(rGroupIdsAll[0]);
    int NearestNodeID = rGroupIdsAll[0];
    double TmpX = TmpNode->GetCoordinates()(0) - rTargetX;
    double TmpY = TmpNode->GetCoordinates()(1) - rTargetY;
    double TmpNodeDistance = sqrt(TmpX * TmpX + TmpY * TmpY);
    double NodeDistance = TmpNodeDistance;



    for (int i=1; i<rGroupIdsAll.GetNumRows(); i++)
    {
        TmpNode = rStructure.NodeGetNodePtr(rGroupIdsAll[i]);
        TmpX = TmpNode->GetCoordinates()(0) - rTargetX;
        TmpY = TmpNode->GetCoordinates()(1) - rTargetY;
        TmpNodeDistance = sqrt(TmpX * TmpX + TmpY * TmpY);
        if(TmpNodeDistance<NodeDistance)
        {
            NodeDistance = TmpNodeDistance;
            NearestNodeID = rGroupIdsAll[i];
        }

    }

    return NearestNodeID;
}


bool ElementIsAggregate(NuTo::FullVector<int,Eigen::Dynamic>& rNodes, NuTo::Structure& rStructure)
{

    NuTo::NodeBase* tmpNode = rStructure.NodeGetNodePtr(rNodes[0]);
    Eigen::Matrix<double, 2, 1> Coordinates = tmpNode->GetCoordinates2D();

    float maxX = Coordinates(0);
    float minX = Coordinates(0);
    float maxY = Coordinates(1);
    float minY = Coordinates(1);

    for(unsigned int i=1; i<rNodes.GetNumRows(); ++i)
    {
        tmpNode = rStructure.NodeGetNodePtr(rNodes[i]);
        Coordinates = tmpNode->GetCoordinates2D();
        if(Coordinates(0)<minX)
        {
            minX = Coordinates(0);
        }
        if(Coordinates(0)>maxX)
        {
            maxX = Coordinates(0);
        }
        if(Coordinates(1)<minY)
        {
            minY = Coordinates(1);
        }
        if(Coordinates(1)>maxY)
        {
            maxY = Coordinates(1);
        }
    }


    if(minX >= 0.06 &&
       maxX <= 0.10 &&
       minY >= 0.01 &&
       maxY <= 0.03 )
    {
        //return true;
    }



    return false;
}



int main(int argc, char* argv[])
{

    try
    {
    /* -------------------------------- *\
    |             Variables              |
    \* -------------------------------- */


        double lX = 0.16/1.0;
        double lY = 0.04/1.0;
        double lZ = 0.04/1.0;

        unsigned int NumElementsX = 160;
        unsigned int NumElementsY = 40;


        double          delta_t                         = 1.0/4.0 *     1.0 * 24.0 * 60.0 * 60.0;
        double          delta_t_write                   = delta_t * 4.0;
        double          SimulationTime                  = 1.0/1.0 *     8.0 * 24.0 * 60.0 * 60.0;
        double          BC_TransitionTime               =                     24.0 * 60.0 * 60.0;


        // Linear Elastic
        double DensityMatrix                                  = 1.0;      // N/mm³
        double PoissonRatioMatrix                             = 0.15;
        double YoungsModulusMatrix                            = 30.0 * 1e9;  //N/m²

        double DensityAggregates                              = 1.0;      // N/mm³
        double PoissonRatioAggregates                         = 0.15;
        double YoungsModulusAggregates                        = 60.0 * 1e9;  //N/m²


        // Moisture Transport


        bool            EnableSorptionHysteresis            = false;
        bool            EnableModiefiedTangentialStiffness  = false;
        double          MassExchangeRate                =    3.42e-7    ;
        double          Porosity                        =    0.25      ;
        double          VaporPhaseDiffusionCoefficient  =    3.9e-12    ;
        double          VaporPhaseDiffusionExponent     =    1.0        ;
        double          VaporPhaseSaturationDensity     =    0.0173     ;
        double          WaterPhaseDensity               =  999.97       ;
        double          WaterPhaseDiffusionCoefficient  =    1.17e-7    ;
        double          WaterPhaseDiffusionExponent     =    2.0        ;
        double          InitialRelativeHumidity         =    1.00;
        double          InitialWaterVolumeFraction      =    0.03;

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



        boost::filesystem::path outputPath = std::string(argv[0]) + "Out/";
        boost::filesystem::path resultPath = outputPath.string() + "TMPBLOCKResults/";
        boost::filesystem::create_directory(outputPath);
        boost::filesystem::create_directory(resultPath);

        std::string gmshFile = outputPath.string() + "geometry";

        //std::cout << "Gmsh File:  " << gmshFile << ".msh" << std::endl;
        //std::cout << "Result dir: " << resultPath.string() << std::endl;





        /* -------------------------------- *\
        |         Fit Sorption Curves        |
        \* -------------------------------- */


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


    /* -------------------------------- *\
    |          Create Structure          |
    \* -------------------------------- */

        NuTo::Structure myStructure(2);
        myStructure.SetNumTimeDerivatives(1);
        myStructure.SetNumProcessors(4);
        myStructure.LoggerOpenFile(outputPath.string() + "Log.dat");
        myStructure.LoggerSetQuiet(false);
        myStructure.SetVerboseLevel(10);
        myStructure.SetShowTime(false);





    /* -------------------------------- *\
    |     Create Meso Scale Geometry     |
    \* -------------------------------- */

        CreateNodes(myStructure,
                    lX,
                    lY,
                    NumElementsX,
                    NumElementsY);

//        auto groupIndices = myStructure.ImportFromGmsh(gmshFile+".msh", NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::STATICDATA);


//        int groupMatrix = groupIndices.GetValue(0,0);
//        int groupAggregates = groupIndices.GetValue(1,0);


//        int interpolationMatrix = groupIndices.GetValue(0,1);
//        int interpolationAggregates = groupIndices.GetValue(1,1);





    /* -------------------------------- *\
    |      Set Interpolation Types       |
    \* -------------------------------- */

        // Matrix
        // ------

        int interpolationMatrix = myStructure.InterpolationTypeCreate("Quad2D");
        myStructure.InterpolationTypeAdd(interpolationMatrix,       NuTo::Node::COORDINATES,            NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(interpolationMatrix,       NuTo::Node::DISPLACEMENTS,          NuTo::Interpolation::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(interpolationMatrix,       NuTo::Node::RELATIVEHUMIDITY,       NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(interpolationMatrix,       NuTo::Node::WATERVOLUMEFRACTION,    NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);



        // Aggregates
        // ----------

        int interpolationAggregates = myStructure.InterpolationTypeCreate("Quad2D");
        myStructure.InterpolationTypeAdd(interpolationAggregates,   NuTo::Node::COORDINATES,            NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(interpolationAggregates,   NuTo::Node::DISPLACEMENTS,          NuTo::Interpolation::EQUIDISTANT1);








        /* -------------------------------- *\
        |          Create Elements           |
        \* -------------------------------- */



        int groupMatrix = myStructure.GroupCreate("ELEMENTS");
        int groupAggregates = myStructure.GroupCreate("ELEMENTS");

        NuTo::FullVector<int,Eigen::Dynamic> Nodes(4);

        unsigned int NumNodesX = NumElementsX + 1;
        unsigned int NumNodesY = NumElementsY + 1;
        int elementNumber = 0;

        for (unsigned int iY=0; iY< NumElementsY; iY++)
        {
            for (unsigned int iX=0; iX < NumElementsX; iX++)
            {
                Nodes(0) = iX    +  iY    *NumNodesX;
                Nodes(1) = iX +1 +  iY    *NumNodesX;
                Nodes(2) = iX +1 + (iY+1) *NumNodesX;
                Nodes(3) = iX    + (iY+1) *NumNodesX;

                if(ElementIsAggregate(Nodes,myStructure))
                {
                    elementNumber = myStructure.ElementCreate(interpolationAggregates, Nodes, "CONSTITUTIVELAWIP","STATICDATA");
                    myStructure.GroupAddElement(groupAggregates,elementNumber);
                }
                else
                {
                    elementNumber = myStructure.ElementCreate(interpolationMatrix, Nodes, "CONSTITUTIVELAWIP","STATICDATA");
                    myStructure.GroupAddElement(groupMatrix, elementNumber);
                }
                // set element section
                //myStructure.ElementSetSection(elementNumber,mySection);

                // set element constitutive law
                //myStructure.ElementSetConstitutiveLaw(elementNumber,ConstLawMultiPhysics);
            }
        }

        myStructure.ElementTotalConvertToInterpolationType();



        /* -------------------------------- *\
        |      Create Constitutive Laws      |
        \* -------------------------------- */



        // Matrix
        // ------

        int myConstitutiveLawMatrixLinearElastic        = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::LINEAR_ELASTIC_ENGINEERING_STRESS);
        int myConstitutiveLawMatrixDryingShrinkage      = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::DRYING_SHRINKAGE);
        int myConstitutiveLawMatrixMoistureTransport    = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::MOISTURE_TRANSPORT);
        int myConstitutiveLawMatrixMultiPhysics         = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::MULTI_PHYSICS);
        myStructure.ConstitutiveLawMultiPhysicsAddConstitutiveLaw(myConstitutiveLawMatrixMultiPhysics,myConstitutiveLawMatrixLinearElastic);
        myStructure.ConstitutiveLawMultiPhysicsAddConstitutiveLaw(myConstitutiveLawMatrixMultiPhysics,myConstitutiveLawMatrixDryingShrinkage);
        myStructure.ConstitutiveLawMultiPhysicsAddConstitutiveLaw(myConstitutiveLawMatrixMultiPhysics,myConstitutiveLawMatrixMoistureTransport);

        myStructure.ConstitutiveLawSetParameterDouble(myConstitutiveLawMatrixLinearElastic, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, YoungsModulusMatrix);
        myStructure.ConstitutiveLawSetParameterDouble(myConstitutiveLawMatrixLinearElastic, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, PoissonRatioMatrix);
        myStructure.ConstitutiveLawSetParameterDouble(myConstitutiveLawMatrixLinearElastic, NuTo::Constitutive::eConstitutiveParameter::DENSITY,        DensityMatrix);

        myStructure.ConstitutiveLawSetParameterBool    (myConstitutiveLawMatrixMoistureTransport,"ENABLE_MODIFIED_TANGENTIAL_STIFFNESS",EnableModiefiedTangentialStiffness);      // sets whether modified tangential stiffness should be used or not
        myStructure.ConstitutiveLawSetParameterBool    (myConstitutiveLawMatrixMoistureTransport,"enable_sorption_hysteresis",EnableSorptionHysteresis);                          // sets whether sorption hysteresis should be used or not


        myStructure.ConstitutiveLawSetParameterDouble  (myConstitutiveLawMatrixMoistureTransport,NuTo::Constitutive::eConstitutiveParameter::BOUNDARY_TRANSPORT_CONSTANT_GAS_PHASE,BC_Surface_Moisture_Transfer_RH);        // set water phase density
        myStructure.ConstitutiveLawSetParameterDouble  (myConstitutiveLawMatrixMoistureTransport,NuTo::Constitutive::eConstitutiveParameter::BOUNDARY_TRANSPORT_CONSTANT_WATER_PHASE,BC_Surface_Moisture_Transfer_WVF);     // set water phase density
        myStructure.ConstitutiveLawSetParameterDouble  (myConstitutiveLawMatrixMoistureTransport,NuTo::Constitutive::eConstitutiveParameter::DENSITY_WATER_PHASE,WaterPhaseDensity);                                        // set water phase density
        myStructure.ConstitutiveLawSetParameterDouble  (myConstitutiveLawMatrixMoistureTransport,NuTo::Constitutive::eConstitutiveParameter::DIFFUSION_CONSTANT_GAS_PHASE,VaporPhaseDiffusionCoefficient);                  // set vapor phase diffusion coefficient
        myStructure.ConstitutiveLawSetParameterDouble  (myConstitutiveLawMatrixMoistureTransport,NuTo::Constitutive::eConstitutiveParameter::DIFFUSION_CONSTANT_WATER_PHASE,WaterPhaseDiffusionCoefficient);                // set water phase diffusion coefficient
        myStructure.ConstitutiveLawSetParameterDouble  (myConstitutiveLawMatrixMoistureTransport,NuTo::Constitutive::eConstitutiveParameter::DIFFUSION_EXPONENT_GAS_PHASE,VaporPhaseDiffusionExponent);                     // set vapor phase diffusion exponent
        myStructure.ConstitutiveLawSetParameterDouble  (myConstitutiveLawMatrixMoistureTransport,NuTo::Constitutive::eConstitutiveParameter::DIFFUSION_EXPONENT_WATER_PHASE,WaterPhaseDiffusionExponent);                   // set water phase diffusion exponent
        myStructure.ConstitutiveLawSetParameterDouble  (myConstitutiveLawMatrixMoistureTransport,NuTo::Constitutive::eConstitutiveParameter::GRADIENT_CORRECTION_ADSORPTION_DESORPTION,Kd);                                 // set gradient correction when changing from adsorption to desorption
        myStructure.ConstitutiveLawSetParameterDouble  (myConstitutiveLawMatrixMoistureTransport,NuTo::Constitutive::eConstitutiveParameter::GRADIENT_CORRECTION_DESORPTION_ADSORPTION,Ka);                                 // set gradient correction when changing from desorption to adsorption
        myStructure.ConstitutiveLawSetParameterDouble  (myConstitutiveLawMatrixMoistureTransport,NuTo::Constitutive::eConstitutiveParameter::MASS_EXCHANGE_RATE,MassExchangeRate);                                          // set mass exchange rate
        myStructure.ConstitutiveLawSetParameterDouble  (myConstitutiveLawMatrixMoistureTransport,NuTo::Constitutive::eConstitutiveParameter::POROSITY,Porosity);                                                            // set porosity
        myStructure.ConstitutiveLawSetParameterDouble  (myConstitutiveLawMatrixMoistureTransport,NuTo::Constitutive::eConstitutiveParameter::SATURATION_DENSITY_GAS_PHASE,VaporPhaseSaturationDensity);                     // set vapor phase saturation density

        myStructure.ConstitutiveLawSetParameterFullVectorDouble    (myConstitutiveLawMatrixMoistureTransport,NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_ADSORPTION,AdsorptionFit.GetPolynomialCoefficients());               // set adsorption coefficients
        myStructure.ConstitutiveLawSetParameterFullVectorDouble    (myConstitutiveLawMatrixMoistureTransport,NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION,DesorptionFit.GetPolynomialCoefficients());               // set desorption coefficients


        // Calculate equilibrium water volume fraction
        InitialWaterVolumeFraction   = myStructure.ConstitutiveLawGetEquilibriumWaterVolumeFraction(myConstitutiveLawMatrixMoistureTransport,
                                                                                                    InitialRelativeHumidity,myStructure.ConstitutiveLawGetParameterFullVectorDouble(myConstitutiveLawMatrixMoistureTransport,NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));

        myStructure.ElementGroupSetConstitutiveLaw(groupMatrix, myConstitutiveLawMatrixMultiPhysics);





        // Aggregates
        // ----------

        int myConstitutiveLawAggregates = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::LINEAR_ELASTIC_ENGINEERING_STRESS);
        myStructure.ConstitutiveLawSetParameterDouble(myConstitutiveLawAggregates, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, YoungsModulusAggregates);
        myStructure.ConstitutiveLawSetParameterDouble(myConstitutiveLawAggregates, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, PoissonRatioAggregates);

        myStructure.ElementGroupSetConstitutiveLaw(groupAggregates, myConstitutiveLawAggregates);




    /* -------------------------------- *\
    |     Create Boundary Elements       |
    \* -------------------------------- */


        // Add Nodes to boundary
        int nodeGroupBoundary = myStructure.GroupCreate("NODES");


        auto GetBoundaryNodesLambda = [lX,lY](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNumCoordinates()>0)
                                    {
                                        double x = rNodePtr->GetCoordinate(0);
                                        double y = rNodePtr->GetCoordinate(1);
                                        if ((x >=        -Tol   && x <=          Tol) ||
                                            (x >= lX     -Tol   && x <= lX     + Tol) ||
                                            (y >=        -Tol   && y <=          Tol) ||
                                            (y >= lY     -Tol   && y <= lY     + Tol) )
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };

        myStructure.GroupAddNodeFunction(nodeGroupBoundary,GetBoundaryNodesLambda);

        // Group all elements with boundary nodes
        int elemGroupBoundary = myStructure.GroupCreate("ELEMENTS");
        myStructure.GroupAddElementsFromNodes(elemGroupBoundary, nodeGroupBoundary, false);



        NuTo::FullVector<double,Eigen::Dynamic> tmpCoord(2);
        tmpCoord[0] = 12.0;
        tmpCoord[1] = 12.0;
        int BoundaryNodeID = myStructure.NodeCreateDOFs("RELATIVEHUMIDITY",tmpCoord);


        myStructure.BoundaryElementsCreate(elemGroupBoundary,nodeGroupBoundary,myStructure.NodeGetNodePtr(BoundaryNodeID));








    /* -------------------------------- *\
    |          Set Static Data           |
    \* -------------------------------- */


        // Loop over all integration points
        NuTo::FullVector<int,Eigen::Dynamic> GroupIdsMatrix = myStructure.GroupGetMemberIds(groupMatrix);
        for (int i=0; i<GroupIdsMatrix.GetNumRows(); i++)
        {
            NuTo::ElementBase* ElementPtr = myStructure.ElementGetElementPtr(GroupIdsMatrix[i]);
            for (int theIP=0; theIP< ElementPtr->GetNumIntegrationPoints(); theIP++)
            {
                NuTo::ConstitutiveStaticDataMultiPhysics *StaticDataMultiPhysics = ElementPtr->GetStaticData(theIP)->AsMultiPhysics();
                StaticDataMultiPhysics->AddNewStaticData(NuTo::Constitutive::eConstitutiveStaticDataType::MOISTURE_TRANSPORT);
                NuTo::ConstitutiveStaticDataMoistureTransport *StaticDataMoistureTransport = StaticDataMultiPhysics->AsMoistureTransport();

                StaticDataMoistureTransport->SetLastSorptionCoeff(myStructure.ConstitutiveLawGetParameterFullVectorDouble(myConstitutiveLawMatrixMoistureTransport,NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
                StaticDataMoistureTransport->SetActualSorptionCoeff(myStructure.ConstitutiveLawGetParameterFullVectorDouble(myConstitutiveLawMatrixMoistureTransport,NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
                StaticDataMoistureTransport->SetLastRelHumValue(InitialRelativeHumidity);
                StaticDataMoistureTransport->SetSorptionHistoryDesorption(SorptionHistoryDesorption);
            }

            for(unsigned int j = 0; j<ElementPtr->GetNumNodes(); ++j)
            {
                NuTo::NodeBase* NodePtr = ElementPtr->GetNode(j);
                if(NodePtr->GetNumRelativeHumidity() != 0)
                {
                    NodePtr->SetRelativeHumidity(0,InitialRelativeHumidity);
                }
                if(NodePtr->GetNumWaterVolumeFraction() != 0)
                {
                    NodePtr->SetWaterVolumeFraction(0,InitialWaterVolumeFraction);
                }

            }
        }

        // Add Nodes
        int nodeGroupAll = myStructure.GroupCreate("NODES");

        auto GetAllNodesLambda = [](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    return true;
                                };
        myStructure.GroupAddNodeFunction(nodeGroupAll,GetAllNodesLambda);

        int elemGroupAll = myStructure.GroupCreate("ELEMENTS");
        myStructure.GroupAddElementsFromNodes(elemGroupAll, nodeGroupAll, false);

        NuTo::FullVector<int,Eigen::Dynamic> GroupIdsAll = myStructure.GroupGetMemberIds(elemGroupAll);

        for (int i=0; i<GroupIdsAll.GetNumRows(); i++)
        {
            NuTo::ElementBase* ElementPtr = myStructure.ElementGetElementPtr(GroupIdsAll[i]);


            for (int theIP=0; theIP< ElementPtr->GetNumIntegrationPoints(); theIP++)
            {
                if(ElementPtr->GetConstitutiveLaw(theIP)->GetType()==NuTo::Constitutive::MULTI_PHYSICS)
                {
                    NuTo::ConstitutiveStaticDataMultiPhysics *StaticDataMultiPhysics = ElementPtr->GetStaticData(theIP)->AsMultiPhysics();
                    StaticDataMultiPhysics->AddNewStaticData(NuTo::Constitutive::eConstitutiveStaticDataType::MOISTURE_TRANSPORT);
                    NuTo::ConstitutiveStaticDataMoistureTransport *StaticDataMoistureTransport = StaticDataMultiPhysics->AsMoistureTransport();

                    StaticDataMoistureTransport->SetLastSorptionCoeff(myStructure.ConstitutiveLawGetParameterFullVectorDouble(myConstitutiveLawMatrixMoistureTransport,NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
                    StaticDataMoistureTransport->SetActualSorptionCoeff(myStructure.ConstitutiveLawGetParameterFullVectorDouble(myConstitutiveLawMatrixMoistureTransport,NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
                    StaticDataMoistureTransport->SetLastRelHumValue(InitialRelativeHumidity);
                    StaticDataMoistureTransport->SetSorptionHistoryDesorption(SorptionHistoryDesorption);
                }
            }
        }

    /* -------------------------------- *\
    |           Create Section           |
    \* -------------------------------- */

        int mySection = myStructure.SectionCreate(NuTo::Section::PLANE_STRESS);
        myStructure.SectionSetThickness(mySection, lZ);
        myStructure.ElementTotalSetSection(mySection);





    /* -------------------------------- *\
    |          Set Constraints           |
    \* -------------------------------- */

        auto RelativeHumidityFunc = [InitialRelativeHumidity]
                                     (double rTime)->double
                                  {
                                        double StartTime = 60.0*60.0*24;
                                        double EndTime = 60*60*24*2;

                                        if(rTime < StartTime)
                                        {
                                            return InitialRelativeHumidity;
                                        }
                                        else
                                        {
                                            if(rTime< EndTime)
                                            {
//                                                double test = InitialRelativeHumidity - sin((rTime - StartTime) / (EndTime-StartTime) * 3.14 /2.0) * (InitialRelativeHumidity-0.45);
                                                return InitialRelativeHumidity - sin((rTime - StartTime) / (EndTime-StartTime) * 3.14 /2.0) * (InitialRelativeHumidity-0.45);
                                            }
                                            {
                                                return 0.45;
                                            }
                                        }
                                  };

        NuTo::FullVector<int,Eigen::Dynamic> NodeGroupIdsAll = myStructure.GroupGetMemberIds(nodeGroupAll);


        int IDNodeLeft   = GetNearestNode(0.0    ,lY/2.0, myStructure,NodeGroupIdsAll);
        int IDNodeCenter = GetNearestNode(lX/2.0 ,lY/2.0, myStructure,NodeGroupIdsAll);



        int groupNodeLeft   = myStructure.GroupCreate(NuTo::Groups::Nodes);
        int groupNodeCenter = myStructure.GroupCreate(NuTo::Groups::Nodes);


        myStructure.GroupAddNode(groupNodeLeft  ,IDNodeLeft);
        myStructure.GroupAddNode(groupNodeCenter,IDNodeCenter);

        myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeLeft,   NuTo::FullVector<double,2>::UnitY(), 0.);
        myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeCenter, NuTo::FullVector<double,2>::UnitX(), 0.);
        myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeCenter, NuTo::FullVector<double,2>::UnitY(), 0.);

        double deltaD = 0.001;

        //int gNodesWest = myStructure.GroupCreate(NuTo::Groups::Nodes);
        //int gNodesEast = myStructure.GroupCreate(NuTo::Groups::Nodes);

        //int iNodeOrigin = myStructure.NodeGetIdAtCoordinate(NuTo::FullVector<double, 2>({0.,0.}), 1.e-6);
        //myStructure.GroupAddNodeCoordinateRange(gNodesWest, 0, 0.,0.);
        //myStructure.GroupAddNodeCoordinateRange(gNodesEast, 0, lX,lX);


        //myStructure.ConstraintLinearSetDisplacementNodeGroup(gNodesWest, NuTo::FullVector<double,2>::UnitX(), 0.);
        //myStructure.ConstraintLinearSetDisplacementNode(iNodeOrigin, NuTo::FullVector<double,2>::UnitX(), 0.);
        //myStructure.ConstraintLinearSetDisplacementNode(iNodeOrigin, NuTo::FullVector<double,2>::UnitY(), 0.);

        //int bc = myStructure.ConstraintLinearSetDisplacementNodeGroup(gNodesEast, NuTo::FullVector<double,2>::UnitY(), 0);

        int BoundaryConstraint = myStructure.ConstraintLinearSetRelativeHumidityNode(BoundaryNodeID,InitialRelativeHumidity);

    /* -------------------------------- *\
    |           Visulaization            |
    \* -------------------------------- */

#ifdef ENABLE_VISUALIZE
        myStructure.AddVisualizationComponent(groupMatrix,NuTo::VisualizeBase::RELATIVE_HUMIDITY);
        myStructure.AddVisualizationComponent(groupMatrix,NuTo::VisualizeBase::WATER_VOLUME_FRACTION);
        myStructure.AddVisualizationComponent(groupMatrix,NuTo::VisualizeBase::DISPLACEMENTS);
        myStructure.AddVisualizationComponent(groupMatrix,NuTo::VisualizeBase::ENGINEERING_STRAIN);
        myStructure.AddVisualizationComponent(groupMatrix,NuTo::VisualizeBase::SECTION);
        myStructure.AddVisualizationComponent(groupMatrix,NuTo::VisualizeBase::ENGINEERING_STRESS);
        myStructure.AddVisualizationComponent(groupMatrix,NuTo::VisualizeBase::PRINCIPAL_ENGINEERING_STRESS);


        myStructure.ExportVtkDataFileElements("test.vtk");
#endif // ENABLE_VISUALIZE


        //**********************************************
        //          Solver
        //**********************************************

        myStructure.NodeBuildGlobalDofs();
        myStructure.CalculateMaximumIndependentSets();
//        NuTo::CrankNicolsonEvaluate myTimeIntegrationScheme(&myStructure);
        NuTo::NewmarkDirect myTimeIntegrationScheme(&myStructure);
        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> dispRHS(2,2);
        dispRHS << 0, 0, SimulationTime, deltaD;

        myTimeIntegrationScheme.AddTimeDependentConstraintFunction(BoundaryConstraint,  RelativeHumidityFunc);

        //myTimeIntegrationScheme.AddTimeDependentConstraint(bc, dispRHS);
        //myTimeIntegrationScheme.SetCheckEquilibriumOnStart(false);
        myTimeIntegrationScheme.SetMaxTimeStep(delta_t);
        myTimeIntegrationScheme.SetMinTimeStep(0.000001);
        myTimeIntegrationScheme.SetTimeStep(delta_t);
        myTimeIntegrationScheme.SetPerformLineSearch(false);
//        myTimeIntegrationScheme.SetToleranceForce(1e-7);
        myTimeIntegrationScheme.SetAutomaticTimeStepping(false);
        myTimeIntegrationScheme.SetVerboseLevel(0);
        myTimeIntegrationScheme.SetShowTime(true);
        myTimeIntegrationScheme.SetMinTimeStepPlot(delta_t_write);
        myTimeIntegrationScheme.SetNewmarkBeta(0.25);
        myTimeIntegrationScheme.SetNewmarkGamma(0.5);
        myTimeIntegrationScheme.SetToleranceResidual(NuTo::Node::eDof::DISPLACEMENTS,1e-7);
        myTimeIntegrationScheme.SetToleranceResidual(NuTo::Node::eDof::WATERVOLUMEFRACTION,1e-7);
        myTimeIntegrationScheme.SetToleranceResidual(NuTo::Node::eDof::RELATIVEHUMIDITY,1e-7);
//        myTimeIntegrationScheme.AddCalculationStep({NuTo::Node::eDof::WATERVOLUMEFRACTION,NuTo::Node::eDof::RELATIVEHUMIDITY});
//        myTimeIntegrationScheme.AddCalculationStep({NuTo::Node::eDof::DISPLACEMENTS});
        //myTimeIntegrationScheme.SetPlotElementGroups();

        bool deleteDirectory = true;
        myTimeIntegrationScheme.SetResultDirectory(resultPath.string(), deleteDirectory);

        myTimeIntegrationScheme.SetMaxNumIterations(40);
//        myTimeIntegrationScheme.Solve(SimulationTime);
        myTimeIntegrationScheme.Solve(SimulationTime);
    } catch (NuTo::Exception& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
        std::cout << "\n\n\n Errors occured! \n\n\n" << std::endl;
    }

    std::cout << "Calculation finished" << std::endl;

    return 0;


}
