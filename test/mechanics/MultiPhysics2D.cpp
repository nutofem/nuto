
#include <fstream>
#include <sys/stat.h>
#include <sys/time.h>

#if defined HAVE_PARDISO
    #include <nuto/math/SparseDirectSolverPardiso.h>
#elif defined HAVE_MUMPS
    #include <nuto/math/SparseDirectSolverMUMPS.h>
#else
    std::cout << "Solver not available - can't solve system of equations " << std::endl;
#endif

#include <boost/progress.hpp>

#include <nuto/math/SparseMatrixCSRGeneral.h>


#include <nuto/mechanics/constitutive/multiPhysics/ConstitutiveStaticDataMultiPhysics.h>
#include <nuto/mechanics/constitutive/moistureTransport/ConstitutiveStaticDataMoistureTransport.h>
#include <nuto/mechanics/nodes/NodeDof.h>
#include <nuto/mechanics/nodes/NodeCoordinates.h>
#include <nuto/mechanics/structures/unstructured/Structure.h>
#include "nuto/mechanics/timeIntegration/CrankNicolsonEvaluate.h"

#include <nuto/metamodel/PolynomialLeastSquaresFitting.h>

#include <nuto/mechanics/structures/StructureOutputFullVectorDouble.h>
#include <nuto/mechanics/structures/StructureOutputSparseMatrix.h>


int main()
{
    try
    {



        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // Declaration of neccessary variables
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        bool            EnableSorptionHysteresis            = false;
        bool            EnableModiefiedTangentialStiffness  = false;

        //unsigned int    NElements   = 16;                             // Number of elements
        unsigned int    NNodes;                                         // Number of nodes

        double Length           = 0.16;
        double Height           = 0.04;

        unsigned int NumElementsX =  32;
        unsigned int NumElementsY =   8;

        //double ElementHeight    = Height / static_cast<double>(NumElementsY);
        double ElementLength    = Length / static_cast<double>(NumElementsX);


        unsigned int NumNodesX = NumElementsX + 1;
        unsigned int NumNodesY = NumElementsY + 1;

        //double ElementHeight    = static_cast<double>(Height)/static_cast<double>(NumNodesY-1);
        //double ElementLength    = static_cast<double>(Length)/static_cast<double>(NumNodesX-1);


        double          Thickness    = 0.04;

        double          delta_t                         = 1.0/1.0 *     1.0 * 24.0 * 60.0 * 60.0;
        double          SimulationTime                  = 1.0/1.0 *     1.0 * 24.0 * 60.0 * 60.0;
        double          BC_TransitionTime               =                     24.0 * 60.0 * 60.0;

        // initial node values
        double          InitialRelativeHumidity         =    0.95;
        double          InitialWaterVolumeFraction      =    0.03;


        // constitutive law values
        // -----------------------

        // Linear Elastic
        double Density                                  = 1.0;      // N/mm³
        double PoissonRatio                             = 0.15;
        double YoungsModulus                            = 30.0 * 10e9;  //N/m²
        // Moisture Transport
        double          MassExchangeRate                =    3.42e-7    ;
        double          Porosity                        =    0.25      ;
        double          VaporPhaseDiffusionCoefficient  =    3.9e-12     ;
        double          VaporPhaseDiffusionExponent     =    1.0        ;
        double          VaporPhaseSaturationDensity     =    0.0173     ;
        double          WaterPhaseDensity               =  999.97       ;
        double          WaterPhaseDiffusionCoefficient  =    1.17e-7    ;
        double          WaterPhaseDiffusionExponent     =    2.0        ;

        // Boundary Condition Values
        double          BC_RelativeHumidity             =    0.15;
        double          BC_Surface_Moisture_Transfer_RH =    1.0e-10 * delta_t;
        double          BC_Surface_Moisture_Transfer_WVF=    1.0e-7 * delta_t;
        bool            SorptionHistoryDesorption       =    true;

        // sorption hysteresis
        double          Ka                              =    0.26       ;
        double          Kd                              =    0.56       ;
        NuTo::PolynomialLeastSquaresFitting AdsorptionFit;
        NuTo::PolynomialLeastSquaresFitting DesorptionFit;


        // max residual
        double          MaxResidual                     =    1.0e-6;





        // %%%%%%%%%%%%%%%%%%%
        // Fit Sorption Curves
        // %%%%%%%%%%%%%%%%%%%


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



        // %%%%%%%%%%%%%%%%
        // Create Structure
        // %%%%%%%%%%%%%%%%


        NuTo::Structure myStructure(2);
        myStructure.SetNumTimeDerivatives(2);

        // disable output of calculation times
        myStructure.SetShowTime(false);
        myStructure.UseMaximumIndependentSets(true);





        // %%%%%%%%%%%%%%
        // Create Section
        // %%%%%%%%%%%%%%


        int mySection = myStructure.SectionCreate("Plane_Stress");
        myStructure.SectionSetThickness(mySection,Thickness);
        myStructure.ElementTotalSetSection(mySection);



        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // Create and set Constitutive Law
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        int ConstLawDryingShrinkage   = myStructure.ConstitutiveLawCreate("DryingShrinkage");
        int ConstLawLinearElastic     = myStructure.ConstitutiveLawCreate("LinearElasticEngineeringStress");
        int ConstLawMoistureTransport = myStructure.ConstitutiveLawCreate("MoistureTransport");
        int ConstLawMultiPhysics      = myStructure.ConstitutiveLawCreate("MultiPhysics");
        myStructure.ConstitutiveLawMultiPhysicsAddConstitutiveLaw(ConstLawMultiPhysics,ConstLawDryingShrinkage);
        myStructure.ConstitutiveLawMultiPhysicsAddConstitutiveLaw(ConstLawMultiPhysics,ConstLawLinearElastic);
        myStructure.ConstitutiveLawMultiPhysicsAddConstitutiveLaw(ConstLawMultiPhysics,ConstLawMoistureTransport);


        // set variables

        myStructure.ConstitutiveLawSetParameterDouble  (ConstLawMultiPhysics,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,YoungsModulus);
        myStructure.ConstitutiveLawSetParameterDouble  (ConstLawMultiPhysics,NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,PoissonRatio);
        myStructure.ConstitutiveLawSetParameterDouble  (ConstLawMultiPhysics,NuTo::Constitutive::eConstitutiveParameter::DENSITY,Density);

        myStructure.ConstitutiveLawSetParameterBool    (ConstLawMultiPhysics,"ENABLE_MODIFIED_TANGENTIAL_STIFFNESS",EnableModiefiedTangentialStiffness);      // sets whether modified tangential stiffness should be used or not
        myStructure.ConstitutiveLawSetParameterBool    (ConstLawMultiPhysics,"enable_sorption_hysteresis",EnableSorptionHysteresis);                          // sets whether sorption hysteresis should be used or not


        myStructure.ConstitutiveLawSetParameterDouble  (ConstLawMultiPhysics,"boundary_TRANSPORT_CONSTANT_GAS_PHASE",BC_Surface_Moisture_Transfer_RH);        // set water phase density
        myStructure.ConstitutiveLawSetParameterDouble  (ConstLawMultiPhysics,"BOUNDARY_TRANSPORT_CONSTANT_WATER_PHASE",BC_Surface_Moisture_Transfer_WVF);     // set water phase density
        myStructure.ConstitutiveLawSetParameterDouble  (ConstLawMultiPhysics,"DENSITY_WATER_PHASE",WaterPhaseDensity);                                        // set water phase density
        myStructure.ConstitutiveLawSetParameterDouble  (ConstLawMultiPhysics,"DIFFUSION_CONSTANT_GAS_PHASE",VaporPhaseDiffusionCoefficient);                  // set vapor phase diffusion coefficient
        myStructure.ConstitutiveLawSetParameterDouble  (ConstLawMultiPhysics,"DIFFUSION_CONSTANT_WATER_PHASE",WaterPhaseDiffusionCoefficient);                // set water phase diffusion coefficient
        myStructure.ConstitutiveLawSetParameterDouble  (ConstLawMultiPhysics,"DIFFUSION_EXPONENT_GAS_PHASE",VaporPhaseDiffusionExponent);                     // set vapor phase diffusion exponent
        myStructure.ConstitutiveLawSetParameterDouble  (ConstLawMultiPhysics,"DIFFUSION_EXPONENT_WATER_PHASE",WaterPhaseDiffusionExponent);                   // set water phase diffusion exponent
        myStructure.ConstitutiveLawSetParameterDouble  (ConstLawMultiPhysics,"GRADIENT_CORRECTION_ADSORPTION_DESORPTION",Kd);                                 // set gradient correction when changing from adsorption to desorption
        myStructure.ConstitutiveLawSetParameterDouble  (ConstLawMultiPhysics,"GRADIENT_CORRECTION_DESORPTION_ADSORPTION",Ka);                                 // set gradient correction when changing from desorption to adsorption
        myStructure.ConstitutiveLawSetParameterDouble  (ConstLawMultiPhysics,"MASS_EXCHANGE_RATE",MassExchangeRate);                                          // set mass exchange rate
        myStructure.ConstitutiveLawSetParameterDouble  (ConstLawMultiPhysics,"POROSITY",Porosity);                                                            // set porosity
        myStructure.ConstitutiveLawSetParameterDouble  (ConstLawMultiPhysics,"SATURATION_DENSITY_GAS_PHASE",VaporPhaseSaturationDensity);                     // set vapor phase saturation density

        myStructure.ConstitutiveLawSetParameterFullVectorDouble    (ConstLawMultiPhysics,"polynomial_COEFFICIENTS_ADSORPTION",AdsorptionFit.GetPolynomialCoefficients());               // set adsorption coefficients
        myStructure.ConstitutiveLawSetParameterFullVectorDouble    (ConstLawMultiPhysics,"POLYNOMIAL_COEFFICIENTS_DESORPTION",DesorptionFit.GetPolynomialCoefficients());               // set desorption coefficients


        // Calculate equilibrium water volume fraction
        InitialWaterVolumeFraction   = myStructure.ConstitutiveLawGetEquilibriumWaterVolumeFraction(ConstLawMoistureTransport,InitialRelativeHumidity,myStructure.ConstitutiveLawGetParameterFullVectorDouble(ConstLawMoistureTransport,NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));



        // %%%%%%%%%%%%
        // Create Nodes
        // %%%%%%%%%%%%




        double DeltaX = Length/(NumNodesX-1);
        double DeltaY = Height/(NumNodesY-1);

        int NodeNum(0);
        NuTo::FullVector<double,Eigen::Dynamic> Coordinates(2);

        for (unsigned int iY=0; iY < NumNodesY; iY++)
        {
            for (unsigned int iX=0; iX < NumNodesX; iX++)
            {
                Coordinates(0) = iX * DeltaX;
                Coordinates(1) = iY * DeltaY;
                myStructure.NodeCreate(NodeNum,Coordinates);
                NodeNum++;
            }
        }



        // %%%%%%%%%%%%%%%%%%
        // Interpolation Type
        // %%%%%%%%%%%%%%%%%%

        int myInterpolationType = myStructure.InterpolationTypeCreate("Quad2D");
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::RELATIVEHUMIDITY, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::WATERVOLUMEFRACTION, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);






        // %%%%%%%%%%%%%%%
        // Create Elements
        // %%%%%%%%%%%%%%%


        NuTo::FullVector<int,Eigen::Dynamic> Nodes(4);

        for (unsigned int iY=0; iY< NumElementsY; iY++)
        {
            for (unsigned int iX=0; iX < NumElementsX; iX++)
            {
                Nodes(0) = iX    +  iY    *NumNodesX;
                Nodes(1) = iX +1 +  iY    *NumNodesX;
                Nodes(2) = iX +1 + (iY+1) *NumNodesX;
                Nodes(3) = iX    + (iY+1) *NumNodesX;
                int elementNumber = myStructure.ElementCreate(myInterpolationType, Nodes, "CONSTITUTIVELAWIP","STATICDATA");

                // set element section
                myStructure.ElementSetSection(elementNumber,mySection);

                // set element constitutive law
                myStructure.ElementSetConstitutiveLaw(elementNumber,ConstLawMultiPhysics);
            }
        }

        myStructure.ElementTotalConvertToInterpolationType();

//        // Apply nodal start values

        NNodes = myStructure.GetNumNodes();

        for (unsigned int i=0; i<NNodes; i++)
        {
            /*
            double NodeX = myStructure.NodeGetNodePtr(i)->GetCoordinates()[0] / Length;
            double NodeY = myStructure.NodeGetNodePtr(i)->GetCoordinates()[1] / Height;
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



        // %%%%%%%%%%
        // Set groups
        // %%%%%%%%%%

        //left boundary
        int GRPNodes_Left = myStructure.GroupCreate("Nodes");
        int Direction = 0;
        double Min = 0. - 0.01 * ElementLength;
        double Max = 0. + 0.01 * ElementLength;;
        myStructure.GroupAddNodeCoordinateRange(GRPNodes_Left,Direction,Min,Max);

        //right boundary
        int GRPNodes_Right = myStructure.GroupCreate("Nodes");
        Direction = 0;
        Min = Length - 0.01 * ElementLength;
        Max = Length + 0.01 * ElementLength;
        myStructure.GroupAddNodeCoordinateRange(GRPNodes_Right,Direction,Min,Max);



        auto PickNodeFunction = [Length,Height](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    if (rNodePtr->GetNumCoordinates()>0)
                                    {
                                        double x = rNodePtr->GetCoordinate(0);
                                        double y = rNodePtr->GetCoordinate(1);
                                        if (x >= Length/2 -0.1 &&
                                            x <= Length/2 +0.1 &&
                                            y >= Height / 2.0)
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };

        int GRPNodes_Center = myStructure.GroupCreate("Nodes");
        myStructure.GroupAddNodeFunction(GRPNodes_Center,PickNodeFunction);



        auto PickSingleNodeLeftFunction = [Length,Height](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    if (rNodePtr->GetNumCoordinates()>0)
                                    {
                                        double x = rNodePtr->GetCoordinate(0);
                                        double y = rNodePtr->GetCoordinate(1);
                                        if (x >= -1e-6 &&
                                            x <= 1e-6 &&
                                            y >= Height/2.0-1e-6 &&
                                            y <= Height/2.0+1e-6)
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };
        int GRPNodes_SingleNodeLeft = myStructure.GroupCreate("Nodes");
        myStructure.GroupAddNodeFunction(GRPNodes_SingleNodeLeft,PickSingleNodeLeftFunction);

        auto PickSingleNodeRightFunction = [Length,Height](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    if (rNodePtr->GetNumCoordinates()>0)
                                    {
                                        double x = rNodePtr->GetCoordinate(0);
                                        double y = rNodePtr->GetCoordinate(1);
                                        if (x >= Length-1e-6 &&
                                            x <= Length+1e-6 &&
                                            y >= Height/2.0-1e-6 &&
                                            y <= Height/2.0+1e-6)
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };
        int GRPNodes_SingleNodeRight = myStructure.GroupCreate("Nodes");
        myStructure.GroupAddNodeFunction(GRPNodes_SingleNodeRight,PickSingleNodeRightFunction);

        // %%%%%%%%%%%%%%%%%%%%%%%%
        // Create Boundary Elements
        // %%%%%%%%%%%%%%%%%%%%%%%%


        // Add Nodes to boundary
        int nodeGroupBoundary = myStructure.GroupCreate("NODES");


        auto GetBoundaryNodesLambda = [Length,Height](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNumCoordinates()>0)
                                    {
                                        double x = rNodePtr->GetCoordinate(0);
                                        double y = rNodePtr->GetCoordinate(1);
                                        if ((x >=        -Tol   && x <=          Tol) ||
                                            (x >= Length -Tol   && x <= Length + Tol) ||
                                            (y >=        -Tol   && y <=          Tol) ||
                                            (y >= Height -Tol   && y <= Height + Tol) )
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



        int BoundaryNodeID = myStructure.NodeCreateDOFs("RELATIVEHUMIDITY");


        myStructure.BoundaryElementsCreate(elemGroupBoundary,nodeGroupBoundary,myStructure.NodeGetNodePtr(BoundaryNodeID));



        // %%%%%%%%%%%%%%%
        // Set Static Data
        // %%%%%%%%%%%%%%%


        // Loop over all integration points
        for (int i=0; i<myStructure.GetNumElements(); i++)
        {
            for (int theIP=0; theIP< myStructure.ElementGetElementPtr(i)->GetNumIntegrationPoints(); theIP++)
            {
                NuTo::ConstitutiveStaticDataMultiPhysics *StaticDataMultiPhysics = myStructure.ElementGetElementPtr(i)->GetStaticData(theIP)->AsMultiPhysics();
                StaticDataMultiPhysics->AddNewStaticData(NuTo::Constitutive::eConstitutiveStaticDataType::MOISTURE_TRANSPORT);
                NuTo::ConstitutiveStaticDataMoistureTransport *StaticDataMoistureTransport = StaticDataMultiPhysics->AsMoistureTransport();

                StaticDataMoistureTransport->SetLastSorptionCoeff(myStructure.ConstitutiveLawGetParameterFullVectorDouble(ConstLawMoistureTransport,NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
                StaticDataMoistureTransport->SetActualSorptionCoeff(myStructure.ConstitutiveLawGetParameterFullVectorDouble(ConstLawMoistureTransport,NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
                StaticDataMoistureTransport->SetLastRelHumValue(InitialRelativeHumidity);
                StaticDataMoistureTransport->SetSorptionHistoryDesorption(SorptionHistoryDesorption);
            }
        }

        // %%%%%%%%%%%%%%%
        // Set Constraints
        // %%%%%%%%%%%%%%%

        auto RelativeHumidityFunc = [BC_RelativeHumidity,
                                     BC_TransitionTime,
                                     InitialRelativeHumidity,
                                     SimulationTime]
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
        //myStructure.ConstraintLinearSetWaterVolumeFractionNode(BoundaryNodeID,1.0);

        /*auto DispFunctionCenter = [SimulationTime](double rTime)->double
                                  {
                                        double QuarterTime = SimulationTime/4.0;
                                        double factor = 0.0;

                                        if (rTime<=1*QuarterTime)
                                        {
                                            factor = rTime/QuarterTime * 0.1;
                                        }
                                        else if(rTime<=3*QuarterTime)
                                        {
                                            factor = 0.1;
                                        }
                                        else
                                        {
                                            factor = 0.1 - (rTime-3*QuarterTime)/QuarterTime * 0.1;
                                        }

                                        return factor * sin(rTime)*0.2;
                                  };*/

        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> DirectionX(2,1);
        DirectionX.SetValue(0,0,1.0);

        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> DirectionY(2,1);
        DirectionY.SetValue(1,0,1.0);


        myStructure.ConstraintLinearSetDisplacementNodeGroup(GRPNodes_SingleNodeLeft           ,DirectionX,0);
        myStructure.ConstraintLinearSetDisplacementNodeGroup(GRPNodes_SingleNodeLeft           ,DirectionY,0);
        myStructure.ConstraintLinearSetDisplacementNodeGroup(GRPNodes_SingleNodeRight          ,DirectionY,0);
        //myStructure.ConstraintLinearSetDisplacementNodeGroup(GRPNodes_Right          ,DirectionY,0);
        //int ConstraintRightDisp = myStructure.ConstraintLinearSetDisplacementNodeGroup(GRPNodes_Right          ,DirectionX,0);
        //int ConstraintCenterDisp = myStructure.ConstraintLinearSetDisplacementNodeGroup(GRPNodes_Center         ,DirectionY,0);

        // %%%%%%%%%%%%%%%%%%%%%
        // Multi processor Setup
        // %%%%%%%%%%%%%%%%%%%%%


#ifdef _OPENMP
        myStructure.SetNumProcessors(4);
        std::cout << "OpenMP enabled" << std::endl;
#else
        myStructure.SetNumProcessors(1);
#endif
        myStructure.CalculateMaximumIndependentSets();
        myStructure.UseMaximumIndependentSets(true);
        myStructure.NodeBuildGlobalDofs();



        //myStructure.NodeInfo(10);

        // %%%%%%%%%%
        // Set Solver
        // %%%%%%%%%%


#if defined HAVE_PARDISO
        NuTo::SparseDirectSolverPardiso Solver(4);
#elif defined HAVE_MUMPS
        NuTo::SparseDirectSolverMUMPS Solver;
#else
        std::cout << "Solver not available - can't solve system of equations " << std::endl;
#endif
        Solver.SetShowTime(false);


        // %%%%%%%%%%%%%%%%
        // Time Integration
        // %%%%%%%%%%%%%%%%

#ifdef ENABLE_VISUALIZE

        int visualizationGroup = myStructure.GroupCreate(NuTo::Groups::eGroupId::Elements);
        myStructure.GroupAddElementsTotal(visualizationGroup);

        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::RELATIVE_HUMIDITY);
        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::WATER_VOLUME_FRACTION);
        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::DISPLACEMENTS);
        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRAIN);
        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRESS);

#endif // ENABLE_VISUALIZE


        NuTo::CrankNicolsonEvaluate myTimeIntegrationScheme(&myStructure);


        //myTimeIntegrationScheme.SetMaxTimeStep(delta_t);
        //myTimeIntegrationScheme.SetMinTimeStep(delta_t);

        myTimeIntegrationScheme.SetPerformLineSearch(false);
        myTimeIntegrationScheme.SetCheckEquilibriumOnStart(false);
        myTimeIntegrationScheme.SetToleranceForce(MaxResidual);
        myTimeIntegrationScheme.SetMaxNumIterations(40);

        myTimeIntegrationScheme.AddTimeDependentConstraintFunction(BoundaryConstraint,  RelativeHumidityFunc);
        //myTimeIntegrationScheme.AddTimeDependentConstraint(ConstraintRightDisp, DispFunctionCenter);

        //set result directory
        bool deleteResultDirectoryFirst(true);
        myTimeIntegrationScheme.SetResultDirectory("./ResultsMultiPhysics2D",deleteResultDirectoryFirst);


        myTimeIntegrationScheme.SetTimeStep(delta_t);
        myTimeIntegrationScheme.SetMinTimeStepPlot(delta_t);
        myTimeIntegrationScheme.Solve(SimulationTime);


    }
    catch(NuTo::Exception e)
    {
        std::cout << e.ErrorMessage() << std::endl;
        return 1;
    }
    catch(std::exception e)
    {
            std::cout << e.what() << std::endl;
            return 1;
    }
}


// paraview calc Code:
// coordsX * iHat + coordsY * jHat + (coordsZ + WaterVolumeFraction) * kHat
