#include <nuto/math/SparseDirectSolverMUMPS.h>
#include <nuto/math/SparseMatrixCSRGeneral.h>

#include <nuto/mechanics/constitutive/moistureTransport/ConstitutiveStaticDataMoistureTransport.h>
#include <nuto/mechanics/nodes/NodeDof.h>
#include <nuto/mechanics/nodes/NodeCoordinates.h>
#include <nuto/mechanics/nodes/NodeCoordinatesDof.h>
#include <nuto/mechanics/structures/unstructured/Structure.h>

#include <nuto/metamodel/PolynomialLeastSquaresFitting.h>

int main()
{
    try
    {



        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // Declaration of neccessary variables
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        unsigned int    NElements   = 50;                               // Number of elements
        unsigned int    NNodes      = NElements+1;                      // Number of nodes
        double          L           = 0.16;                             // Length of the specimen
        double          Area        = 0.04*0.04;

        double          delta_t_fix = 1.0/3.0  * 24 * 3600.0;
        double          t           = 0.0;
        double          t_final     = 1.0/1.0 * 1 * 40 * 24 *3600.0;
        double          t_write     = 1.0/1.0 * 1 * 24 * 3600.0;

        // initial node values
        double          InitialRelativeHumidity         =    0.95;
        double          InitialWaterPhaseFraction       =    0.03;

        // Boundary Condition Values
        double          BC_RelativeHumidity             =    0.7;
        double          BC_WaterVolumeFraction;
        double          BC_Surface_Moisture_Transfer    =    1.0e-4;
        bool            SorptionHistoryDesorption       =    true;

        // constitutive law values
        double          MassExchangeRate                =    3.42e-7    ;
        double          Porosity                        =    0.25      ;
        double          VaporPhaseDiffusionCoefficient  =    3.9e-10     ;
        double          VaporPhaseDiffusionExponent     =    1.0        ;
        double          VaporPhaseSaturationDensity     =    0.0173     ;
        double          WaterPhaseDensity               =  999.97       ;
        double          WaterPhaseDiffusionCoefficient  =    1.17e-7    ;
        double          WaterPhaseDiffusionExponent     =    2.0        ;

        // sorption hysteresis
        double          Ka                              =    0.26       ;
        double          Kd                              =    0.56       ;
        NuTo::PolynomialLeastSquaresFitting AdsorptionFit;
        NuTo::PolynomialLeastSquaresFitting DesorptionFit;


        // max residual
        double          MaxResidual                     =    1.0e-6 * Area;



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


        NuTo::Structure MTStructure1D(1);



        // %%%%%%%%%%%%%%
        // Create Section
        // %%%%%%%%%%%%%%


        int Section1 = MTStructure1D.SectionCreate("Truss");
        MTStructure1D.SectionSetArea(Section1, Area);
        MTStructure1D.SectionSetDOF(Section1,"WaterPhaseFraction RelativeHumidity");



        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // Create and set Constitutive Law
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        int ConstLaw = MTStructure1D.ConstitutiveLawCreate("MoistureTransport");

        // set variables
        MTStructure1D.ConstitutiveLawSetAdsorptionCoefficients                      (ConstLaw,AdsorptionFit.GetPolynomialCoefficients());               // set adsorption coefficients
        MTStructure1D.ConstitutiveLawSetBoundarySurfaceMoistureTransportCoefficient (ConstLaw,BC_Surface_Moisture_Transfer);                            // set boundary surface moisture transport coefficient
        MTStructure1D.ConstitutiveLawSetDesorptionCoefficients                      (ConstLaw,DesorptionFit.GetPolynomialCoefficients());               // set desorption coefficients
        MTStructure1D.ConstitutiveLawSetKa                                          (ConstLaw,Ka);                                                      // set Ka
        MTStructure1D.ConstitutiveLawSetKd                                          (ConstLaw,Kd);                                                      // set Kd
        MTStructure1D.ConstitutiveLawSetMassExchangeRate                            (ConstLaw,MassExchangeRate);                                        // set mass exchange rate
        MTStructure1D.ConstitutiveLawSetPorosity                                    (ConstLaw,Porosity);                                                // set porosity
        MTStructure1D.ConstitutiveLawSetVaporPhaseDiffusionCoefficient              (ConstLaw,VaporPhaseDiffusionCoefficient);                          // set vapor phase diffusion coefficient
        MTStructure1D.ConstitutiveLawSetVaporPhaseDiffusionExponent                 (ConstLaw,VaporPhaseDiffusionExponent);                             // set vapor phase diffusion exponent
        MTStructure1D.ConstitutiveLawSetVaporPhaseSaturationDensity                 (ConstLaw,VaporPhaseSaturationDensity);                             // set vapor phase saturation density
        MTStructure1D.ConstitutiveLawSetWaterPhaseDensity                           (ConstLaw,WaterPhaseDensity);                                       // set water phase density
        MTStructure1D.ConstitutiveLawSetWaterPhaseDiffusionCoefficient              (ConstLaw,WaterPhaseDiffusionCoefficient);                          // set water phase diffusion coefficient
        MTStructure1D.ConstitutiveLawSetWaterPhaseDiffusionExponent                 (ConstLaw,WaterPhaseDiffusionExponent);                             // set water phase diffusion exponent

        // Calculate equilibrium water volume fraction
        InitialWaterPhaseFraction   = MTStructure1D.ConstitutiveLawGetEquilibriumWaterVolumeFraction(ConstLaw,InitialRelativeHumidity,MTStructure1D.ConstitutiveLawGetDesorptionCoefficients(ConstLaw));
        BC_WaterVolumeFraction      = MTStructure1D.ConstitutiveLawGetEquilibriumWaterVolumeFraction(ConstLaw,BC_RelativeHumidity,MTStructure1D.ConstitutiveLawGetDesorptionCoefficients(ConstLaw));



        // %%%%%%%%%%%%
        // Create Nodes
        // %%%%%%%%%%%%


        // Calculating Coordinates
        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> NodeCoordinates(1,NNodes);

        for (unsigned int i=0; i<NNodes; i++)
        {
            NodeCoordinates(0,i)= i*(L/(NNodes-1));
        }

        // Create Nodes
        MTStructure1D.NodesCreate("waterphasefraction relativehumidity", NodeCoordinates);

        for (unsigned int i=0; i<NNodes; i++)
        {
            MTStructure1D.NodeGetNodePtr(i)->SetRelativeHumidity(0,InitialRelativeHumidity);
            MTStructure1D.NodeGetNodePtr(i)->SetWaterPhaseFraction(0,InitialWaterPhaseFraction);
        }



        // %%%%%%%%%%%%%%%
        // Create Elements
        // %%%%%%%%%%%%%%%


        NuTo::FullVector<int, Eigen::Dynamic> ElementNodeNumbers(2);

        for (unsigned int i=0; i<NElements; i++)
        {
            // set node numbers
            ElementNodeNumbers(0) = i;
            ElementNodeNumbers(1) = i+1;

            // create elements
            MTStructure1D.ElementCreate(i,"Truss1D2N",ElementNodeNumbers,"CONSTITUTIVELAWIP","STATICDATA");

            // set element section
            MTStructure1D.ElementSetSection(i,Section1);

            // set element constitutive law
            MTStructure1D.ElementSetConstitutiveLaw(i,ConstLaw);

            // set static data for the elements integration points
            for (int theIP=0; theIP< MTStructure1D.ElementGetElementPtr(i)->GetNumIntegrationPoints(); theIP++)
            {
                NuTo::ConstitutiveStaticDataMoistureTransport *StaticData = MTStructure1D.ElementGetElementPtr(i)->GetStaticData(theIP)->AsMoistureTransport();
                StaticData->SetLastSorptionCoeff(MTStructure1D.ConstitutiveLawGetDesorptionCoefficients(ConstLaw));
                StaticData->SetActualSorptionCoeff(MTStructure1D.ConstitutiveLawGetDesorptionCoefficients(ConstLaw));
                StaticData->SetLastRelHumValue(InitialRelativeHumidity);
                StaticData ->SetSorptionHistoryDesorption(SorptionHistoryDesorption);
            }
        }



        // %%%%%%%%%%%%%%%%%%%%%%%%
        // Create Boundary Elements
        // %%%%%%%%%%%%%%%%%%%%%%%%


        // Add Nodes to boundary
        int nodeGroupBoundary = MTStructure1D.GroupCreate("NODES");
        MTStructure1D.GroupAddNodeCoordinateRange(nodeGroupBoundary, 0,         -1.e-6,         1.e-6);
        MTStructure1D.GroupAddNodeCoordinateRange(nodeGroupBoundary, 0,         L-1.e-6,        L + 1.e-6);

        // Group all elements with boundary nodes
        int elemGroupBoundary = MTStructure1D.GroupCreate("ELEMENTS");
        MTStructure1D.GroupAddElementsFromNodes(elemGroupBoundary, nodeGroupBoundary, false);

        // Create boundary elements
        ElementNodeNumbers.Resize(1);
        ElementNodeNumbers(0) = 0;
        MTStructure1D.BoundaryElementsCreate("BoundaryMoistureTransport1D", elemGroupBoundary, nodeGroupBoundary, NuTo::BoundaryCondition::NEUMANN_HOMOGENEOUS ,"ConstitutiveLawIp","StaticData");



        // %%%%%%%%%%%%%%%%%%%%%%
        // Test Getter and Setter
        // %%%%%%%%%%%%%%%%%%%%%%


        // Sorption Curves
        if (AdsorptionFit.GetDegree() != 3)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: AdsorptionFit.GetDegree() != 3");
        if (DesorptionFit.GetDegree() != 3)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: DesorptionFit.GetDegree() != 3");

        // Section
        if (MTStructure1D.SectionGetArea(Section1)!=Area)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.SectionGetArea(Section1)!=Area");

        // Constitutive Law
        if (MTStructure1D.ConstitutiveLawGetBoundarySurfaceMoistureTransportCoefficient(ConstLaw)!=BC_Surface_Moisture_Transfer)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetBoundarySurfaceMoistureTransportCoefficient(ConstLaw)!=BC_Surface_Moisture_Transfer");
        if (MTStructure1D.ConstitutiveLawGetKa(ConstLaw)!=Ka)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetKa(ConstLaw)!=Ka");
        if (MTStructure1D.ConstitutiveLawGetKd(ConstLaw)!=Kd)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetKd(ConstLaw)!=Kd");
        if (MTStructure1D.ConstitutiveLawGetMassExchangeRate(ConstLaw)!=MassExchangeRate)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetMassExchangeRate(ConstLaw)!=MassExchangeRate");
        if (MTStructure1D.ConstitutiveLawGetPorosity(ConstLaw)!=Porosity)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetPorosity(ConstLaw)!=Porosity");
        if (MTStructure1D.ConstitutiveLawGetVaporPhaseDiffusionCoefficient(ConstLaw)!=VaporPhaseDiffusionCoefficient)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetVaporPhaseDiffusionCoefficient(ConstLaw)!=VaporPhaseDiffusionCoefficient");
        if (MTStructure1D.ConstitutiveLawGetVaporPhaseDiffusionExponent(ConstLaw)!=VaporPhaseDiffusionExponent)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetVaporPhaseDiffusionExponent(ConstLaw)!=VaporPhaseDiffusionExponent");
        if (MTStructure1D.ConstitutiveLawGetVaporPhaseSaturationDensity(ConstLaw)!=VaporPhaseSaturationDensity)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetVaporPhaseSaturationDensity(ConstLaw)!=VaporPhaseSaturationDensity");
        if (MTStructure1D.ConstitutiveLawGetWaterPhaseDensity(ConstLaw)!=WaterPhaseDensity)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetWaterPhaseDensity(ConstLaw)!=WaterPhaseDensity");
        if (MTStructure1D.ConstitutiveLawGetWaterPhaseDiffusionCoefficient(ConstLaw)!=WaterPhaseDiffusionCoefficient)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetWaterPhaseDiffusionCoefficient(ConstLaw)!=WaterPhaseDiffusionCoefficient");
        if (MTStructure1D.ConstitutiveLawGetWaterPhaseDiffusionExponent(ConstLaw)!=WaterPhaseDiffusionExponent)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetWaterPhaseDiffusionExponent(ConstLaw)!=WaterPhaseDiffusionExponent");

        // Nodes
        for (unsigned int i=0; i<NNodes; i++)
        {
            if (MTStructure1D.NodeGetNodePtr(i)->GetRelativeHumidity(0) != InitialRelativeHumidity)
                throw NuTo::Exception(std::string("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.NodeGetNodePtr("+std::to_string(i)+")->GetRelativeHumidity(0) != InitialRelativeHumidity"));
            if (MTStructure1D.NodeGetNodePtr(i)->GetWaterPhaseFraction(0) != InitialWaterPhaseFraction)
                throw NuTo::Exception(std::string("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.NodeGetNodePtr("+std::to_string(i)+")->GetRelativeHumidity(0) != InitialRelativeHumidity"));
        }

    }
    catch(NuTo::Exception e)
    {
        std::cout << e.ErrorMessage() << std::endl;
        return 1;
    }
}
