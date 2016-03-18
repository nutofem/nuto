// $Id: Brick8NcycleJump.cpp 2016-02-02 11:10:59Z vkindrac $

#include <iostream>
#include "nuto/base/Exception.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/FullVector.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/mechanics/timeIntegration/JumpDirect.h"
#include "nuto/mechanics/timeIntegration/TimeIntegrationBase.h"
#include <boost/filesystem.hpp>

// AmplFraction = AmplValue/rMeanValue
// returns a vector [TimeToLoadTillMeanValue,ValueAt(Time=rTime)]
std::pair<double,double> Loading(double rTime, double rMeanValue, double rAmplFraction, double rFrequency) {

	std::pair<double,double> Output;
	double Period, Value;
	const double pi = boost::math::constants::pi<double>();

	try {
		Period = (1/rAmplFraction)/(pi*rFrequency);
		if (Period <= 0) {
			throw "[NuTo::myNutoExamples::Brick8NcycleJump] negative loading time.";
		}
		Value = rMeanValue/pow(Period,4.);
		Value *= (2.*Period - rTime)*pow(rTime,3.);

		Output.first = Period; Output.second = Value;

		return Output;

	} catch (NuTo::Exception& e) {
        e.AddMessage("[NuTo::myNutoExamples::Brick8NcycleJump] something wrong in the Loading() function.");
        throw e;
	}

}

int main(int argc, char* argv[])
{
try
{
		// create the path to the result directory
    	boost::filesystem::path resultDir = boost::filesystem::current_path();
    	resultDir += std::string("/ResultsBrick8NcycleJump");

    	//delete result directory if it exists
    	if (boost::filesystem::exists(resultDir))    // does p actually exist?
    	{
    		if (boost::filesystem::is_directory(resultDir))      // is p a directory?
    		{
    			boost::filesystem::remove_all(resultDir);
    		}
    	}
    	// create result directory
    	boost::filesystem::create_directory(resultDir);

        // definitions
		// size of a cube
		double Width  = 1000.;
		double Height = 1000.;
		double Length = 1000.;
		// # of elements per side
		int NumElementsX = 2;
		int NumElementsY = 2;
		int NumElementsZ = 2;

        // applied boundary displacement: mean displacement, frequency, amplitude as a fraction of the mean displacement
        double BoundaryDisplacement = -2.;
        double Frequency(2.);
        double AmplitudeFraction(0.1);		// AmplValue/MeanValue

        // other values
        double ExtrapolationTolerance(0.2);
        double Period(Loading(0.,1.,AmplitudeFraction,Frequency).first);
        int NumberOfIncrementsToLoadTillMeanValue(20);

        // create structure
        NuTo::Structure myStructure(3);
        myStructure.SetNumTimeDerivatives(2);	//comment if static, uncomment if dynamic

        //**********************************************
        //          Material 1 (local viscoplastic)
        //**********************************************
        // create material law
    	int Material1 = myStructure.ConstitutiveLawCreate("DamageViscoPlasticityHardeningEngineeringStress");

    	double Density, youngsModulus, nu, TensileStrength, CompressiveStrength,
    		BiaxialCompressiveStrength, Viscosity, ViscosityExponent, DamageDistribution,
    		ViscoplasticYieldSurfaceOffset, FractureEnergy, HardeningValue, HardeningExponent;

    	Density = 2.3e-5;
    	youngsModulus = 29000.;
    	nu = 0.2;
    	TensileStrength = 6.74;
    	CompressiveStrength = 65.4;
    	BiaxialCompressiveStrength = 74;
    	Viscosity = 1.8;
    	ViscosityExponent = 2.6;
    	DamageDistribution = 0.5;
    	ViscoplasticYieldSurfaceOffset = -21.2;
    	FractureEnergy = 0.024;
    	HardeningValue = 5.8;
    	HardeningExponent = 1800;


    	myStructure.ConstitutiveLawSetParameterDouble(Material1,NuTo::Constitutive::eConstitutiveParameter::DENSITY, Density);
    	myStructure.ConstitutiveLawSetParameterDouble(Material1,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    	myStructure.ConstitutiveLawSetParameterDouble(Material1,NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, nu);
        myStructure.ConstitutiveLawSetParameterDouble(Material1,NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH, TensileStrength);
        myStructure.ConstitutiveLawSetParameterDouble(Material1,NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH, CompressiveStrength);
    	myStructure.ConstitutiveLawSetParameterDouble(Material1,NuTo::Constitutive::eConstitutiveParameter::BIAXIAL_COMPRESSIVE_STRENGTH, BiaxialCompressiveStrength);
    	myStructure.ConstitutiveLawSetParameterDouble(Material1,NuTo::Constitutive::eConstitutiveParameter::VISCOSITY, Viscosity);
    	myStructure.ConstitutiveLawSetParameterDouble(Material1,NuTo::Constitutive::eConstitutiveParameter::VISCOSITY_EXPONENT, ViscosityExponent);
    	myStructure.ConstitutiveLawSetParameterDouble(Material1,NuTo::Constitutive::eConstitutiveParameter::DAMAGE_DISTRIBUTION, DamageDistribution);
    	myStructure.ConstitutiveLawSetParameterDouble(Material1,NuTo::Constitutive::eConstitutiveParameter::VISCOPLASTIC_YIELD_SURFACE_OFFSET, ViscoplasticYieldSurfaceOffset);
    	myStructure.ConstitutiveLawSetParameterDouble(Material1,NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY, FractureEnergy);
    	myStructure.ConstitutiveLawSetParameterDouble(Material1,NuTo::Constitutive::eConstitutiveParameter::HARDENING_VALUE, HardeningValue);
    	myStructure.ConstitutiveLawSetParameterDouble(Material1,NuTo::Constitutive::eConstitutiveParameter::HARDENING_EXPONENT, HardeningExponent);
//		Set fatigue flag
    	myStructure.ConstitutiveLawSetParameterBool(Material1,NuTo::Constitutive::eConstitutiveParameter::FATIGUE_EXTRAPOLATION, true);

        int mySection1 = myStructure.SectionCreate("VOLUME");

        // create nodes
        NuTo::FullVector<double,Eigen::Dynamic> nodeCoordinates(3);
        int node = 0;
        for(int zCount = 0; zCount < NumElementsZ + 1; zCount++)
        {
            nodeCoordinates(2) = (double)zCount * Height/(double)NumElementsZ;
            for(int yCount = 0; yCount < NumElementsY + 1; yCount++)
            {
                nodeCoordinates(1) = (double)yCount * Width/(double)NumElementsY;
                for(int xCount = 0; xCount < NumElementsX + 1; xCount++)
                {
                    nodeCoordinates(0) = (double)xCount * Length/(double)NumElementsX;
                    //std::cout << "node: " << node << " coordinates: " << nodeCoordinates.GetValue(0,0) << "," << nodeCoordinates.GetValue(1,0) << "," << nodeCoordinates.GetValue(2,0) << std::endl;
                    myStructure.NodeCreate(node, nodeCoordinates);
                    node ++;
                }
            }
        }

        // interpolation type
        int InterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::BRICK3D);
        myStructure.InterpolationTypeAdd(InterpolationType, NuTo::Node::COORDINATES,    NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(InterpolationType, NuTo::Node::DISPLACEMENTS,  NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);


        // create elements
        NuTo::FullVector<int,Eigen::Dynamic> elementIncidence(8);
        int element = 0;
        for(int zCount = 0; zCount < NumElementsZ; zCount++)
        {
            for(int yCount = 0; yCount < NumElementsY; yCount++)
            {
                for(int xCount = 0; xCount < NumElementsX; xCount++)
                {
                    int node1 = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1) + xCount;
                    elementIncidence(0) = node1;
                    elementIncidence(1) = node1 + 1;
                    elementIncidence(2) = node1 + NumElementsX + 2;
                    elementIncidence(3) = node1 + NumElementsX + 1;
                    elementIncidence(4) = node1 + (NumElementsX + 1) * (NumElementsY + 1);
                    elementIncidence(5) = node1 + (NumElementsX + 1) * (NumElementsY + 1) + 1;
                    elementIncidence(6) = node1 + (NumElementsX + 1) * (NumElementsY + 1) + NumElementsX + 2;
                    elementIncidence(7) = node1 + (NumElementsX + 1) * (NumElementsY + 1) + NumElementsX + 1;
                    //std::cout << "element: " << element << " incidence: " << std::endl;
                    //elementIncidence.Info();
                    myStructure.ElementCreate(InterpolationType, elementIncidence, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP,NuTo::IpData::eIpDataType::STATICDATA);
                   	myStructure.ElementSetConstitutiveLaw(element,Material1);
            	    myStructure.ElementSetSection(element,mySection1);
                    element ++;
                }
            }
        }

        myStructure.ElementTotalConvertToInterpolationType();

        // boundary conditions
        NuTo::FullVector<double,Eigen::Dynamic> direction(3);
        direction(0)= 1;
        direction(1)= 0;
        direction(2)= 0;
        for(int zCount = 0; zCount < NumElementsZ + 1; zCount++)
        {
            for(int yCount = 0; yCount < NumElementsY + 1; yCount++)
            {
                int node = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1);
                //std::cout << "node: " << node << std::endl;
                myStructure.ConstraintLinearSetDisplacementNode(node, direction, 0.0);
            }
        }
        direction(0)= 0;
        direction(1)= 0;
        direction(2)= 1;
        myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);											// dynamic
        myStructure.ConstraintLinearSetDisplacementNode(NumElementsY * (NumElementsX + 1), direction, 0.0);   		// dynamic
        direction(0)= 0;
        direction(1)= 1;
        direction(2)= 0;
        myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);											// dynamic

        // define integration
        NuTo::JumpDirect myIntegrationScheme(&myStructure);			//	JumpDirect

        // apply nodes

        std::cout << "Displacement control" << std::endl;
        // create the set of the nodes experiencing the time-dependent active constraints
        int ActiveNodesBC = myStructure.GroupCreate("Nodes");

        for(int zCount = 0; zCount < NumElementsZ + 1; zCount++)
        {
        	//std::cout << zCount << std::endl;
        	for(int yCount = 0; yCount < NumElementsY + 1; yCount++)
        	{
        		int node = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1) + NumElementsX;
        		myStructure.GroupAddNode(ActiveNodesBC,node);  // myChange
            }
        }
        std::cout << "#ActiveNodeBC = " << myStructure.GroupGetNumMembers(ActiveNodesBC) << std::endl;

        // direction of the time-dependent constraint
        direction(0)= 1;
        direction(1)= 0;
        direction(2)= 0;

        // displacement time dependent factor
        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> displacementsBC(NumberOfIncrementsToLoadTillMeanValue + 1,2);
//    	NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> displacementsBC(2,2);

        for (int Incr = 0; Incr <= NumberOfIncrementsToLoadTillMeanValue; ++Incr) {
        	displacementsBC(Incr,0) = Period*Incr/NumberOfIncrementsToLoadTillMeanValue;	// time
        	displacementsBC(Incr,1)	= Loading(Period*Incr/NumberOfIncrementsToLoadTillMeanValue,
        			BoundaryDisplacement,AmplitudeFraction,Frequency).second;				// displ
		}
//      displacementsBC(0,0) = 0.;   // t=0
//      displacementsBC(0,1) = 0.;   // u(0)=0
//      displacementsBC(1,0) = 1.;   // t=1
//      displacementsBC(1,1) = 1.;   // u(1)= 1
//
//      displacementsBC.col(1) *= BoundaryDisplacement;

        // apply time-dependent constraints to the node set ActiveNodeBC
        int constraintLoading = myStructure.ConstraintLinearSetDisplacementNodeGroup(ActiveNodesBC,direction,0.);
        myIntegrationScheme.SetTimeDependentConstraint(constraintLoading, displacementsBC);

        // output and visualization for the reaction force
        myIntegrationScheme.AddResultGroupNodeForce("Out_NodeForceBC",ActiveNodesBC);


        // pseudo time-dependent analysis

        // start analysis
//        myIntegrationScheme.SetDynamic(true);
        myIntegrationScheme.SetDampingCoefficientMass(0.0);

        // set incrementation of a single cycle
   		myIntegrationScheme.SetHarmonicIncrementation(16);

        // apply harmonic excitation
        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> HarmonicConstraintFactor(1,3);	//	JumpDirect

       	HarmonicConstraintFactor(0,0) = AmplitudeFraction*BoundaryDisplacement;		//	set displacement amplitude
        HarmonicConstraintFactor(0,1) = Frequency;									//	set frequency
        HarmonicConstraintFactor(0,2) = 3;											// 	set number of cycles

        myIntegrationScheme.SetHarmonicExcitation(HarmonicConstraintFactor);		//	JumpDirect
   		myIntegrationScheme.SetHarmonicExtrapolation(true,&ExtrapolationTolerance);	// 	if .false., then straight forward integration of all cycles

        // specify solver parameters
        myIntegrationScheme.SetToleranceForce(1.e-4);
        myStructure.SetShowTime(true);
        myIntegrationScheme.SetAutomaticTimeStepping(true);
        myIntegrationScheme.SetTimeStep(0.05);
//      myIntegrationScheme.SetTimeStep(Period/NumberOfIncrementsToLoadTillMeanValue);
        myIntegrationScheme.SetMaxTimeStep(0.1);

        // output and visualization for the results
#ifdef ENABLE_VISUALIZE
        // visualize results
        int visualizationGroup = myStructure.GroupCreate(NuTo::Groups::eGroupId::Elements);
        myStructure.GroupAddElementsTotal(visualizationGroup);

        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::DISPLACEMENTS);
        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRAIN);
        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRESS);
        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::DAMAGE);
        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::TOTAL_INELASTIC_EQ_STRAIN);

    	boost::filesystem::path VTKFile;
    	VTKFile = resultDir;

    	VTKFile /= std::string("Brick8NcycleJumpElementOutput")+std::string(".vtk");
    	myStructure.ExportVtkDataFileElements(VTKFile.string());

#endif // ENABLE_VISUALIZE

        myIntegrationScheme.AddResultNodeDisplacements("Out_NodeDispl1",1);
        myIntegrationScheme.AddResultNodeDisplacements("Out_NodeDisplBC",14);
//        myIntegrationScheme.AddResultNodeAccelerations("Out_NodeAcc1,1);
//        myIntegrationScheme.AddResultNodeAccelerations("Out_NodeAccBC",14);
        myIntegrationScheme.AddResultElementIpData("Out_ElementStress",4,NuTo::IpData::eIpStaticDataType::ENGINEERING_STRESS);
        myIntegrationScheme.AddResultTime("Out_Time");
        myIntegrationScheme.SetResultDirectory(resultDir.string(),true);

        // solving

        myStructure.CalculateMaximumIndependentSets();
        try{
            myIntegrationScheme.Solve(5.); // solve untill time = 7
        }
        	catch(NuTo::Exception& e)
    	{
    		std::cout << e.ErrorMessage() << std::endl;
    	}

        myStructure.Info();

}
catch (NuTo::Exception& e)
{
	std::cout << "Error executing Brick8NcycleJump "<< std::endl;
	std::cout << e.ErrorMessage() << std::endl;
	return(-1);
}
return 0;
}

