#include <boost/filesystem.hpp>

#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/MechanicsEnums.h"
#include "mechanics/groups/Group.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"

#include "math/EigenCompanion.h"
#include "mechanics/constraints/ConstraintCompanion.h"

#include "visualize/VisualizeEnum.h"

//Example: linear elastic 2D plane stress example solved with PGD (extra coordinates: E-factor (E_0=10000), load factor (F_0=500))
//Theory and results were presented in group meeting 01.02.2017 by Annika

struct Properties
{
    double youngsModulus = 10000.0;	//YoungsModulus fix value is regulated by E-modul factor
    double poissonRatio = 0.2;		//fix
    double density = 1.0;			//need mass matrix for norm
    double mt = 0.1;				//thickness fix

}properties;

struct SimpleMesh1D
{
	double ValueMin;
	double ValueMax;

	int numNodes;

	double Valuedelta;

}MeshE,MeshL;

struct SimpleMesh2D
{
	double ValueMin1;
	double ValueMax1;
	double ValueMin2;
	double ValueMax2;

	int numNodes;
	int numdofs;

}MeshX;

//Compute Integration of discrete function given in vector Data over 1D mesh given in Mesh numerically with simpsons rule
double SimpsonsIntegration(SimpleMesh1D Mesh, Eigen::VectorXd Data)
{
	double Value = 0.0;
	Value = Data(0)+Data(Mesh.numNodes-1);
	int k=1;
	for (int i=0; i<(Mesh.numNodes-2)/2; i++)
	{
		Value=Value+4*Data(k)+2*Data(k+1);
		k=k+2;

	}
	//if Mesh.numNodes not even one more
	if (Mesh.numNodes-2 != k-2)
	{
		Value=Value+4*Data(Mesh.numNodes-2);
	}

	Value=(Mesh.ValueMax - Mesh.ValueMin)/(Mesh.numNodes-1)*Value/3.0;

    return Value;
}

//set up NuTo structure for physical space (x,y,z)
void SetStructureX(NuTo::Structure &structure, Properties properties, SimpleMesh2D MeshX, int nodeGroupX, int nodeGroupConst, int nodeGroupLoad, int elemGroupLoad)
{

	structure.SetNumTimeDerivatives(0);
	structure.SetShowTime(false);
	structure.SetNumProcessors(1);

	//load mesh from gmesh file
	auto createdGroupIds=structure.ImportFromGmsh("./PGDPlateMeso2D.msh");

	int groupId = createdGroupIds[0].first;
	int myInterpolationType = createdGroupIds[0].second;

    //Interpolation type linear
    structure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeSetIntegrationType(myInterpolationType,NuTo::eIntegrationType::IntegrationType2D3NGauss3Ip);

	structure.SetVerboseLevel(10);
	structure.ElementGroupSetInterpolationType(groupId, myInterpolationType);
	structure.ElementTotalConvertToInterpolationType();

    //section
	auto mySectionX = NuTo::SectionPlane::Create(properties.mt, false);
	structure.ElementTotalSetSection(mySectionX);

    //create constitutive law
    int myMaterialX = structure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    structure.ConstitutiveLawSetParameterDouble(myMaterialX,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, properties.youngsModulus);
    structure.ConstitutiveLawSetParameterDouble(myMaterialX,NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, properties.poissonRatio);
    structure.ConstitutiveLawSetParameterDouble(myMaterialX,NuTo::Constitutive::eConstitutiveParameter::DENSITY, properties.density);
    structure.ElementTotalSetConstitutiveLaw(myMaterialX);
    structure.ConstitutiveLawInfo(10);

    //node groups (constraints and loads will be set in the PGD loop)
    //constraint at x=0
	int dir=0; //X
	double min=MeshX.ValueMin1-0.01;
	double max=MeshX.ValueMin1+0.01;
	structure.GroupAddNodeCoordinateRange(nodeGroupConst,dir,min,max);

	//Load at x=L
	min=MeshX.ValueMax1-0.01;
	max=MeshX.ValueMax1+0.01;
	structure.GroupAddNodeCoordinateRange(nodeGroupLoad,dir,min,max);
	structure.GroupAddElementsFromNodes(elemGroupLoad, nodeGroupLoad, false);

	//all nodes group
	int visualizationGroup = structure.GroupCreate(NuTo::eGroupId::Elements);
	structure.GroupAddElementsTotal(visualizationGroup);
	structure.GroupAddNodeFromElementGroupCoordinateRange(nodeGroupX,visualizationGroup,0,MeshX.ValueMin1,MeshX.ValueMax1);

}

int main()
{

	std::cout.precision(8);

	// results directory
    boost::filesystem::path resultDirectory = boost::filesystem::initial_path().string()+std::string("/ResultsPlateMeso2DPGD");
    //delete result directory
    if (boost::filesystem::exists(resultDirectory))// does p actually exist?
    {
        if (boost::filesystem::is_directory(resultDirectory))      // is p a directory?
        {
            boost::filesystem::remove_all(resultDirectory);
        }
    }
    // create result directory
    boost::filesystem::create_directory(resultDirectory);


	//INPUT:
    //structure X
    MeshX.ValueMin1=0.0; //xdirection // needed in SetStructure to create constraint groups
    MeshX.ValueMax1=32.0;
    MeshX.ValueMin2=0.0; //ydirection
    MeshX.ValueMax2=16;

    MeshX.numNodes=1061; //from gmesh input
    MeshX.numdofs=1061*2; //from gmesh input


    //structure E (E-factor, simple 1D Mesh)
    MeshE.ValueMin=0.5;
    MeshE.ValueMax=1.5;
    MeshE.numNodes=11;
    MeshE.Valuedelta=(MeshE.ValueMax-MeshE.ValueMin)/(MeshE.numNodes-1);
    std::vector<int> NodesE(MeshE.numNodes);
    Eigen::VectorXd CoordinatesE(MeshE.numNodes);
	for (int i=0; i<MeshE.numNodes; i++)
	{
		NodesE[i]=i;
		CoordinatesE(i) = i*MeshE.Valuedelta+MeshE.ValueMin;
	}
//    std::cout << "MeshE->Nodes\n";
//    for(int i : NodesE) std::cout << i << '\t';
//    std::cout << std::endl;
//    std::cout << "MeshE->Coordinates\n" << CoordinatesE.transpose() << std::endl;

    //structure L (load factor, simple 1D Mesh)
    MeshL.ValueMin=0.0;
    MeshL.ValueMax=1.0;
    MeshL.numNodes=11;
    MeshL.Valuedelta=(MeshL.ValueMax-MeshL.ValueMin)/(MeshL.numNodes-1);
    std::vector<int> NodesL(MeshL.numNodes);
    Eigen::VectorXd CoordinatesL(MeshL.numNodes);
	for (int i=0; i<MeshL.numNodes; i++)
	{
		NodesL[i]=i;
		CoordinatesL(i) = i*MeshL.Valuedelta+MeshL.ValueMin;
	}
//    std::cout << "MeshL->Nodes\n";
//    for(int i : NodesL) std::cout << i << '\t';
//    std::cout << std::endl;
//    std::cout << "MeshL->Coordinates\n" << CoordinatesL.transpose() << std::endl;

    NuTo::EigenCompanion::WriteToFile(CoordinatesE,(resultDirectory/"MeshE.dat").string(), "  ");
    NuTo::EigenCompanion::WriteToFile(CoordinatesL,(resultDirectory/"MeshL.dat").string(), "  ");

    //define load in separated form
    //n_1(x)=constant loadvalue in x direction, zero in y direction (exactly only at boundary, but will be integrated only over boundary, so doesn't matter)
    double loadvalue = 500;
    int counter =0;
    Eigen::VectorXd LoadX(MeshX.numdofs);
    for (int count=0; count<MeshX.numdofs/2; count++)
    {
    	LoadX(counter)=loadvalue; //x dof
    	LoadX(counter+1)=0;	//y dof
    	counter=counter+2;
    }
//    std::cout << "load in separated form" << std::endl;
//    std::cout << "n_1(x): " << LoadX.transpose() << std::endl;

    //n_2(la)=linear from ValueMin to ValueMax
    Eigen::VectorXd LoadL(MeshL.numNodes);
    for (int count=0; count<MeshL.numNodes; count++)
    {
    	LoadL(count)=MeshL.ValueMin+MeshL.Valuedelta*count;
    }
//    std::cout << "n_2(la): " << LoadL.transpose() << std::endl;

    //n_3(E)=constant(1.0) independent
    Eigen::VectorXd LoadE(MeshE.numNodes);
    for (int count=0; count<MeshE.numNodes; count++)
    {
    	LoadE(count)=1.0;
    }
//    std::cout << "n_3(E): " << LoadE.transpose() << std::endl;


    //PGD solver
  	int Nmax=3; //maximum number of PGD basis functions
   	int FPmax=5;//maximum number of fixed point iteration
   	double normX=1.0; //initialize norms
   	double normL=1.0;
   	double normE=1.0;

   	//PGD solution structure initialization
    Eigen::MatrixXd PGDX = Eigen::MatrixXd::Zero(MeshX.numdofs,Nmax);
    Eigen::MatrixXd PGDL = Eigen::MatrixXd::Zero(MeshL.numNodes,Nmax);
    Eigen::MatrixXd PGDE = Eigen::MatrixXd::Zero(MeshE.numNodes,Nmax);

    //new PGD basis functions initialization
    Eigen::VectorXd NewX = Eigen::VectorXd::Zero(MeshX.numdofs);
    Eigen::VectorXd NewL = Eigen::VectorXd::Zero(MeshL.numNodes);
    Eigen::VectorXd NewE = Eigen::VectorXd::Zero(MeshE.numNodes);

    //create structure ones to extract constraint nodes and mass matrix (do not change)
    NuTo::Structure myStructureTemp(2);
    int nodeGroupTemp = myStructureTemp.GroupCreate("Nodes");
    int nodeGroupConstTemp = myStructureTemp.GroupCreate("Nodes");
    int nodeGroupLoadTemp = myStructureTemp.GroupCreate("Nodes");
    int elemGroupLoadTemp = myStructureTemp.GroupCreate("Elements");
    properties.density=1.0;
    SetStructureX(myStructureTemp, properties, MeshX, nodeGroupTemp, nodeGroupConstTemp, nodeGroupLoadTemp, elemGroupLoadTemp);

    std::vector<int> NodesConst = myStructureTemp.GroupGetMemberIds(nodeGroupConstTemp);
    int NumNodesConst = myStructureTemp.GroupGetNumMembers(nodeGroupConstTemp);
//    std::cout << "Node with x=0: total " << NumNodesConst << " are\n";
//    for (int i : NodesConst) std::cout << i << '\t';
//    std::cout << std::endl;

//	compute mass matrix for normX
    auto hessian2 = myStructureTemp.BuildGlobalHessian2();
//  std::cout << "mass matrix \n " << hessian2 << std::endl;
    Eigen::MatrixXd mass = hessian2.JJ.ExportToFullMatrix();
//    std::cout << "mass matrix \n " << mass << std::endl;

    //PGD loop counter
    int countN=0;

    std::cout << "########### begin PGD loop ########### "<< std::endl;
    //PGD loop
    for (countN=0; countN<Nmax; countN++)
    {
    	std::cout << "########### enrichment step: "<< countN << std::endl;
        //initialize PGD basis functions with ones
    	Eigen::VectorXd EnrX(MeshX.numdofs);
    	Eigen::VectorXd EnrL(MeshL.numNodes);
    	Eigen::VectorXd EnrE(MeshE.numNodes);
       	EnrX.setOnes();
        EnrL.setOnes();
        EnrE.setOnes();
        //constraints have to be fulfilled for X
    	// constraint conditions in EnrX Node numbering according unconstraint stiffness matrix
    	for (int i=0;i<NumNodesConst;i++)
    	{
    		EnrX(NodesConst[i]*2)=0.; //x dof
    		EnrX(NodesConst[i]*2+1)=0.; //y dof
    	}
//        std::cout << "initialized PGD basis:\n EnrX: "<< EnrX.transpose()<< "\n EnrL: " << EnrL.transpose() << "\n EnrE: " << EnrE.transpose() <<std::endl;


//        //fix point iteration
        int countFP=0;
    	for (countFP=0; countFP<FPmax; countFP++)
    	{
    		std::cout << "########### fix point iteration step: " << countFP << std::endl;


    		/////////////////////////
    		//Step 1: X
    		/////////////////////////
    		//compute scalar factors
    		//\int_E S² E dE
    		Eigen::VectorXd ValuesE(MeshE.numNodes);
    		for (int count=0; count<MeshE.numNodes; count++)
    	    {
    	    	ValuesE(count)=EnrE(count)*EnrE(count)*CoordinatesE(count);
    	    }
//    		std::cout << "ValuesE: "<< ValuesE.transpose()<< std::endl;
    		double alpha1 = SimpsonsIntegration(MeshE, ValuesE);

    		//\int_E S n_3(E) dE
    		for (int count=0; count<MeshE.numNodes; count++)
    	    {
    	    	ValuesE(count)=EnrE(count)*LoadE(count);
    	    }
//    		std::cout << "ValuesE: "<< ValuesE.transpose()<< std::endl;
    		double beta1 = SimpsonsIntegration(MeshE,ValuesE);

    		//\int_L T T dL
    		Eigen::VectorXd ValuesL(MeshL.numNodes);
    		for (int count=0; count<MeshL.numNodes; count++)
    	    {
    	    	ValuesL(count)=EnrL(count)*EnrL(count);
    	    }
//    		std::cout << "ValuesL: "<< ValuesL.transpose()<< std::endl;
    		double alpha2 = SimpsonsIntegration(MeshL,ValuesL);

    		//\int_L T n_2(L) dL
    		for (int count=0; count<MeshL.numNodes; count++)
    	    {
    	    	ValuesL(count)=EnrL(count)*LoadL(count);
    	    }
//    		std::cout << "ValuesL: "<< ValuesL.transpose()<< std::endl;
    		double beta2 = SimpsonsIntegration(MeshL,ValuesL);

    		std::cout << "calculate scalar factors step1 (X): " << std::endl;
    		std::cout << "alpha1: " << alpha1  << " beta1: " << beta1 << " alpha2: " << alpha2 << " beta2: " << beta2 << std::endl;



    		//from last steps get initial displacements Uinit  (not really tested because problem convergend after first enrichment step!!)
    		Eigen::VectorXd Uinit = Eigen::VectorXd::Zero(MeshX.numdofs);
    		if (alpha1*alpha2>1e-3) //otherwise solution is zero anyway
    		{
				if (countN>0)
				{
					std::cout << "N>0: " << std::endl;
					for (int countNc=0; countNc<countN; countNc++)
					{
						for (int count=0; count<MeshE.numNodes; count++)
						{
							ValuesE(count)=EnrE(count)*PGDE(count,countNc)*CoordinatesE(count);
						}
						double alpha1i = SimpsonsIntegration(MeshE, ValuesE); // \int_E S E F3i dE

						for (int count=0; count<MeshL.numNodes; count++)
						{
							ValuesL(count)=EnrL(count)*PGDL(count,countNc);
						}
						double alpha2i = SimpsonsIntegration(MeshL, ValuesL); // \int_L T F2i dL

						std::cout << "alpha1i: " << alpha1i << " alpha2i: " << alpha2i << std::endl;
						Eigen::VectorXd PGDXNc = PGDX.col(countNc);
	//    				std::cout << "PGDXNc " << PGDXNc.transpose() << std::endl;
						Uinit=Uinit+PGDXNc*alpha1i*alpha2i/(alpha1*alpha2); //division of whole equation by alpha1*alpha2
//						std::cout << "Uinit: "<< Uinit.transpose()<< std::endl;
					}
				}
    		}

    		//create structureX and solve
    		properties.density=0.0;

    		NuTo::Structure myStructureX(2);
    		int nodeGroupX = myStructureX.GroupCreate("Nodes");
        	int nodeGroupConst = myStructureX.GroupCreate("Nodes");
        	int nodeGroupLoad = myStructureX.GroupCreate("Nodes");
        	int elemGroupLoad = myStructureX.GroupCreate("Elements");
    		SetStructureX(myStructureX, properties, MeshX, nodeGroupX,nodeGroupConst,nodeGroupLoad,elemGroupLoad);

    		//set constraints x=0 fixed in X and Y direction
            const auto& groupX0 = *myStructureX.GroupGetGroupPtr(nodeGroupConst)->AsGroupNode();
            myStructureX.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
                    NuTo::Constraint::Component(groupX0, {NuTo::eDirection::X, NuTo::eDirection::Y}));

    	    myStructureX.NodeBuildGlobalDofs();
    	    myStructureX.CalculateMaximumIndependentSets();
//    	    myStructureX.Info();

    		//set load conditions
    	    NuTo::NewmarkDirect myIntegrationScheme(&myStructureX);
    	    myIntegrationScheme.SetResultDirectory(resultDirectory.string(),false);

    	    double simulationTime = 1.0;

            double loadvaluePGD=0;
            if (alpha1*alpha2>1e-3)
            {
            	loadvaluePGD=loadvalue*beta1*beta2/(alpha1*alpha2);
            }
            std::cout << "loadvaluePGD " << loadvaluePGD << std::endl;

            Eigen::Vector2d Load({1.0,0.0});
            myStructureX.LoadSurfaceConstDirectionCreate2D(elemGroupLoad,nodeGroupLoad,Load);
            Eigen::Matrix<double, 2, 2> forceRHS;
//            forceRHS << 0, loadvaluePGD, simulationTime, loadvaluePGD;
            forceRHS << 0, 0, simulationTime, loadvaluePGD;
            myIntegrationScheme.SetTimeDependentLoadCase(0, forceRHS);

    		//solve FE problem for dimension X
			myIntegrationScheme.SetTimeStep(1.0);
			myIntegrationScheme.SetToleranceForce(1e-6);
			myIntegrationScheme.SetAutomaticTimeStepping(false);
			myIntegrationScheme.SetMaxTimeStep(1.0);
			myIntegrationScheme.SetVerboseLevel(0);
			myIntegrationScheme.SetShowTime(true);
			myIntegrationScheme.SetPerformLineSearch(false);

			myIntegrationScheme.Solve(simulationTime);

			//extract PGD function NewX from FE results
			Eigen::MatrixXd displacementsX;
			myStructureX.NodeGroupGetDisplacements(nodeGroupX, displacementsX);
//			std::cout << "DISP: " << displacementsX.transpose() << std::endl;
			int numrows=displacementsX.rows();
//			std::cout << "numrows: " << numrows << " numdofs/2 " << MeshX.numdofs/2<< std::endl;
			//reorganizate displacment matrix in vector according numbering of stiffness matrix without constraints!!
			counter = 0;
			for (int i=0; i<numrows; i++)
			{
				NewX(counter)=displacementsX(i,0)-Uinit(counter);
				NewX(counter+1)=displacementsX(i,1)-Uinit(counter+1);
				counter=counter+2;
			}
//			std::cout << "NewX: " << NewX.transpose() << std::endl;

			normX=NewX.transpose()* (mass * NewX);
			normX=sqrt(normX);
			std::cout << "normX: " << normX << std::endl;

			/////////////////////////
    		//Step 2: L
			/////////////////////////
    		//compute scalar factors (alpha1 and beta1 unchanged)
			//\int_X R'^2 dx		//
			//get stiffness without constraints to calculate integral (it's maybe possible to directly use only the unconstrained problem form the beginning)
            myStructureX.Constraints().RemoveAll();

			auto hessian0X = myStructureX.BuildGlobalHessian0();
//			std::cout << "After delete Constraints node numbering " << std::endl;
//			myStructureX.NodeInfo(10);
//			std::cout << "hessian0X" << hessian0X << std::endl;
			Eigen::MatrixXd stiffX = hessian0X.JJ.ExportToFullMatrix();
//			std::cout << "stiff\n "<< stiffX << std::endl;
			double alpha3=NewX.transpose()* (stiffX * NewX);

			//\int_boundaryX R n_1(x) dboundary_x
			auto ExternalLoad = myStructureX.BuildGlobalExternalLoadVector();
//			std::cout << "External Load Vector " << ExternalLoad.J << std::endl;
			Eigen::VectorXd ExLoad = ExternalLoad.J.Export();
			double beta3 = loadvalue* NewX.transpose()* ExLoad;


    		std::cout << "calculate scalar factors step2 L: " << std::endl;
    		std::cout << "alpha1: " << alpha1  << " beta1: " << beta1 << " alpha3: " << alpha3 << " beta3: " << beta3 << std::endl;

    		//compute new PGD basis function
			NewL=LoadL * beta3 * beta1;
//			std::cout << "NewL: " << NewL.transpose() << std::endl;

			//minus parts from last steps
    		if (countN>0)
    		{
    			std::cout << "N>0: " << std::endl;
    			for (int countNc=0; countNc<countN; countNc++)
    			{
    	    		for (int count=0; count<MeshE.numNodes; count++)
    	    	    {
    	    	    	ValuesE(count)=EnrE(count)*PGDE(count,countNc)*CoordinatesE(count);
    	    	    }
    				double alpha1i = SimpsonsIntegration(MeshE, ValuesE); // \int_E S E F3i dE

    				Eigen::VectorXd PGDXNc = PGDX.col(countNc);
    				double alpha3i = NewX.transpose() * (stiffX * PGDXNc);
    				alpha3i = alpha3i/properties.youngsModulus; // \int_X R' F1i' dx
    				Eigen::VectorXd PGDLNc = PGDL.col(countNc);

    				std::cout << "alpha1i: " << alpha1i << "alpha3i: " << alpha3i << std::endl;

    				NewL=NewL - PGDLNc * alpha3i *alpha1i;
    			}
    		}

    		if (alpha1*alpha3 > 0.001)
    		{
    			NewL = 1.0/(alpha1*alpha3) * NewL;
    		}
//    		std::cout << "NewL: " << NewL.transpose() << std::endl;
    		normL=MeshL.Valuedelta * NewL.transpose() * NewL;
    		normL=sqrt(normL);
    		std::cout << "normL: " << normL << std::endl;


    		/////////////////////////
    		//Step 3: E
    		/////////////////////////
    		//compute scalar factors (alpha3 beta 3 unchanged)
    		//\int_L T T dL
     		for (int count=0; count<MeshL.numNodes; count++)
    		{
    		  	ValuesL(count)=NewL(count)*NewL(count);
    		}
//    		std::cout << "ValuesL: "<< ValuesL.transpose()<< std::endl;
    		alpha2 = SimpsonsIntegration(MeshL,ValuesL);

    		//\int_L T n_2(L) dL
    		for (int count=0; count<MeshL.numNodes; count++)
    		{
    			ValuesL(count)=NewL(count)*LoadL(count);
    		}
//    		std::cout << "ValuesL: "<< ValuesL.transpose()<< std::endl;
    		beta2 = SimpsonsIntegration(MeshL,ValuesL);

    		std::cout << "calculate scalar factors step3 E: " << std::endl;
    		std::cout << "alpha2: " << alpha2  << " beta2: " << beta2 << " alpha3: " << alpha3 << " beta3: " << beta3 << std::endl;

    		//compute new PGD basis function
			NewE=LoadE * beta3 * beta2;
//			std::cout << "NewE: " << NewE.transpose() << std::endl;

			//minus parts from last steps
    		if (countN>0)
    		{
    			std::cout << "N>0: " << std::endl;
    			for (int countNc=0; countNc<countN; countNc++)
    			{
    	    		for (int count=0; count<MeshL.numNodes; count++)
    	    	    {
    	    	    	ValuesL(count)=NewL(count)*PGDL(count,countNc);
    	    	    }
    				double alpha2i = SimpsonsIntegration(MeshL, ValuesL); // \int_L T F2i dE

    				Eigen::VectorXd PGDXNc = PGDX.col(countNc);
    				double alpha3i = NewX.transpose()* (stiffX * PGDXNc);

    				Eigen::VectorXd PGDENc = PGDE.col(countNc);

    				std::cout << "alpha2i: " << alpha2i << "alpha3i: " << alpha3i << std::endl;

    				NewE=NewE - PGDENc * alpha2i *alpha3i;
    			}
    		}

    		if (alpha2*alpha3 > 0.001)
    		{
    			NewE = 1.0/(alpha2*alpha3) * NewE;
    		}
    		for (int count=0; count < MeshE.numNodes; count++)
    		{
    			NewE(count)=1/CoordinatesE(count) * NewE(count);
    		}
//    		std::cout << "NewE: " << NewE.transpose() << std::endl;
    		normE=MeshE.Valuedelta * NewE.transpose() * NewE;
    		normE=sqrt(normE);
    		std::cout << "normE: " << normE << std::endl;

    		//check convergence fix point iteration
    		//break if |NewX-EnrX| < tol and |NewL-EnrL|<tol and |newE-EnrE|<tol
    	    auto Delta1=(NewX-EnrX).cwiseAbs();
    	    auto Delta2=(NewL-EnrL).cwiseAbs();
    	    auto Delta3=(NewE-EnrE).cwiseAbs();
//    	    std::cout << "Delta1 " << Delta1.maxCoeff() << std::endl;
//    	    std::cout << "Delta2 " << Delta2.maxCoeff() << std::endl;
//    	    std::cout << "Delta3 " << Delta3.maxCoeff() << std::endl;
    	    if((Delta1.maxCoeff() > 1e-3) || (Delta2.maxCoeff() > 1e-3) || (Delta3.maxCoeff() > 1e-3))
    	    {
    	    	std::cout << "######### fix point iteration not converged" << std::endl;
				//update
				EnrX=NewX;
				EnrL=NewL;
				EnrE=NewE;
    	    }
    	    else
    	    {
    	    	std::cout << "######### fix point iteration converged " << std::endl;
    	    	std::cout << "step: "<< countFP << " error: "<< Delta1.maxCoeff()<< ", "<<Delta2.maxCoeff()<< ", "<<Delta3.maxCoeff() << std::endl;
    	    	break;
    	    }

    	} //END FIXPOINT (countFP)

    	//compute weighting factors
    	double normW=pow(normX*normL*normE,1.0/3.0);
    	double factorX=1.0;
    	double factorL=1.0;
    	double factorE=1.0;
    	if(normW>0.001)
    	{
    		factorX = normW/normX;
    		factorL = normW/normL;
    		factorE = normW/normE;
    	}

    	//update PGD functions
    	PGDX.col(countN)=NewX*factorX;
    	PGDL.col(countN)=NewL*factorL;
    	PGDE.col(countN)=NewE*factorE;


    	//check residuum here displacement based
    	double normU=normX*normX*normL*normL*normE*normE;
    	std::cout << "######### displacement norm: " << normU << std::endl;
    	if (normU < 1e-8)
    	{
    		std::cout << "######### Convergence reached (normU= " << normU << ")" << std::endl;
    		break;
    	}


    }

    std::cout << "########### PGD modes saved! computed number=" << countN << " modes total=" << Nmax << std::endl;
//    std::cout << "PGDX: \n" << PGDX << std::endl;
//    std::cout << "PGDL: \n" << PGDL << std::endl;
//    std::cout << "PGDE: \n" << PGDE << std::endl;


    //POSTPROCESSING
    std::cout << "########### POSTPROCESSING ########### "<< std::endl;
    // for PGDX change format to x y z values times number PGD modes
    Eigen::MatrixXd PGDXneu = Eigen::MatrixXd::Zero(MeshX.numNodes,Nmax*3);
    for (int k=0; k<countN; k++)
    {
    	 int counter1=0;
    	 for (int i=0; i<MeshX.numNodes; i++)
    	 {
    		 PGDXneu(i,(k+1)*3-3)=PGDX(counter1,k); //xValue
    		 PGDXneu(i,(k+1)*3-2)=PGDX(counter1+1,k); //yValue
    		 PGDXneu(i,(k+1)*3-1)=0.0; //zValue
    		 counter1=counter1+2;

    	 }
    }

//    save computed PGD modes
    NuTo::EigenCompanion::WriteToFile(PGDXneu,(resultDirectory/"PGDX.dat").string(), "  ");
    NuTo::EigenCompanion::WriteToFile(PGDL,(resultDirectory/"PGDL.dat").string(), "  ");
    NuTo::EigenCompanion::WriteToFile(PGDE,(resultDirectory/"PGDE.dat").string(), "  ");

//    post processing visualization of strains and stresses (cell data!!)
//    epsilon = sum PGDeps^i PGDL^i PGDE^i (cell data)
//    stress = sum PGDstress^i PGDL^i PGDE2^i (cell data)

//    PGDE2 = E * PGDE
    Eigen::MatrixXd PGDE2 = Eigen::MatrixXd::Zero(MeshE.numNodes,Nmax);
    for (int i=0; i<countN;i++)
    {
    	for (int k=0;k<MeshE.numNodes;k++)
        {
        	PGDE2(k,i)=CoordinatesE(k)*PGDE(k,i);
        }
	}
//    std::cout << "PGDE2: \n" << PGDE2 << std::endl; //in this case constant because PGDE ~ 1/E --> PGDE2 ~ E 1/E
    NuTo::EigenCompanion::WriteToFile(PGDE2,(resultDirectory/"PGDE2.dat").string(), "  ");

	//transform PGDE PGDE2 PGDL into cell data (one possible decomposition)
    Eigen::MatrixXd PGDEnew = Eigen::MatrixXd::Zero(MeshE.numNodes*2-1,Nmax); //for U
    Eigen::VectorXd CoordinatesEnew(MeshE.numNodes*2-1);
    Eigen::MatrixXd PGDE2cell = Eigen::MatrixXd::Zero(MeshE.numNodes*2-2,Nmax); //for sig
    Eigen::MatrixXd PGDEcell = Eigen::MatrixXd::Zero(MeshE.numNodes*2-2,Nmax); //for eps
    for (int i=0; i<countN;i++)
    {
    	int kn=1;
        PGDEnew(0,i)=PGDE(0,i);
        CoordinatesEnew(0)=CoordinatesE(0);
        for (int k=1;k<MeshE.numNodes;k++)
        {
        	PGDEnew(kn,i)=(PGDE(k,i)-PGDE(k-1,i))/2+PGDE(k-1,i);
           	PGDEnew(kn+1,i)=PGDE(k,i);
           	CoordinatesEnew(kn,i)=(CoordinatesE(k,i)-CoordinatesE(k-1,i))/2+CoordinatesE(k-1,i);
           	CoordinatesEnew(kn+1,i)=CoordinatesE(k,i);
           	kn=kn+2;
        }
        PGDEcell(0,i)=PGDE(0,i);
        PGDEcell(MeshE.numNodes*2-3,i)=PGDE(MeshE.numNodes-1,i);
        PGDE2cell(0,i)=PGDE2(0,i);
        PGDE2cell(MeshE.numNodes*2-3,i)=PGDE2(MeshE.numNodes-1,i);
        kn=1;
        for (int k=1;k<MeshE.numNodes-1;k++)
        {
        	PGDEcell(kn,i)=PGDE(k,i);
        	PGDEcell(kn+1,i)=PGDE(k,i);
        	PGDE2cell(kn,i)=PGDE2(k,i);
        	PGDE2cell(kn+1,i)=PGDE2(k,i);
        	kn=kn+2;
        }

	}
    NuTo::EigenCompanion::WriteToFile(PGDEnew,(resultDirectory/"PGDEnodeData.dat").string(), "  ");
    NuTo::EigenCompanion::WriteToFile(PGDEcell,(resultDirectory/"PGDEcellData.dat").string(), "  ");
    NuTo::EigenCompanion::WriteToFile(PGDE2cell,(resultDirectory/"PGDE2cellData.dat").string(), "  ");
    NuTo::EigenCompanion::WriteToFile(CoordinatesEnew,(resultDirectory/"CoordinatesE.dat").string(), "  ");

    Eigen::MatrixXd PGDLnew = Eigen::MatrixXd::Zero(MeshL.numNodes*2-1,Nmax); //for U
    Eigen::VectorXd CoordinatesLnew(MeshL.numNodes*2-1);
    Eigen::MatrixXd PGDLcell = Eigen::MatrixXd::Zero(MeshL.numNodes*2-2,Nmax); //for sig and eps
    for (int i=0; i<countN;i++)
    {
    	int kn=1;
        PGDLnew(0,i)=PGDL(0,i);
        CoordinatesLnew(0)=CoordinatesL(0);
        for (int k=1;k<MeshL.numNodes;k++)
        {
        	PGDLnew(kn,i)=(PGDL(k,i)-PGDL(k-1,i))/2+PGDL(k-1,i);
           	PGDLnew(kn+1,i)=PGDL(k,i);
           	CoordinatesLnew(kn,i)=(CoordinatesL(k,i)-CoordinatesL(k-1,i))/2+CoordinatesL(k-1,i);
           	CoordinatesLnew(kn+1,i)=CoordinatesL(k,i);
           	kn=kn+2;
        }
        PGDLcell(0,i)=PGDL(0,i);
        PGDLcell(MeshL.numNodes*2-3,i)=PGDL(MeshL.numNodes-1,i);
        kn=1;
        for (int k=1;k<MeshL.numNodes-1;k++)
        {
        	PGDLcell(kn,i)=PGDL(k,i);
            PGDLcell(kn+1,i)=PGDL(k,i);
            kn=kn+2;
        }
    }
//    std::cout << "PGDLnew " << PGDLnew << std::endl;
//    std::cout << "CoordinatesLnew " << CoordinatesLnew << std::endl;
//    std::cout << "PGDLcell " << PGDLcell << std::endl;
    NuTo::EigenCompanion::WriteToFile(PGDLnew,(resultDirectory/"PGDLnodeData.dat").string(), "  ");
    NuTo::EigenCompanion::WriteToFile(PGDLcell,(resultDirectory/"PGDLcellData.dat").string(), "  ");
    NuTo::EigenCompanion::WriteToFile(CoordinatesLnew,(resultDirectory/"CoordinatesL.dat").string(), "  ");



//    PGDeps und PGDstress by computing structure with given "displacements" PGDX

    for(int PGDmodenum=0;PGDmodenum<countN;PGDmodenum++)
	{
    	properties.density=0;

        NuTo::Structure myStructureX(2);
        int nodeGroupX = myStructureX.GroupCreate("Nodes");
        int nodeGroupConst = myStructureX.GroupCreate("Nodes");
        int nodeGroupLoad = myStructureX.GroupCreate("Nodes");
        int elemGroupLoad = myStructureX.GroupCreate("Elements");
        SetStructureX(myStructureX, properties, MeshX, nodeGroupX,nodeGroupConst,nodeGroupLoad,elemGroupLoad);

        //set displacement at each node to PGDX value
        for (int count=0;count<MeshX.numNodes;count++)
        {
        	Eigen::VectorXd disp(2);
        	for (int i=0;i<2;i++)
        	{
        		disp(i)=PGDXneu(count,(PGDmodenum+1)*3-3+i);
        	}
        	myStructureX.NodeSetDisplacements(count, disp);
//   		std::cout << "set disp at node " << count << " to " << disp.transpose() << std::endl;
        }


       	//add visualization
  		int visualizationGroup = myStructureX.GroupCreate(NuTo::eGroupId::Elements);
   		myStructureX.GroupAddElementsTotal(visualizationGroup);

   		myStructureX.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
  		myStructureX.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
   		myStructureX.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);


//    	Eigen::MatrixXd displacementsX;
//    	myStructureX.NodeGroupGetDisplacements(nodeGroupX, displacementsX);
//    	std::cout << "disp " << displacementsX << std::endl;

       	boost::filesystem::path resultFile;
       	resultFile = resultDirectory;
       	resultFile /= std::string("ElementOutput_END-")+std::to_string(PGDmodenum)+std::string(".vtk");
       	myStructureX.ExportVtkDataFileElements(resultFile.string());

	}

    //Check if solution is ok
    if (countN == 1 and (PGDE2(0,0) - 1.7909852) < 1e-7)
    {
    	std::cout << "All right!!" << std::endl;
    	return EXIT_SUCCESS;
    }
    else
    {
    	std::cout << "[PlateMeso2DPGD] result is not correct!" << std::endl;
    	return EXIT_FAILURE;
    }

}

