
#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "AbstractElement.hpp"
#include "Node.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "ChastePoint.hpp"
#include "Debug.hpp"
#include "AbstractConductivityModifier.hpp"

#include "PropagationPropertiesCalculator.hpp"

//#include "ChasteCuboid.hpp" //uncommend if using a circle shaped activation
//#include "ChasteEllipsoid.hpp"

#include <fstream>
#include <vector>
#include <set>
#include <cassert> // standard debugging tool (evaluation assertion)
#include <cmath> // for sqrt
#include "../src/CellICCBioPhy.hpp"
#include "../src/DummyCell.hpp"

// *************************** CELL FACTORY ************************************* //
class ICCCellFactory : public AbstractCardiacCellFactory<3>
{
private:
	std::set<unsigned> setICCNode;
public:
    ICCCellFactory(std::set<unsigned> iccNodes) //read list of ICC from the test
        : AbstractCardiacCellFactory<3>(), setICCNode(iccNodes)
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)

    {
        // Define pacemaker region
        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];

        double r = 200e-4; // set size of the radius
				double scale = 2e-4; // Set the same like for mesh.scale

        // Find all nodes which are ICC by reading the list
        //".find" tries to locate an element in a set. Returns iterator point to the sought-after element
        // otherwise .end() if not found
        unsigned index = pNode->GetIndex();
        if (setICCNode.find(index) != setICCNode.end())
        {
            // ICC Cells NOT IN pacemaker region
            CellICCBioPhy* cell = new CellICCBioPhy(mpSolver, mpZeroStimulus);
            // Parameters for CellICCBioPhy: DON'T change any of them
            cell->SetParameter("V_excitation", -55); // excitation value (threshold) default: -55mV
            cell->SetParameter("live_time", 10000); // time of resting potential
            cell->SetParameter("ode_time_step", 0.1); // Set the same as defined in HeartHandler
            cell->SetParameter("IP3Par", 0.00069); //
            cell->SetParameter("t_start", 600000); // Set larger than total simulation time

            // Active ICC Cells inside the pacemaker region (circle shaped over the whole z-depth)
            if  (((x-302.25*scale)*(x-302.25*scale)+(y-302.25*scale)*(y-302.25*scale)) < r*r)
            //  for a quadratic pacemaker region: if (x>=(306.968*scale-100e-4) && y>=(307.1210*scale-100e-4))
            {
                cell->SetParameter("t_start", 0); //Overwrites t_start
            }
            //Sets ICC
            return cell;
         }
        // All other cells which are not ICC or Bath
         else
         {
             CellDummyCellFromCellML* i_cell = new CellDummyCellFromCellML(mpSolver, mpZeroStimulus);
             return i_cell;
         }
    }
};

// *************************** CONDUCTIVITY MODIFIER ************************************* //
class ICCConductivityModifier : public AbstractConductivityModifier<3,3>
{
private:
    std::set<unsigned> setICCElement; //Copy list of Indexes with ICC attributes to this class
    c_matrix<double, 3,3> mSpecialMatrix;
    c_matrix<double, 3,3> mSpecialMatrix2;
    c_matrix<double, 3,3> mSpecialMatrix3;
    c_matrix<double, 3,3> mSpecialMatrix4;

public:
    ICCConductivityModifier(std::set<unsigned> elementIndexesICC)
        : AbstractConductivityModifier<3,3>(),
          setICCElement(elementIndexesICC),
          mSpecialMatrix( zero_matrix<double>(3,3) ),
          mSpecialMatrix2( zero_matrix<double>(3,3) ),
          mSpecialMatrix3( zero_matrix<double>(3,3) ),
          mSpecialMatrix4( zero_matrix<double>(3,3) )
          {
            //Conductivities for ICC and Dummy (Bath conductivity is set at HeartConfig::Instance())
            double intraICC = 0.024; // Intracellular ICC [mS/cm] = [S/m]*e-2
            double extraICC = 0.036; // Extracellular ICC
            double intraDummy = 0.02; // Intracellular Dummy
            double extraDummy = 0.02; // Extracellular Dummy

			// Introducing faster propagation in x direction set mSpecialMatrix(0,0) higher
			// e.g. intraICC + 0.05 and extraICC + 0.05
              mSpecialMatrix(0,0) = intraICC; // Intracellular ICC Conductivities
              mSpecialMatrix(1,1) = intraICC;
              mSpecialMatrix(2,2) = intraICC;

              mSpecialMatrix2(0,0) = extraICC; // Extracellular ICC Conductivities
              mSpecialMatrix2(1,1) = extraICC;
              mSpecialMatrix2(2,2) = extraICC;

              mSpecialMatrix3(0,0) = intraDummy; // Intracellular Dummy Conductivities
              mSpecialMatrix3(1,1) = intraDummy;
              mSpecialMatrix3(2,2) = intraDummy;

              mSpecialMatrix4(0,0) = extraDummy; // Extracellular Dummy Conductivities
              mSpecialMatrix4(1,1) = extraDummy;
              mSpecialMatrix4(2,2) = extraDummy;

          }

    c_matrix<double,3,3>& rCalculateModifiedConductivityTensor(unsigned elementIndex,
                                                               const c_matrix<double,3,3>& rOriginalConductivity,
                                                               unsigned domainIndex)
    {
        // Conductivities for ICC
        // depending on the `domainIndex` (intra/extracellular).
        if (setICCElement.find(elementIndex) != setICCElement.end())
        {
            if (domainIndex == 0) // domainIndex==0 implies intracellular
            {
                return mSpecialMatrix;
            }
            else //domain Index==1 implies extracellular
            {
                return mSpecialMatrix2;
            }

        }
        // Conductivities for Dummy cells
        else
        {
            if (domainIndex == 0) // domainIndex==0 implies intracellular
            {
                return mSpecialMatrix3;
            }
            else //domain Index==1 implies extracellular
            {
                return mSpecialMatrix4;
            }

        }
    }
};

// *************************** TEST ************************************* //
class TestICC3D : public CxxTest::TestSuite
{
public:
    void TestMesh3D() //throw(Exception)
    {
        // Read mesh created by TetGen
        // 'COARSENED' Mesh:
        //  TrianglesMeshReader<3,3> reader("projects/ICC_3D/test/Mesh/6465antrum_m05_ICC-MY_s8-29_surface_resempling10-box2.5-10pq1.5-10AV");
        // 'ORIGINAL' Mesh:
        //TrianglesMeshReader<3,3> reader("projects/mesh/ICC3D/icc-nonicc-bath.1");
	TrianglesMeshReader<3,3> reader("projects/mesh/ICC3D/6465antrum_m05_ICC-MY_s8-29_surface_resempling10-box2.5-10pq1.5-10AV");

        DistributedTetrahedralMesh<3,3> mesh; // Data shared among processes if run in parallel
        mesh.ConstructFromMeshReader(reader);

        // Scale mesh from 512x512x22pix to 303.64x303.64x22um by 1e-4
        mesh.Scale(2e-4, 2e-4, 2e-4); // double of original image size (because of the Jacobian and converging errors)
        mesh.SetMeshHasChangedSinceLoading();

        // Create list with ICC indexes
        std::set<unsigned> iccNodes;
        std::set<unsigned> elementIndexesICC;

        // Iterating trough all Elements in the mesh and assigning attributes, conductivities and saving all ICC nodes
        for (DistributedTetrahedralMesh<3,3>::ElementIterator iter = mesh.GetElementIteratorBegin();
                        iter != mesh.GetElementIteratorEnd(); ++iter)
        {
            // Read Attributes
            double attribute = iter->GetAttribute();

            // Set all small islands in the mesh to Dummy cells
            if (attribute != 1 && attribute != 2 && attribute != 3) // ICC=1(red), Dummy=2(green), Bath=3(blue)
            {
                iter->SetAttribute(2);
                attribute = iter->GetAttribute();
            }

            // Copy all nodes of the element to the elementIndexesICC list
            if (attribute == 1) // Check if ICC node
            {
                elementIndexesICC.insert(iter->GetIndex());
                for(int j = 0; j<=3; ++j)
                {
                    iccNodes.insert(iter->GetNodeGlobalIndex(j));
                }
            }
        }

        // Check that if we're in parallel no single process owns every element (to ensure that the conductivities
        // really are distributed).
        if (PetscTools::IsParallel())
               {
                   TS_ASSERT_DIFFERS( mesh.GetNumElements(), mesh.GetNumLocalElements() );
               }

 // Setting Attributes
        std::set<unsigned> background_ids;
        static unsigned background_id1 = 3; // Set Bath
        background_ids.insert(background_id1);
        // Uncomment for only ICC and Bath:
        //static unsigned background_id2 = 2; // Set Dummy as Bath
        //background_ids.insert(background_id2);

		// Set Attributes for ICC and Dummy
        std::set<unsigned> ICC_ids;
        static unsigned ICC_id1 = 1; // Set ICCs
        ICC_ids.insert(ICC_id1);
        static unsigned ICC_id2 = 2; // Set Dummy (Comment for ICC-Bath simulation without Dummy)
        ICC_ids.insert(ICC_id2);

        HeartConfig::Instance()->SetTissueAndBathIdentifiers(ICC_ids, background_ids); // tissue and bath ids

    	// Set Information for simulation
        HeartConfig::Instance()->SetSimulationDuration(100); //ms for one cycle 10,000
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1, 1, 5); //timesteps: ode, pde, printing
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(2000); // Ratio for each cell
        HeartConfig::Instance()->SetUseAbsoluteTolerance(2e-3); //Changed to get around the DIVERGED_ITS error default:2e-4
        HeartConfig::Instance()->SetCapacitance(3); // Membrane Capacitance
        HeartConfig::Instance()->SetBathConductivity(0.02); // Bath capacitance

		// Set outputfile name
        HeartConfig::Instance()->SetOutputDirectory("TestMesh3D_ICCDummyBath_100ms_dt5ms_v1");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer(true); // Set for visualizing with Meshlab
				//HeartConfig::Instance()->SetVisualizeWithCmgui(true);

        // Initialize cell_factory
        ICCCellFactory cell_factory(iccNodes);

        // Declare the problem class, `BidomainProblem<3>`
        BidomainProblem<3> bidomain_problem( &cell_factory, true); // true indicates we are solving a bath problem

        // When not used 'HeartConfig' for reading the mesh. Has to be called before Initialise
        bidomain_problem.SetMesh(&mesh);

        //  min/max voltage is printed as the simulation runs, useful for verifying that cells are stimulated
        bidomain_problem.SetWriteInfo();

        //Initialise
        bidomain_problem.Initialise(); // Has to be initialised before setting the conductivities

        // Set conductivities with the Conductivity Modifier for every element
        BidomainTissue<3>* p_bidomain_tissue = bidomain_problem.GetBidomainTissue();
        ICCConductivityModifier modifier(elementIndexesICC); // Initialise Conductivity Modifier
        p_bidomain_tissue->SetConductivityModifier( &modifier );

        //Solve
        bidomain_problem.Solve();

        // Switch on the output into log. For writing into file use "scons test-suite=<whatever> | tee output.txt"
        HeartEventHandler::Headings();
        HeartEventHandler::Report();
    }
};
