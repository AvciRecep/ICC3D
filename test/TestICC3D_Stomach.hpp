
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

#include <fstream>
#include <vector>
#include <set>
#include <cassert> // standard debugging tool (evaluation assertion)
#include <cmath> // for sqrt
#include "../src/CellICCBioPhy.hpp"

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
	
        double r = 0.1; // set size of the radius
	
        // Find all nodes which are ICC by reading the list
        //".find" tries to locate an element in a set. Returns iterator point to the sought-after element
        // otherwise .end() if not found
        unsigned index = pNode->GetIndex();
        if (setICCNode.find(index) != setICCNode.end())
        {
            // ICC Cells NOT IN pacemaker region
            CellICCBioPhy* cell = new CellICCBioPhy(mpSolver, mpZeroStimulus);
            // Parameters for CellICCBioPhy: DON'T change any of them
            cell->SetParameter("V_excitation", -65); // excitation value (threshold) default: -55mV
            cell->SetParameter("live_time", 10000); // time of resting potential
            cell->SetParameter("ode_time_step", 0.1); // Set the same as defined in HeartHandler
            cell->SetParameter("IP3Par", 0.00069); //
            cell->SetParameter("t_start", 600000); // Set larger than total simulation time

            // Active ICC Cells inside the pacemaker region (circle shaped over the whole z-depth)
            if  ( (x-0.392)*(x-0.392)+(y+15.072)*(y+15.072) < r*r)
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
	TrianglesMeshReader<3,3> reader("projects//mesh/ICC3D/stom_bath.1");
        DistributedTetrahedralMesh<3,3> mesh; // Data shared among processes if run in parallel
        mesh.ConstructFromMeshReader(reader);

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
            if (attribute != 1)
            {
                iter->SetAttribute(2);
                attribute = iter->GetAttribute();
            }

            // Copy all nodes of the element to the elementIndexesICC list
            if (attribute == 2) // Check if ICC node
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
        static unsigned background_id1 = 1; // Set Bath
        background_ids.insert(background_id1);

	// Set Attributes for ICC
        std::set<unsigned> ICC_ids;
        static unsigned ICC_id1 = 2; // Set ICCs
        ICC_ids.insert(ICC_id1);

        HeartConfig::Instance()->SetTissueAndBathIdentifiers(ICC_ids, background_ids); // tissue and bath ids

    	// Set Information for simulation
        HeartConfig::Instance()->SetSimulationDuration(10000); //ms for one cycle 10,000
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1, 1, 100); //timesteps: ode, pde, printing
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(2000); // Ratio for each cell
        HeartConfig::Instance()->SetUseAbsoluteTolerance(2e-3); //Changed to get around the DIVERGED_ITS error default:2e-4
        HeartConfig::Instance()->SetCapacitance(3); // Membrane Capacitance
        HeartConfig::Instance()->SetBathConductivity(0.02); // Bath capacitance

	// Set outputfile name
        HeartConfig::Instance()->SetOutputDirectory("TestMesh3D_Stomach_10s_dt100ms_v1");
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
