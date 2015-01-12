#ifndef TESTCROSSSECTIONMODEL_HPP_
#define TESTCROSSSECTIONMODEL_HPP_

#include <cxxtest/TestSuite.h>          // Allows use of certain methods (include in any test)

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

// Must be included before any other cell_based headers
#include "CrossSectionModelArchiver.hpp"

#include "CylindricalHoneycombMeshGenerator.hpp"   // Helper class for generating mesh
#include "UniformDistributionSimpleWntCellCycleModel.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "CrossSectionModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CryptModelInteractionForce.hpp"
#include "BasementMembraneForce.hpp"
#include "CrossSectionCellsGenerator.hpp"
#include "CellPropertyRegistry.hpp"
#include "StromalCellMutationState.hpp"
#include "CellLabel.hpp"
#include "CrossSectionModelRandomCellKiller.hpp"
#include "Debug.hpp"
#include "PetscSetupAndFinalize.hpp"    // Needed for all tests that use Petsc (suite of data
										// structures and routines used in FE PDE solvers)

class TestCrossSectionModel : public AbstractCellBasedTestSuite
{	

public:

    void TestCrossSectionModelPeriodic() throw (Exception)
    {             	
    	RandomNumberGenerator::Instance()->Reseed(56);
    		
    	// Create mesh
    	unsigned cells_across = 10;      // Number of cells in each row
    	unsigned cells_up = 29;          // Number of rows of cells
    	unsigned ghosts = 2;             // Number of layers of ghost nodes
    	
    	// (unsigned numNodesAlongWidth, unsigned numNodesAlongLength, unsigned ghosts, double scaleFactor)
    	CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, ghosts);
    	Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
    	
    	// Aim to get 20-25 cells down each side (See Van Leeuwen 2005)
    	
    	double width = (double) cells_across - 0.5;
    	double height = ((double) cells_up - 1.0)*sqrt(3.0)*0.5;
    	
    	double centre_line = width*0.5;                         // The axis of symmetry down the centre of the crypt
    	double radius = width*0.25;                              // Distance from centre line to edge of crypt
    	    	
    	double crypt_edge_left = centre_line - radius;
    	double crypt_edge_right = centre_line + radius;
    	double crypt_base = height*0.25;
    	
    	c_vector<double,2> circle_centre;
    	circle_centre(0) = centre_line;
    	circle_centre(1) = crypt_base;
    	
    	// Getting the curvature for the crypt base    	      
//    	double curvature_crypt_base = 1.0/radius; 
    	    	
    	// Get location indices corresponding to real cells: eliminate any nodes that lie
    	// in the crypt lumen or outside the box and define them to be ghost nodes
    		
    	std::vector<unsigned> real_node_indices;
    	std::vector<unsigned> ghost_node_indices;
    	
    	for (unsigned i=0; i<p_mesh->GetNumAllNodes(); i++)
    	{    	
    		c_vector<double, 2> node_location;
    	    node_location(0) = p_mesh->GetNode(i)->rGetLocation()[0];
	    	node_location(1) = p_mesh->GetNode(i)->rGetLocation()[1];
	    	
	    	// Need to calculate the length of the vector between the node and the circle centre
	    	c_vector<double, 2> vector_to_circle_centre = node_location - circle_centre;
	    	
	    	double distance_to_centre = norm_2(vector_to_circle_centre);
	    	assert(distance_to_centre > 0);
	    	assert(!isnan(distance_to_centre));
	    	
	    	if ( ( (crypt_edge_left <= node_location(0)) && (node_location(0) <= crypt_edge_right)
	    	        && (node_location(1) >= crypt_base) ) || (distance_to_centre <= radius) || (node_location(0) < -1e-6) ||
	    	        (node_location(0) > width) || (node_location(1) < -1e-6) || (node_location(1) > height) )
	    	{
	    		ghost_node_indices.push_back(i);
	    	}
	    	else
	    	{
	    		real_node_indices.push_back(i);
	    	}
    	}
    	
    	// Create cells
    	std::vector<CellPtr> cells;
    	
    	CrossSectionCellsGenerator<UniformDistributionSimpleWntCellCycleModel> cells_generator;
    	cells_generator.Generate(cells, *p_mesh, crypt_edge_right, crypt_edge_left, height,
    			real_node_indices, ghost_node_indices, true);
    	
    	// Create tissue
    	MeshBasedCellPopulationWithGhostNodes<2> tissue(*p_mesh, cells, real_node_indices);
    	assert(tissue.GetNumRealCells() != 0);
		
    	// Make sure we have a Voronoi tessellation to begin with
		tissue.CreateVoronoiTessellation();
		
    	// Don't give mutant cells different properties to normal ones
    	tissue.SetDampingConstantMutant(tissue.GetDampingConstantNormal());
    		
    	// Create an instance of a Wnt concentration
    	WntConcentration<2>::Instance()->SetType(LINEAR);  // decreases from 1 to zero at height specified by CellPopulationConfig::mWntConcentrationParameter    	
    	WntConcentration<2>::Instance()->SetCellPopulation(tissue);
    	// Note: As the crypt equilibriates, it shrinks, so don't set this too high else you end up with far too many cells dividing
        WntConcentration<2>::Instance()->SetCryptLength(height + 3.0);	// Increase amount in () to lower the Wnt threshold for dividing cells
        																//  (i.e. increase the the height up the crypt that dividing cells reach)
        
        TS_ASSERT_EQUALS(WntConcentration<2>::Instance()->IsWntSetUp(), true);
                  
        tissue.SetOutputVoronoiData(true);
        tissue.SetOutputCellMutationStates(true);
        tissue.SetOutputCellProliferativeTypes(true);
        tissue.SetOutputCellAncestors(true);
        tissue.SetOutputCellVolumes(true);			// In 2D will output the cell areas
        tissue.SetOutputCellIdData(true);

    	// Create simulation from tissue
    	CrossSectionModel simulator(tissue, false, true, width, true);
    	// Width defined earlier, last true indicates to pin bottom cells
      
    	// Set maximum height for cells that should be pinned (use to pin just the bottom cells)
    	simulator.SetMaxHeightForPinnedCells(0.5);//(1.0-(crypt_base-radius-1));         // Adjust for the translation
    	simulator.SetRightHandBoundaryForPinnedCells(width);
    	simulator.SetLeftHandBoundaryForPinnedCells(0.0);
        simulator.SetOutputNodeVelocities(true);

    	// Create spring force law
    	CryptModelInteractionForce<2> meineke_force;  // Variable spring strengths
    	meineke_force.SetEpithelialStromalCellDependentSprings(true, 3, 3, 3, 3);

    	// Spring strengths dependent on the labelled states (healthy/labelled) of neighbouring cells: H-H, L-L, H-L, APC2-L
    	meineke_force.SetEdgeBasedSpringConstant(false);		// Spring constant dependent on edge length between two cells    	  	
    	meineke_force.SetPositionDependentSpringConstants(false);
    	meineke_force.SetUseOneWaySprings(true);

    	// Create basement membrane force law
    	BasementMembraneForce basement_membrane;
    	basement_membrane.SetBasementMembraneParameter(12.0);
    	basement_membrane.SetCryptBaseCurvature(0.3);
    	basement_membrane.SetPositionDependentMultiplier(false);	// Multiplier for basement membrane force at base of crypt

        simulator.AddForce(&meineke_force);
    	simulator.AddForce(&basement_membrane);
    	
    	// Create cell killer and pass in to crypt simulation
    	CrossSectionModelRandomCellKiller killing_cells(&tissue, true, 0.1);    	// (tissue, bool sloughOrifice)
    	simulator.AddCellKiller(&killing_cells);
    	
    	// Set output directory and end time
    	simulator.SetOutputDirectory("CrossSectionModelSimulation");
    	std::ofstream out("CrossSectionModelSimulationData.dat");

    	// Need a small timestep to keep the simulations stable
    	double normal_timestep = simulator.GetDt();		// This is 1/120 => every 30 seconds
    	double timestep_used = normal_timestep*0.5;
    	
    	simulator.SetDt(timestep_used);
    	simulator.SetSamplingTimestepMultiple(240);		// Every hour
    	simulator.SetEndTime(50.0);
    	
    	// Run simulation
    	simulator.Solve();

    	// Tidy up
    	WntConcentration<2>::Destroy();
        RandomNumberGenerator::Destroy();
    }

};

#endif /*TESTCROSSSECTIONMODEL_HPP_*/
