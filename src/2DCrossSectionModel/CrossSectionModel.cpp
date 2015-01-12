#include "CrossSectionModel.hpp"
#include "Debug.hpp"

//bool SortVectorAccordingToXCoordinate(const c_vector<double, 2> lhs, const c_vector<double, 2> rhs)
//{
//    return lhs[0] < rhs[0];
//}
//
//bool SortVectorAccordingToYCoordinate(const c_vector<double, 2> lhs, const c_vector<double, 2> rhs)
//{
//    return lhs[1] < rhs[1];
//}

bool CrossSectionModel::StoppingEventHasOccurred()
{

//	// /todo: Static cast the mesh in the constructor
//	Cylindrical2dMesh* p_mesh = dynamic_cast<Cylindrical2dMesh*>(&(mpStaticCastCellPopulation->rGetMesh()));
//
//	bool mismatched_elements = p_mesh->GetInstanceOfMismatchedBoundaryNodes();
//
//	return mismatched_elements;
	return false;
}


CrossSectionModel::CrossSectionModel(AbstractCellPopulation<2>& rCellPopulation,
                  bool deleteCellPopulationAndForceCollection,
                  bool initialiseCells,
                  double tissueWidth,
                  bool pinBottomCells)
    : CellBasedSimulation<2>(rCellPopulation,
                          deleteCellPopulationAndForceCollection,
                          initialiseCells),
      mUseJiggledBottomCells(false),
      mCellPopulationWidth(tissueWidth),
      mPinBottomCells(pinBottomCells)
{
  	assert(tissueWidth > 0.0);
    mpStaticCastCellPopulation = static_cast<MeshBasedCellPopulationWithGhostNodes<2>*>(&mrCellPopulation);
}


void CrossSectionModel::UseJiggledBottomCells()
{
    mUseJiggledBottomCells = true;
}


/*
 * Pinning the bottom tissue cells and applying rigid walls - overridden ApplyCellPopulationBoundaryConditions() method.
 */
void CrossSectionModel::ApplyCellPopulationBoundaryConditions(const std::vector< c_vector<double, 2> >& rOldLocations)
{
    // Iterate over all nodes associated with real cells to update their positions
    // according to any tissue boundary conditions

	if (mPinBottomCells == true)
	{
	    for (AbstractCellPopulation<2>::Iterator cell_iter = mrCellPopulation.Begin();
	    	 cell_iter != mrCellPopulation.End();
	    	 ++cell_iter)
	    {
	    	assert(mpStaticCastCellPopulation->IsCellAssociatedWithADeletedLocation(*cell_iter) == false);
	    	
	        // Get index of node associated with cell
	        unsigned node_index = mpStaticCastCellPopulation->GetNodeCorrespondingToCell(*cell_iter)->GetIndex();
	        assert(!(mpStaticCastCellPopulation->IsGhostNode(node_index)));
	        
	        // Get pointer to this node
	        Node<2>* p_node = mpStaticCastCellPopulation->GetNodeCorrespondingToCell(*cell_iter);
	
	        // Pin the bottom cells (only labelled)
	        
	        if ( (mpStaticCastCellPopulation->GetCellUsingLocationIndex(node_index)->GetMutationState()->IsType<StromalCellMutationState>())
	        		&& (p_node->rGetLocation()[1] < mMaxHeightForPinnedCells) )
	        {
	        	// Get old node location
	        	c_vector<double, 2> old_node_location = rOldLocations[node_index];
	
	        	// Return node to old location
	        	p_node->rGetModifiableLocation()[0] = old_node_location[0];
	        	p_node->rGetModifiableLocation()[1] = old_node_location[1];
	        }
	
	//        assert(p_node->rGetLocation()[0] >= mLeftHandBoundaryForPinnedCells);
	//    	assert(p_node->rGetLocation()[1] >= 0.0);
	//        assert(p_node->rGetLocation()[0] <=  mRightHandBoundaryForPinnedCells);

	         //Boundary conditions on the vertical edges for a test:
//			 if (p_node->rGetLocation()[0] < 0.0)
//			 {
//				p_node->rGetModifiableLocation()[0] = 0.0;
//			 }
//
//			 if (p_node->rGetLocation()[0] > 3.5)
//			 {
//				 p_node->rGetModifiableLocation()[0] = 4.5;
//			 }
	    }
	}
}


/* Defining the node locations when a cell divides. As we are considering the 2D cross section, we define
 * a cell to divide parallel to the vector that connects its immediate neighbours. So far, this is dependent
 * on a cell only having two healthy neighbours - will lead to problems if we get clustering.
 */

c_vector<double, 2> CrossSectionModel::CalculateCellDivisionVector(CellPtr pParentCell)
{		
    assert(pParentCell->GetMutationState()->IsType<StromalCellMutationState>()== false); // Should be epithelial
//    assert(pParentCell->GetMutationState()->IsType<WildTypeCellMutationState>()== true); // Should be epithelial or APC-2hit, for example
 
	double separation = mpStaticCastCellPopulation->GetMeinekeDivisionSeparation();

    unsigned parent_node_index = mpStaticCastCellPopulation->GetLocationIndexUsingCell(pParentCell);

    c_vector<double, 2> parent_coords = mpStaticCastCellPopulation->GetLocationOfCellCentre(pParentCell); //p_node->rGetLocation();
    c_vector<double, 2> daughter_coords, division_direction;

    std::set<unsigned> neighbours = GetNeighbouringNodeIndices(parent_node_index);

    std::vector<c_vector<double, 2> > epithelial_neighbours; 	// Vector of coordinates of only the neighbouring epithelial nodes
    
    for(std::set<unsigned>::iterator neighbour_iter = neighbours.begin();
			         neighbour_iter != neighbours.end();
				     ++neighbour_iter)
    {
    	// Want to check if the node is deleted, as it should be if that cell has been killed, and not consider that node in this case
    	if (!mpStaticCastCellPopulation->rGetMesh().GetNode(*neighbour_iter)->IsDeleted())
    	{
	    	CellPtr p_neighbour_cell = mpStaticCastCellPopulation->GetCellUsingLocationIndex(*neighbour_iter);
	
	    	// Only want neighbouring epithelial nodes that are not associated to dead cells
	    	if( (!mpStaticCastCellPopulation->IsGhostNode(*neighbour_iter)) 
	    			&& (p_neighbour_cell->GetMutationState()->IsType<StromalCellMutationState>()== false)
	    			&& (!(*p_neighbour_cell).IsDead()) )
	    	{
		   		c_vector<double,2> coord_of_cell = mpStaticCastCellPopulation->GetNode(*neighbour_iter)->rGetLocation();
	
		   		epithelial_neighbours.push_back(coord_of_cell);
	    	}
    	}
    }
    
    // Under normal circumstances, the epithelial cells in the monolayer should have only two epithelial 
    // neighbours. If they have more than this, something is wrong (e.g. you've stopped killing cells by anoikis) and so they will just
    // divide in a random direction. Update - mutant cells divide in a random direction too.
    
    if ( (epithelial_neighbours.size() != 2)  ||  (pParentCell->GetMutationState()->IsType<ApcTwoHitCellMutationState>()) )
    {
	    // Make a random direction vector of the required length
	    c_vector<double, 2> random_vector;

	    /* Pick a random direction and move the parent cell backwards by 0.5*separation
	     * in that direction and return the position of the daughter cell 0.5*separation
	     * forwards in that direction.
	     */
	    double random_angle = RandomNumberGenerator::Instance()->ranf();
	    random_angle *= 2.0*M_PI;

	    random_vector(0) = 0.5*separation*cos(random_angle);
	    random_vector(1) = 0.5*separation*sin(random_angle);

	    c_vector<double, 2> proposed_new_parent_coords = parent_coords - random_vector;
	    c_vector<double, 2> proposed_new_daughter_coords = parent_coords + random_vector;
	    
	    proposed_new_daughter_coords = parent_coords + 2.0*random_vector;
	    
	    while (proposed_new_daughter_coords(1) < 0.0)
	    {
	        random_angle = RandomNumberGenerator::Instance()->ranf();
	        random_angle *= 2.0*M_PI;

	        random_vector(0) = separation*cos(random_angle);
	        random_vector(1) = separation*sin(random_angle);
	        proposed_new_daughter_coords = parent_coords + random_vector;
	    }
	    
	    daughter_coords = proposed_new_daughter_coords;

	    assert(daughter_coords(1) >= 0.0); // to make sure dividing cells stay in the tissue
	    assert(parent_coords(1) >= 0.0);   // to make sure dividing cells stay in the tissue

	    // Set the parent to use this location
	    ChastePoint<2> parent_coords_point(parent_coords);

	    unsigned node_index = mpStaticCastCellPopulation->GetLocationIndexUsingCell(pParentCell);
	    mrCellPopulation.SetNode(node_index, parent_coords_point);    	
    }

    else 
    {    	
    	assert(epithelial_neighbours.size() == 2);

		// The direction of cell division shall be according to the vector that connects the neighbouring healthy, proliferating nodes
		// and will act in the upwards direction, i.e. the parent node will be placed above the daughter node
	   	// (if the nodes happen to be horizontal then it will just stick the daughter on one side or the other).

//		if(epithelial_neighbours[1][1] > epithelial_neighbours[0][1])
//		{
//			division_direction = epithelial_neighbours[1] - epithelial_neighbours[0];
//		}
//		else
//		{
//			division_direction = epithelial_neighbours[0] - epithelial_neighbours[1];
//		}

    	// Get the vector between those nodes (doesn't matter which way round)
		division_direction = epithelial_neighbours[0] - epithelial_neighbours[1];
    	
		double distance_between_nodes = norm_2(division_direction);
		assert(distance_between_nodes > 0);
		assert(!isnan(distance_between_nodes));

		division_direction /= distance_between_nodes;		// Normalise

	    // Move the parent cell forwards by 0.5*sep in the division direction and return the position of
	    // the daughter cell (0.5*sep backwards)

	    // Make a direction vector of the required length
	    c_vector<double, 2> direction_vector;

	    direction_vector = 0.5*separation*division_direction;

	    daughter_coords = parent_coords - direction_vector;
	    parent_coords += direction_vector;

	    assert(daughter_coords(1)>=0.0); // to make sure dividing cells stay in the tissue
	    assert(parent_coords(1)>=0.0);   // to make sure dividing cells stay in the tissue

	    // Set the parent to use this location
	    ChastePoint<2> parent_coords_point(parent_coords);

	    unsigned node_index = mpStaticCastCellPopulation->GetNodeCorrespondingToCell(pParentCell)->GetIndex();
	    mrCellPopulation.SetNode(node_index, parent_coords_point);   	
    }
   	
    return daughter_coords;
}


/* 
 * The following option is used to have random direction division.
 */
//
//c_vector<double, 2> CrossSectionModel::CalculateCellDivisionVector(CellPtr pParentCell)
//{
//	assert(pParentCell->GetMutationState()->IsType<WildTypeCellMutationState>());
//	
//    // Location of parent and daughter cells
//    c_vector<double, 2> parent_coords = mpStaticCastCellPopulation->GetLocationOfCellCentre(pParentCell);
//    c_vector<double, 2> daughter_coords;
//
//    // Get separation parameter
//    double separation = CellBasedConfig::Instance()->GetDivisionSeparation();
//
//    // Make a random direction vector of the required length
//    c_vector<double, 2> random_vector;
//
//    /*
//     * Pick a random direction and move the parent cell backwards by 0.5*separation
//     * in that direction and return the position of the daughter cell 0.5*separation
//     * forwards in that direction.
//     */
//
//    double random_angle = RandomNumberGenerator::Instance()->ranf();
//    random_angle *= 2.0*M_PI;
//
//    random_vector(0) = 0.5*separation*cos(random_angle);
//    random_vector(1) = 0.5*separation*sin(random_angle);
//
//    c_vector<double, 2> proposed_new_parent_coords = parent_coords - random_vector;
//    c_vector<double, 2> proposed_new_daughter_coords = parent_coords + random_vector;
//
////    if (   (proposed_new_parent_coords(1) >= 0.0)
////        && (proposed_new_daughter_coords(1) >= 0.0))
////    {
////        // We are not too close to the bottom of the tissue, so move parent
////        parent_coords = proposed_new_parent_coords;
////        daughter_coords = proposed_new_daughter_coords;
////    }
////    else
////    {
//    
//    proposed_new_daughter_coords = parent_coords + 2.0*random_vector;
//    
//    while (proposed_new_daughter_coords(1) < 0.0)
//    {
//        random_angle = RandomNumberGenerator::Instance()->ranf();
//        random_angle *= 2.0*M_PI;
//
//        random_vector(0) = separation*cos(random_angle);
//        random_vector(1) = separation*sin(random_angle);
//        proposed_new_daughter_coords = parent_coords + random_vector;
//    }
//    
//    daughter_coords = proposed_new_daughter_coords;
//
////    }
//
//    assert(daughter_coords(1) >= 0.0); // to make sure dividing cells stay in the tissue
//    assert(parent_coords(1) >= 0.0);   // to make sure dividing cells stay in the tissue
//
//    // Set the parent to use this location
//    ChastePoint<2> parent_coords_point(parent_coords);
//
//    unsigned node_index = mpStaticCastCellPopulation->GetLocationIndexUsingCell(pParentCell);
//    mrCellPopulation.SetNode(node_index, parent_coords_point);
//
//    return daughter_coords;
//}

/* 
 * Method to return the nodes connected to a particular node via the Delaunay
 * triangulation, excluding ghost nodes.
 */
std::set<unsigned> CrossSectionModel::GetNeighbouringNodeIndices(unsigned nodeIndex)
{
	assert(!(mpStaticCastCellPopulation->IsGhostNode(nodeIndex)));
	
	// Create a set of neighbouring node indices
	std::set<unsigned> neighbouring_node_indices;

    // Find the indices of the elements owned by this node
	std::set<unsigned> containing_elem_indices = mpStaticCastCellPopulation->GetNode(nodeIndex)->rGetContainingElementIndices();

    // Iterate over these elements
    for (std::set<unsigned>::iterator elem_iter = containing_elem_indices.begin();
         elem_iter != containing_elem_indices.end();
         ++elem_iter)
    {
        // Get all the nodes contained in this element
        // Don't want to include the current node
        unsigned neighbour_global_index;

        for (unsigned local_index=0; local_index<3; local_index++)
        {
            neighbour_global_index = mpStaticCastCellPopulation->rGetMesh().GetElement(*elem_iter)->GetNodeGlobalIndex(local_index);

            if( (neighbour_global_index != nodeIndex) && (!mpStaticCastCellPopulation->IsGhostNode(neighbour_global_index)) )
            {
            	neighbouring_node_indices.insert(neighbour_global_index);
            }
        }
    }

    return neighbouring_node_indices;
}

void CrossSectionModel::SetMaxHeightForPinnedCells(double maxHeightForPinnedCells)
{
	mMaxHeightForPinnedCells = maxHeightForPinnedCells;
}

double CrossSectionModel::GetMaxHeightForPinnedCells()
{
	return mMaxHeightForPinnedCells;
}

void CrossSectionModel::SetRightHandBoundaryForPinnedCells(double rightHandBoundaryForPinnedCells)
{
	mRightHandBoundaryForPinnedCells = rightHandBoundaryForPinnedCells;
}

void CrossSectionModel::SetLeftHandBoundaryForPinnedCells(double leftHandBoundaryForPinnedCells)
{
	mLeftHandBoundaryForPinnedCells = leftHandBoundaryForPinnedCells;
}

void CrossSectionModel::WriteVisualizerSetupFile()
{
    *mpVizSetupFile << "MeshWidth\t" << mpStaticCastCellPopulation->rGetMesh().GetWidth(0u) << "\n";// get farthest distance between nodes in the x-direction
}

std::string CrossSectionModel::GetDataOutputFile()
{
	return mDataOutputFile;
}

void CrossSectionModel::SetupSolve()
{
    // Sets up output file - name differently
	OutputFileHandler output_file_handler(this->mSimulationOutputDirectory + "/", false);
	mNumCellsResultsFile = output_file_handler.OpenOutputFile("results.numcells");
    *mpVizSetupFile << "NumCells\n";
}

void CrossSectionModel::PostSolve()
{
    SimulationTime* p_time = SimulationTime::Instance();
    *mNumCellsResultsFile << p_time->GetTime() << " " << mrCellPopulation.GetNumRealCells() << "\n";

}


void CrossSectionModel::AfterSolve()
{
    if ( mrCellPopulation.Begin() != mrCellPopulation.End() )  // there are any cells
    {
    	mNumCellsResultsFile->close();
    }

    CellBasedSimulation<2>::AfterSolve();
}

c_vector<double,2> CrossSectionModel::OutputVectorOfBasalLaminaResults(double basementMembraneParameter)
{	
    double num_transit_cells = 0.0;

    // We use an iterator provided by the tissue to loop over cells
    for (AbstractCellPopulation<2>::Iterator cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {
        if (cell_iter->GetCellCycleModel()->GetCellProliferativeType() == TRANSIT)
        {
            num_transit_cells++;
        }
    }
        
    c_vector<double,2> results;
    
    results[0] = basementMembraneParameter;
    results[1] = num_transit_cells;

    return results;
}

std::vector<c_vector<double,4> > CrossSectionModel::OutputVectorOfForceResults()
{
	std::vector<c_vector<double,4> > force_results;
	c_vector<double,4> vector_for_this_cell;
	double x, y, cell_area, node_index;
	
    // We use an iterator provided by the tissue to loop over cells
    for (AbstractCellPopulation<2>::Iterator cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {
    	// Only interested in the epithelial cells in the monolayer
    	
        if ( (cell_iter->GetCellCycleModel()->GetCellProliferativeType() == TRANSIT)
        		&& (cell_iter->GetMutationState()->IsType<StromalCellMutationState>()== false) )
        {
			node_index = mpStaticCastCellPopulation->GetNodeCorrespondingToCell(*cell_iter)->GetIndex();
			x = mpStaticCastCellPopulation->GetLocationOfCellCentre(*cell_iter)[0];
			y = mpStaticCastCellPopulation->GetLocationOfCellCentre(*cell_iter)[1];
			cell_area = mpStaticCastCellPopulation->GetVolumeOfVoronoiElement(node_index);
			// force = need to figure out how to calculate this for each node 

			vector_for_this_cell[0] = (double) node_index;
			vector_for_this_cell[1] = x;
			vector_for_this_cell[2] = y;
			vector_for_this_cell[3] = cell_area;
	//		vector_for_this_cell[4] = force;
								
			force_results.push_back(vector_for_this_cell);
        }
    }

    return force_results;
}


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CrossSectionModel)
