#include "CrossSectionModelRandomCellKiller.hpp"
#include "Debug.hpp"

/* A cell killer class that applies
 * - Random apoptosis at the crypt orifice only for differentiated cells (to model sloughing due to waste / harsh environment)
 * - Anoikis (when cells pop out of the layer)
 */

CrossSectionModelRandomCellKiller::CrossSectionModelRandomCellKiller(AbstractCellPopulation<2>* pCrypt,
		bool sloughOrifice, double probabilityOfDeathInAnHour)
    : AbstractCellKiller<2>(pCrypt),
    mSloughOrifice(sloughOrifice),
    mProbabilityOfDeathInAnHour(probabilityOfDeathInAnHour),
    mCellsRemovedByAnoikis(0),
    mCellsRemovedRandomly(0)
{
    // Sets up output file
//	OutputFileHandler output_file_handler(mOutputDirectory + "AnoikisData/", false);
//	mAnoikisOutputFile = output_file_handler.OpenOutputFile("results.anoikis");
}

CrossSectionModelRandomCellKiller::~CrossSectionModelRandomCellKiller()
{
//    mAnoikisOutputFile->close();
}

bool CrossSectionModelRandomCellKiller::GetSloughOrifice() const
{
    return mSloughOrifice;
}

double CrossSectionModelRandomCellKiller::GetDeathProbabilityInAnHour() const
{
    return mProbabilityOfDeathInAnHour;
}

void CrossSectionModelRandomCellKiller::SetOutputDirectory(std::string outputDirectory)
{
	mOutputDirectory = outputDirectory;
}

std::string CrossSectionModelRandomCellKiller::GetOutputDirectory()
{
	return mOutputDirectory;
}

/*
 * Method to get the neighbouring nodes (excluding ghost nodes) of a particular node
 * Can then be used to identify the type of cells that surround a particular cell
 * Useful when you want to identify an epithelial cell that has popped out of the monolayer (no contact with stromal cells)
 */
std::set<unsigned> CrossSectionModelRandomCellKiller::GetNeighbouringNodeIndices(unsigned nodeIndex)
{
	// Create a set of neighbouring node indices
	std::set<unsigned> neighbouring_node_indices;

	// Need access to the mesh but can't get to it because the cell killer only owns a
	// pointer to an AbstractCellPopulation
    MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*> (this->mpCellPopulation);

	// Find the indices of the elements owned by this node
    std::set<unsigned> containing_elem_indices = p_tissue->rGetMesh().GetNode(nodeIndex)->rGetContainingElementIndices();

    // Iterate over these elements
    for (std::set<unsigned>::iterator elem_iter=containing_elem_indices.begin();
         elem_iter != containing_elem_indices.end();
	     ++elem_iter)
    {
	    // Get all the nodes contained in this element
	    unsigned neighbour_global_index;

	    for (unsigned local_index=0; local_index<3; local_index++)
	    {
	    	neighbour_global_index = p_tissue->rGetMesh().GetElement(*elem_iter)->GetNodeGlobalIndex(local_index);
	    	// Don't want to include the original node or ghost nodes
	    	if( (neighbour_global_index != nodeIndex) && (!p_tissue->IsGhostNode(neighbour_global_index)) )
            {
            	neighbouring_node_indices.insert(neighbour_global_index);
            }
	    }
    }
    return neighbouring_node_indices;
}

/** Method to determine if an epithelial cell has lost all contacts with the stromal cells below
 * TRUE if cell has popped up
 * FALSE if cell remains in the monolayer
 */
bool CrossSectionModelRandomCellKiller::HasCellPoppedUp(unsigned nodeIndex)
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*> (this->mpCellPopulation);

	bool has_cell_popped_up = false;	// Initialising

   	std::set<unsigned> neighbours = GetNeighbouringNodeIndices(nodeIndex);

   	unsigned num_stromal_neighbours = 0;

   	// Iterate over the neighbouring cells to check the number of differentiated cell neighbours

   	for(std::set<unsigned>::iterator neighbour_iter=neighbours.begin();
   							neighbour_iter != neighbours.end();
   							++neighbour_iter)
	{
		if ( (!p_tissue->IsGhostNode(*neighbour_iter)) && (p_tissue->GetCellUsingLocationIndex(*neighbour_iter)->GetMutationState()->IsType<StromalCellMutationState>()==true) )
   		{
			num_stromal_neighbours += 1;
		}
   	}

   	if(num_stromal_neighbours < 1)
   	{
   		has_cell_popped_up = true;
   	}

	return has_cell_popped_up;
}

/** A method to return a vector that indicates which cells should be killed by anoikis
 * and which by random apoptosis
 */
std::vector<c_vector<unsigned,3> > CrossSectionModelRandomCellKiller::RemoveByAnoikisOrRandomApoptosis()
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*> (this->mpCellPopulation);
    assert(p_tissue->GetVoronoiTessellation()!=NULL);	// This fails during archiving of a simulation as Voronoi stuff not archived yet

    std::vector<c_vector<unsigned,3> > cells_to_remove;
    c_vector<unsigned,3> individual_node_information;	// Will store the node index and whether to remove or not (1 or 0)

	for (AbstractCellPopulation<2>::Iterator cell_iter = p_tissue->Begin();
    	 cell_iter != p_tissue->End();
    	 ++cell_iter)
	{
		unsigned node_index = p_tissue->GetNodeCorrespondingToCell(*cell_iter)->GetIndex();
		assert((!p_tissue->IsGhostNode(node_index)));

		// Initialise
		individual_node_information[0] = node_index;
		individual_node_information[1] = 0;
		individual_node_information[2] = 0;

		// Examine each epithelial node to see if it should be removed by anoikis and then if it
		// should be removed by compression-driven apoptosis
		if (cell_iter->GetMutationState()->IsType<StromalCellMutationState>()==false)
		{
			// Determining whether to remove this cell by anoikis

			if(this->HasCellPoppedUp(node_index))
			{
				individual_node_information[1] = 1;
			}
		}

	    // We assume a constant time step
	    double death_prob_this_timestep = 1.0 - pow((1.0 - mProbabilityOfDeathInAnHour), SimulationTime::Instance()->GetTimeStep());
	    // Get height of this cell
	    double y = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter)[1];
		// Get the current height of the crypt
		double height = GetCryptHeightExtremes()[0];
		double base = GetCryptHeightExtremes()[1];
		// The height above which cells can be randomly killed
		double height_threshold_for_apoptosis = base + (height - base)*0.9;

	    if ((mSloughOrifice) && (!cell_iter->HasApoptosisBegun()
	    		&& (cell_iter->GetMutationState()->IsType<StromalCellMutationState>()==false)
	    		&& RandomNumberGenerator::Instance()->ranf() < death_prob_this_timestep)
	    		&& (y > height_threshold_for_apoptosis) )
	    {
	    	cell_iter->StartApoptosis();
//	    	individual_node_information[2] = 1;
	    }

		cells_to_remove.push_back(individual_node_information);
	}

	return cells_to_remove;
}


/*
 * Cell Killer that kills healthy cells that pop outwards and become detached from
 * the labelled tissue cells, i.e. removal by anoikis
 *
 * Also will remove differentiated cells at the orifice if mSloughOrifice is true
 */
void CrossSectionModelRandomCellKiller::TestAndLabelCellsForApoptosisOrDeath()
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*> (this->mpCellPopulation);
    assert(p_tissue->GetVoronoiTessellation()!=NULL);	// This fails during archiving of a simulation as Voronoi stuff not archived yet

    // Get the information at this timestep for each node index that says whether to remove by anoikis or random apoptosis
    std::vector<c_vector<unsigned,3> > cells_to_remove = this->RemoveByAnoikisOrRandomApoptosis();

    // Keep a record of how many cells have been removed at this timestep
    this->SetNumberCellsRemoved(cells_to_remove);
    this->SetLocationsOfCellsRemovedByAnoikis(cells_to_remove);

    // Need to avoid trying to kill any cells twice (i.e. both by anoikis or sloughing)
    // Loop over these vectors individually and kill any cells that they tell you to

    for (unsigned i=0; i<cells_to_remove.size(); i++)
    {
    	if ( (cells_to_remove[i][1] == 1) || (cells_to_remove[i][2] == 1) )
    	{
    		// Get cell associated to this node
    		CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(cells_to_remove[i][0]);
    		p_cell->Kill();
    	}
    }
}


void CrossSectionModelRandomCellKiller::SetNumberCellsRemoved(std::vector<c_vector<unsigned,3> > cellsRemoved)
{
	unsigned num_removed_by_anoikis = 0;
	unsigned num_removed_randomly = 0;

    for (unsigned i=0; i<cellsRemoved.size(); i++)
    {
    	if(cellsRemoved[i][1]==1)
    	{
    		num_removed_by_anoikis+=1;
    	}
    	if(cellsRemoved[i][2]==1)
    	{
    		num_removed_randomly+=1;
    	}
    }

    mCellsRemovedByAnoikis += num_removed_by_anoikis;
    mCellsRemovedRandomly += num_removed_randomly;
}

c_vector<unsigned,2> CrossSectionModelRandomCellKiller::GetNumberCellsRemoved()
{
	c_vector<unsigned,2> number_cells_removed;
	number_cells_removed[0] = mCellsRemovedByAnoikis;
	number_cells_removed[1] = mCellsRemovedRandomly;
	return number_cells_removed;
}

void CrossSectionModelRandomCellKiller::SetLocationsOfCellsRemovedByAnoikis(std::vector<c_vector<unsigned,3> > cellsRemoved)
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*> (this->mpCellPopulation);
	double x_location, y_location;
	c_vector<double, 3> time_and_location;

	// Need to use the node indices to store the locations of where cells are removed
    for (unsigned i=0; i<cellsRemoved.size(); i++)
    {
    	if (cellsRemoved[i][1] == 1)		// This cell has been removed by anoikis
    	{
			time_and_location[0] = SimulationTime::Instance()->GetTime();

			unsigned node_index = cellsRemoved[i][0];

			CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(node_index);
			x_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[0];
			y_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[1];

			time_and_location[1] = x_location;
			time_and_location[2] = y_location;

			mLocationsOfAnoikisCells.push_back(time_and_location);
    	}
    }
}

std::vector<c_vector<double,3> > CrossSectionModelRandomCellKiller::GetLocationsOfCellsRemovedByAnoikis()
{
	return mLocationsOfAnoikisCells;
}


/* Method to return the current coordinates of the crypt orifice and
 * crypt base - these can be used to accurately define the region of the
 * crypt base. (This will be the y-coordinate in 2D, or the z coordinate in 3D)
 * [0] - y-coordinate of orifice
 * [1] - y-coordinate of base
 */
c_vector<double,2> CrossSectionModelRandomCellKiller::GetCryptHeightExtremes()
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*> (this->mpCellPopulation);

    // Create a vector to store the y-coordinates of the lowest point of the crypt base and the highest point of the
    // crypt orifice
    c_vector<double,2> height_extremes;

    double max_height = 0.0;
    double min_height = DBL_MAX;

    double current_height_coordinate;

    // We iterate over all cells in the tissue, and deal only with those that are epithelial cells
	for (AbstractCellPopulation<2>::Iterator cell_iter = p_tissue->Begin();
    	 cell_iter != p_tissue->End();
    	 ++cell_iter)
    {
    	boost::shared_ptr<AbstractCellMutationState> p_state = cell_iter->GetMutationState();

	   	// Need these to not be labelled cells
    	if (p_state->IsType<StromalCellMutationState>()==false)
	   	{
	   		Node<2>* p_node = p_tissue->GetNodeCorrespondingToCell(*cell_iter);

   			current_height_coordinate = p_node->rGetLocation()[1];

	    	if (current_height_coordinate > max_height)
	    	{
	    		max_height = current_height_coordinate;
	    	}
	    	else if (current_height_coordinate < min_height)
	    	{
	    		min_height = current_height_coordinate;
	    	}
	    }
    }

    height_extremes[0] = max_height;
    height_extremes[1] = min_height;

    return height_extremes;
}

void CrossSectionModelRandomCellKiller::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<SloughOrifice>" << mSloughOrifice << "</SloughOrifice> \n";
    *rParamsFile << "\t\t\t<ProbabilityOfDeathInAnHour>" << mProbabilityOfDeathInAnHour << "</ProbabilityOfDeathInAnHour> \n";
    *rParamsFile << "\t\t\t<CellsRemovedByAnoikis>" << mCellsRemovedByAnoikis << "</CellsRemovedByAnoikis> \n";
    *rParamsFile << "\t\t\t<CellsRemovedRandomly>" << mCellsRemovedRandomly << "</CellsRemovedRandomly> \n";

    // Call direct parent class
    AbstractCellKiller<2>::OutputCellKillerParameters(rParamsFile);
}




#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CrossSectionModelRandomCellKiller)
