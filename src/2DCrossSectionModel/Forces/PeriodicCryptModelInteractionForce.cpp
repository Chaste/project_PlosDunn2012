#include "PeriodicCryptModelInteractionForce.hpp"
#include "Debug.hpp"
#include <cmath>
#include <list>
#include <fstream>

/**
 * To avoid warnings on some compilers, C++ style initialization of member
 * variables should be done in the order they are defined in the header file.
 */
template<unsigned DIM>
PeriodicCryptModelInteractionForce<DIM>::PeriodicCryptModelInteractionForce()
   : LinearSpringWithVariableSpringConstantsForce<DIM>(),
   mUseCellTypeDependentSprings(false),
   mTransitTransitMultiplier(DOUBLE_UNSET),
   mDifferentiatedDifferentiatedMultiplier(DOUBLE_UNSET),
   mTransitDifferentiatedMultiplier(DOUBLE_UNSET),
   mUseEpithelialStromalCellDependentSprings(false),
   mEpithelialEpithelialMultiplier(DOUBLE_UNSET),
   mStromalStromalMultiplier(DOUBLE_UNSET),
   mEpithelialStromalMultiplier(DOUBLE_UNSET),
   mApcTwoHitStromalMultiplier(DOUBLE_UNSET),
   mUseEdgeBasedSpringConstant(false),
   mUseOneWaySprings(false),
   mUsePositionDependentSpringConstants(false),
   mSpringConstantsMultiplier(DOUBLE_UNSET)
{
    // Sets up output file
//	OutputFileHandler output_file_handler("CurvatureData/", false);
//	mMeinekeOutputFile = output_file_handler.OpenOutputFile("results.curvature");
}

template<unsigned DIM>
PeriodicCryptModelInteractionForce<DIM>::~PeriodicCryptModelInteractionForce()
{
//    mMeinekeOutputFile->close();
}

template<unsigned DIM>
void PeriodicCryptModelInteractionForce<DIM>::SetCellTypeDependentSprings(bool useCellTypeDependentSprings,
		double transitTransitMultiplier,
		double differentiatedDifferentiatedMultiplier,
		double transitDifferentiatedMultiplier)
{
    mUseCellTypeDependentSprings = useCellTypeDependentSprings;
    mTransitTransitMultiplier = transitTransitMultiplier;
    mDifferentiatedDifferentiatedMultiplier = differentiatedDifferentiatedMultiplier;
    mTransitDifferentiatedMultiplier = transitDifferentiatedMultiplier;
}

template<unsigned DIM>
void PeriodicCryptModelInteractionForce<DIM>::SetEpithelialStromalCellDependentSprings(bool useEpithelialStromalCellDependentSprings,
		double epithelialEpithelialMultiplier,
		double stromalStromalMultiplier,
		double epithelialStromalMultiplier,
		double apcTwoHitStromalMultiplier)
{
	mUseEpithelialStromalCellDependentSprings = useEpithelialStromalCellDependentSprings;
	mEpithelialEpithelialMultiplier = epithelialEpithelialMultiplier;
    mStromalStromalMultiplier = stromalStromalMultiplier;
    mEpithelialStromalMultiplier = epithelialStromalMultiplier;
    mApcTwoHitStromalMultiplier = apcTwoHitStromalMultiplier;
}

template<unsigned DIM>
void PeriodicCryptModelInteractionForce<DIM>::SetEdgeBasedSpringConstant(bool useEdgeBasedSpringConstant)
{
    assert(DIM == 2);
    mUseEdgeBasedSpringConstant = useEdgeBasedSpringConstant;
}

template<unsigned DIM>
void PeriodicCryptModelInteractionForce<DIM>::SetUseOneWaySprings(bool useOneWaySprings)
{
	mUseOneWaySprings = useOneWaySprings;
}

template<unsigned DIM>
void PeriodicCryptModelInteractionForce<DIM>::SetPositionDependentSpringConstants(bool usePositionDependentSpringConstants, double springConstantsMultiplier)
{
	mUsePositionDependentSpringConstants = usePositionDependentSpringConstants;
	mSpringConstantsMultiplier = springConstantsMultiplier;
}

template<unsigned DIM>
double PeriodicCryptModelInteractionForce<DIM>::GetPositionDependentSpringConstants()
{
	return mSpringConstantsMultiplier;
}

template<unsigned DIM>
double PeriodicCryptModelInteractionForce<DIM>::VariableSpringConstantMultiplicationFactor(
											unsigned nodeAGlobalIndex,
											unsigned nodeBGlobalIndex,
											AbstractCellPopulation<DIM>& rCellPopulation,
											bool isCloserThanRestLength)
{
    double multiplication_factor = LinearSpringWithVariableSpringConstantsForce<DIM>::VariableSpringConstantMultiplicationFactor(nodeAGlobalIndex,
																																nodeBGlobalIndex,
																																rCellPopulation,
																																isCloserThanRestLength);

    MeshBasedCellPopulation<DIM>* p_tissue = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);
    assert(!(p_tissue->IsGhostNode(nodeAGlobalIndex)));
    assert(!(p_tissue->IsGhostNode(nodeBGlobalIndex)));

    CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
    CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (mUseCellTypeDependentSprings)
    {
        CellProliferativeType cell_A_type = p_cell_A->GetCellCycleModel()->GetCellProliferativeType();
        CellProliferativeType cell_B_type = p_cell_B->GetCellCycleModel()->GetCellProliferativeType();

        if (cell_A_type == cell_B_type)
        {
            if (cell_A_type == TRANSIT)
            {
            	multiplication_factor *= mTransitTransitMultiplier;
            }

            if (cell_A_type == DIFFERENTIATED)
            {
            	multiplication_factor *= mDifferentiatedDifferentiatedMultiplier;
            }
        }
        else
        {
        	multiplication_factor *= mTransitDifferentiatedMultiplier;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (mUseEpithelialStromalCellDependentSprings)
    {
        if ( (p_cell_A->GetMutationState()->IsType<StromalCellMutationState>()==false)			// If both not stromal => epithelial
			&& (p_cell_B->GetMutationState()->IsType<StromalCellMutationState>()==false) )
		{
			multiplication_factor *= mEpithelialEpithelialMultiplier;
		}
		else if ( (p_cell_A->GetMutationState()->IsType<StromalCellMutationState>()==true)		// If both stromal
				&& (p_cell_B->GetMutationState()->IsType<StromalCellMutationState>()==true) )
		{
			multiplication_factor *= mStromalStromalMultiplier;
		}
        else if ( ( (p_cell_A->GetMutationState()->IsType<StromalCellMutationState>()==false) && (p_cell_B->GetMutationState()->IsType<StromalCellMutationState>()==true) )
        	|| ( (p_cell_A->GetMutationState()->IsType<StromalCellMutationState>()==true) && (p_cell_B->GetMutationState()->IsType<StromalCellMutationState>()==false) ) )
        {
        	multiplication_factor *= mEpithelialStromalMultiplier;
        }
        else if ( ( (p_cell_A->GetMutationState()->IsType<ApcTwoHitCellMutationState>()) && (p_cell_A->GetMutationState()->IsType<StromalCellMutationState>()==false) && (p_cell_B->GetMutationState()->IsType<StromalCellMutationState>()==true) )
        		||  ( (p_cell_A->GetMutationState()->IsType<StromalCellMutationState>()==true) && (p_cell_B->GetMutationState()->IsType<ApcTwoHitCellMutationState>()) && (p_cell_B->GetMutationState()->IsType<StromalCellMutationState>()==false) ) )
        {
        	multiplication_factor *= mApcTwoHitStromalMultiplier;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (mUsePositionDependentSpringConstants)
    {
    	// Only need to worry about the connections between epithelial cells, as there is zero
    	// attractive force between epithelial and stromal cells

    	// Get the y-coordinate for the top of the crypt base
    	c_vector<double, 2> height_extremes = GetCryptHeightExtremes(rCellPopulation);

		double top_of_crypt_base = height_extremes(1) + (height_extremes(0) - height_extremes(1))*0.2;

		c_vector<double, 2> node_A_location = p_tissue->GetNode(nodeAGlobalIndex)->rGetLocation();
		c_vector<double, 2> node_B_location = p_tissue->GetNode(nodeBGlobalIndex)->rGetLocation();

        if ( (p_cell_A->GetMutationState()->IsType<StromalCellMutationState>()==false)			// If both not labelled => healthy
			&& (p_cell_B->GetMutationState()->IsType<StromalCellMutationState>()==false)
			&& ( (node_A_location[1] < top_of_crypt_base) || (node_B_location[1] < top_of_crypt_base) ))
        {
        	multiplication_factor *= mSpringConstantsMultiplier;
        }
    }

    return multiplication_factor;
}

/* Method to return the current coordinates of the crypt orifice and
 * crypt base - these can be used to accurately define the region of the
 * crypt base. (This will be the y-coordinate in 2D, or the z coordinate in 3D)
 * [0] - y or z-coordinate of orifice
 * [1] - y or z-coordinate of base
 */
template<unsigned DIM>
c_vector<double,2> PeriodicCryptModelInteractionForce<DIM>::GetCryptHeightExtremes(AbstractCellPopulation<DIM>& rCellPopulation)
{
    MeshBasedCellPopulation<DIM>* p_tissue = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

    // Create a vector to store the y-coordinates of the lowest point of the crypt base and the highest point of the
    // crypt orifice
    c_vector<double,2> height_extremes;

    double max_height = 0.0;
    double min_height = DBL_MAX;

    double current_height_coordinate;

    // We iterate over all cells in the tissue, and deal only with those that are epithelial cells
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
    	boost::shared_ptr<AbstractCellMutationState> p_state = cell_iter->GetMutationState();

	   	// Need these to not be labelled cells
	   	if ( (p_state->IsType<StromalCellMutationState>()==false) )
	   	{
	   		Node<DIM>* p_node = p_tissue->GetNodeCorrespondingToCell(*cell_iter);

   			current_height_coordinate = p_node->rGetLocation()[DIM-1];

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

template<unsigned DIM>
void PeriodicCryptModelInteractionForce<DIM>::RemoveDuplicates1D(std::vector<unsigned>& rVectorWithDuplicates)
{
    std::sort(rVectorWithDuplicates.begin(), rVectorWithDuplicates.end());
    rVectorWithDuplicates.erase(std::unique(rVectorWithDuplicates.begin(), rVectorWithDuplicates.end()), rVectorWithDuplicates.end());
}

/*
 * Method to determine whether an element contains ghost nodes
 */
template<unsigned DIM>
bool PeriodicCryptModelInteractionForce<DIM>::DoesElementContainGhostNodes(AbstractCellPopulation<DIM>& rCellPopulation, unsigned elementIndex)
{
	MeshBasedCellPopulation<DIM>* p_tissue = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

	bool element_contains_ghost_nodes = false;

	// Get a pointer to the element
	Element<DIM,DIM>* p_element = p_tissue->rGetMesh().GetElement(elementIndex);

	// ITERATE OVER NODES owned by this element
	for (unsigned local_index=0; local_index<DIM+1; local_index++)
	{
		if (p_tissue->IsGhostNode(p_element->GetNodeGlobalIndex(local_index)) == true)
		{
			element_contains_ghost_nodes = true;

		}
	}

	return element_contains_ghost_nodes;
}

/*
 * A method to return the number of elements that contain a particular node,
 * excluding those elements that have ghost nodes
 */
template<unsigned DIM>
unsigned PeriodicCryptModelInteractionForce<DIM>::GetNumContainingElementsWithoutGhostNodes(AbstractCellPopulation<DIM>& rCellPopulation, unsigned nodeIndex)
{
	MeshBasedCellPopulation<DIM>* p_tissue = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

    // Get pointer to the node
    Node<DIM>* p_node = p_tissue->GetNode(nodeIndex);
    assert(!(p_tissue->IsGhostNode(nodeIndex)));

    unsigned num_elements_with_no_ghost_nodes = 0;		// Initialise

    // Iterate over containing elements
    for (typename Node<DIM>::ContainingElementIterator iter = p_node->ContainingElementsBegin();
         iter != p_node->ContainingElementsEnd(); ++iter)
    {
        bool element_contains_ghost_nodes = false;
        Element<DIM,DIM>* p_element = p_tissue->rGetMesh().GetElement(*iter);

        // Iterate over nodes owned by this element
        for (unsigned local_index=0; local_index<DIM+1; local_index++)
        {
            if (p_tissue->IsGhostNode(p_element->GetNodeGlobalIndex(local_index)) == true)
            {
                element_contains_ghost_nodes = true;
                break; // I think this should break out of the inner for loop
            }
        }

        if (element_contains_ghost_nodes==false)
        {
            // This element contains no ghost nodes
            num_elements_with_no_ghost_nodes++;
        }
    }

    return num_elements_with_no_ghost_nodes;
}

/*
 * Method to return the nodes connected to a particular node via the Delaunay
 * triangulation, excluding ghost nodes.
 */
template<unsigned DIM>
std::set<unsigned> PeriodicCryptModelInteractionForce<DIM>::GetNeighbouringNodeIndices(AbstractCellPopulation<DIM>& rCellPopulation, unsigned nodeIndex)
{
	MeshBasedCellPopulation<DIM>* p_tissue = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

	assert(!(p_tissue->IsGhostNode(nodeIndex)));

	// Create a set of neighbouring node indices
	std::set<unsigned> neighbouring_node_indices;

    // Find the indices of the elements owned by this node
	std::set<unsigned> containing_elem_indices = p_tissue->GetNode(nodeIndex)->rGetContainingElementIndices();

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
            neighbour_global_index = p_tissue->rGetMesh().GetElement(*elem_iter)->GetNodeGlobalIndex(local_index);

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
template<unsigned DIM>
bool PeriodicCryptModelInteractionForce<DIM>::HasEpithelialCellDetachedFromBasementMembrane(AbstractCellPopulation<DIM>& rCellPopulation, unsigned nodeIndex)
{
	MeshBasedCellPopulation<DIM>* p_tissue = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

	bool has_cell_detached = false;	// Initialising

   	std::set<unsigned> neighbours = GetNeighbouringNodeIndices(rCellPopulation, nodeIndex);

   	unsigned num_stromal_neighbours = 0;

   	// Iterate over the neighbouring cells to check the number of differentiated cell neighbours

   	for(std::set<unsigned>::iterator neighbour_iter=neighbours.begin();
   							neighbour_iter != neighbours.end();
   							++neighbour_iter)
	{
   		boost::shared_ptr<AbstractCellMutationState> p_state = p_tissue->GetCellUsingLocationIndex(*neighbour_iter)->GetMutationState();
		if ( (!p_tissue->IsGhostNode(*neighbour_iter)) && (p_state->IsType<StromalCellMutationState>()==true) )
   		{
			num_stromal_neighbours += 1;
		}
   	}

   	if(num_stromal_neighbours < 1)
   	{
   		has_cell_detached = true;
   	}

	return has_cell_detached;
}

template<unsigned DIM>
void PeriodicCryptModelInteractionForce<DIM>::AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                                                                   AbstractCellPopulation<DIM>& rCellPopulation)
{
    // This method currently works only in 2d
    assert(DIM == 2);

    // Create a helper pointer
	MeshBasedCellPopulation<DIM>* p_tissue = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

	// Rather than iterating over the springs here, we firstly want to extend the mesh by creating
	// the image nodes, then loop over all edges to work out the forces acting on the real nodes.

	// Create a new, extended mesh by copying real nodes to form image nodes on either side, and create an accompanying extended cell population

	// These vectors will contain the real nodes and image nodes, and real cells and image cells, respectively
	std::vector<Node<DIM>*> extended_node_set;
    std::vector<CellPtr > extended_cell_set;

	// This vector will contain only the image nodes, and image cells, respectively
	std::vector<Node<DIM>*> image_node_set;
    std::vector<CellPtr> image_cells;

    // This vector will contain only the real cells
	std::vector<CellPtr> real_cells;

	unsigned num_real_nodes = rCellPopulation.GetNumRealCells();
	unsigned new_image_node_index = num_real_nodes;

	// The width of the extended mesh
	double extended_mesh_width = rCellPopulation.GetWidth(0); // + rCellPopulation.GetWidth(0)*0.5;

	PRINT_VARIABLE(extended_mesh_width);

    // We iterate over all cells in the population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Create the 'real' cell corresponding to this cell and add it to the vector of real cells
        CellPtr p_real_cell(new Cell(cell_iter->GetMutationState(), cell_iter->GetCellCycleModel(), false));
        real_cells.push_back(p_real_cell);

        // Get the node corresponding to this cell
    	Node<DIM>* p_node = p_tissue->GetNodeCorrespondingToCell(*cell_iter);
    	c_vector<double,DIM> real_node_location = p_node->rGetLocation();
    	unsigned real_node_index = p_node->GetIndex();

    	// Create a copy of this node and add it to the vector of all nodes
    	// (we push back all the original nodes first so that they are all kept together)
    	Node<DIM>* p_real_node = new Node<DIM>(real_node_index, real_node_location);
    	extended_node_set.push_back(p_real_node);

    	///\todo The code block below would need to be amended for 3d to cope with the z direction 
    	///      too and to make sure we copy from left to right and from front to back

    	// Compute the location of the image node corresponding to this node
        c_vector<double,DIM> image_node_location = real_node_location;
		if (real_node_location[0] >= rCellPopulation.GetWidth(0)*0.5) // Right-hand boundary node
		{
			image_node_location[0] -= extended_mesh_width;
		}
		else if (real_node_location[0] < rCellPopulation.GetWidth(0)*0.5)
		{
			image_node_location[0] += extended_mesh_width;
		}

		// Create the image node corresponding to this node and add it to the vector of image nodes
		Node<DIM>* p_image_node = new Node<DIM>(new_image_node_index, image_node_location);
        image_node_set.push_back(p_image_node);

		// Create the image cell corresponding to this cell and add it to the vector of image cells
		CellPtr p_image_cell(new Cell(cell_iter->GetMutationState(), cell_iter->GetCellCycleModel(), false));
		image_cells.push_back(p_image_cell);

		///\todo No set age method, which may matter when it comes to working out the rest length 
		///      of the spring for the force calculation
//		p_new_image_cell->SetAge();

        // Start from total number of real nodes and increment upwards
		new_image_node_index++;
    }

    // Now construct the vectors extended_node_set and extended_cell_set so that
    // the image nodes/cells are together at the end of each vector
    for (unsigned i=0; i<image_node_set.size(); i++)
    {
    	extended_node_set.push_back(image_node_set[i]);
    }
    for (unsigned i=0; i<real_cells.size(); i++)
    {
    	extended_cell_set.push_back(real_cells[i]);
    }
    for (unsigned i=0; i<image_cells.size(); i++)
    {
    	extended_cell_set.push_back(image_cells[i]);
    }

    // Check that the vectors extended_node_set and extended_cell_set are the correct size
    assert(extended_cell_set.size() == 2*num_real_nodes);
    assert(extended_node_set.size() == 2*num_real_nodes);

    // We now construct a mesh using extended_node_set...
    MutableMesh<DIM,DIM> extended_mesh(extended_node_set);

    // ...and, with this mesh and extended_cell_set, we create a MeshBasedCellPopulation
    MeshBasedCellPopulation<DIM>* p_extended_cell_population = new MeshBasedCellPopulation<DIM>(extended_mesh, extended_cell_set);

    ///\todo might we need to call Update() on extended_cell_population to ensure
    ///      that mMarkedSprings is correct?

	// Now loop over the extended mesh and calculate the force acting on real nodes
	// (using the edge iterator ensures that each edge is visited only once)
    for (typename MutableMesh<DIM,DIM>::EdgeIterator edge_iterator = extended_mesh.EdgesBegin();
         edge_iterator != extended_mesh.EdgesEnd();
         ++edge_iterator)
    {
        unsigned nodeA_global_index = edge_iterator.GetNodeA()->GetIndex();
        unsigned nodeB_global_index =  edge_iterator.GetNodeB()->GetIndex();

        c_vector<double, DIM> force = CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index, *p_extended_cell_population);

        // Now we make sure that we only apply the force to the real node and not the image node
        if ( (nodeA_global_index < num_real_nodes) &&  (nodeB_global_index < num_real_nodes) )
        {
			rForces[nodeB_global_index] -= force;
			rForces[nodeA_global_index] += force;
        }
        else if ( (nodeA_global_index >= num_real_nodes) &&  (nodeB_global_index < num_real_nodes) )
        {
			rForces[nodeB_global_index] -= force;
        }
        else if ( (nodeA_global_index < num_real_nodes) &&  (nodeB_global_index >= num_real_nodes) )
        {
        	rForces[nodeA_global_index] += force;
        }
    }
}


template<unsigned DIM>
c_vector<double, DIM> PeriodicCryptModelInteractionForce<DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
																						 unsigned nodeBGlobalIndex,
                                                                                         AbstractCellPopulation<DIM>& rCellPopulation)
{
    // We should only ever calculate the force between two distinct nodes
    assert(nodeAGlobalIndex != nodeBGlobalIndex);
    assert(!((static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))->IsGhostNode(nodeAGlobalIndex)));
    assert(!((static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))->IsGhostNode(nodeBGlobalIndex)));

    // Get the node locations
    c_vector<double, DIM> node_a_location = rCellPopulation.GetNode(nodeAGlobalIndex)->rGetLocation();
    c_vector<double, DIM> node_b_location = rCellPopulation.GetNode(nodeBGlobalIndex)->rGetLocation();

    // Get the unit vector parallel to the line joining the two nodes (assuming we don't have a cylindrical geometry here)
    c_vector<double, DIM> unit_difference = node_b_location - node_a_location;

    // Calculate the distance between the two nodes
    double distance_between_nodes = norm_2(unit_difference);
    assert(distance_between_nodes > 0);
    assert(!std::isnan(distance_between_nodes));

    unit_difference /= distance_between_nodes;

    /*
     * If mUseCutoffLength has been set, then there is zero force between
     * two nodes located a distance apart greater than mUseCutoffPoint.
     */
    if (this->mUseCutOffLength)
    {
        if (distance_between_nodes >= this->GetCutOffLength())
        {
            return zero_vector<double>(DIM); // c_vector<double,DIM>() is not guaranteed to be fresh memory
        }
    }

    // Calculate the rest length of the spring connecting the two nodes

    double rest_length = 1.0;

    CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
    CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

    double ageA = p_cell_A->GetAge();
    double ageB = p_cell_B->GetAge();

    assert(!std::isnan(ageA));
    assert(!std::isnan(ageB));

    /*
     * If the cells are both newly divided, then the rest length of the spring
     * connecting them grows linearly with time, until 1 hour after division.
     */
    if (ageA < this->mMeinekeSpringGrowthDuration && ageB < this->mMeinekeSpringGrowthDuration)
    {
        if (rCellPopulation.IsMeshBasedCellPopulation())
        {
            MeshBasedCellPopulation<DIM>* p_static_cast_cell_population = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

            std::pair<CellPtr,CellPtr> cell_pair = p_static_cast_cell_population->CreateCellPair(p_cell_A, p_cell_B);

            if (p_static_cast_cell_population->IsMarkedSpring(cell_pair))
            {
                // Spring rest length increases from a small value to the normal rest length over 1 hour
                double lambda = this->GetMeinekeDivisionRestingSpringLength();
                rest_length = lambda + (1.0 - lambda) * ageA/this->mMeinekeSpringGrowthDuration;
            }
            if (ageA + SimulationTime::Instance()->GetTimeStep() >= this->mMeinekeSpringGrowthDuration)
            {
                // This spring is about to go out of scope
                p_static_cast_cell_population->UnmarkSpring(cell_pair);
            }
        }
        else
        {
            // Spring rest length increases from mDivisionRestingSpringLength to normal rest length, 1.0, over 1 hour
            double lambda = this->GetMeinekeDivisionRestingSpringLength();
            rest_length = lambda + (1.0 - lambda) * ageA/this->mMeinekeSpringGrowthDuration;
        }
    }

    double a_rest_length = rest_length*0.5;
    double b_rest_length = a_rest_length;

    /*
     * If either of the cells has begun apoptosis, then the length of the spring
     * connecting them decreases linearly with time.
     */
    if (p_cell_A->HasApoptosisBegun())
    {
        double time_until_death_a = p_cell_A->GetTimeUntilDeath();
        a_rest_length = a_rest_length * time_until_death_a / p_cell_A->GetApoptosisTime();
    }
    if (p_cell_B->HasApoptosisBegun())
    {
        double time_until_death_b = p_cell_B->GetTimeUntilDeath();
        b_rest_length = b_rest_length * time_until_death_b / p_cell_B->GetApoptosisTime();
    }

    rest_length = a_rest_length + b_rest_length;
    assert(rest_length <= 1.0+1e-12);

    bool is_closer_than_rest_length = (distance_between_nodes - rest_length <= 0);

    // Although in this class the 'spring constant' is a constant parameter, in
    // subclasses it can depend on properties of each of the cells
    double multiplication_factor = VariableSpringConstantMultiplicationFactor(nodeAGlobalIndex, nodeBGlobalIndex, rCellPopulation, is_closer_than_rest_length);
    double spring_stiffness = this->GetMeinekeSpringStiffness();
    double overlap = distance_between_nodes - rest_length;

    /* Want to have one-way springs between epithelial and stromal nodes, so that there is only repulsion due to compression
     * of the spring, but no attraction due to extension
     */
    if ( (mUseOneWaySprings) && ( ( (p_cell_A->GetMutationState()->IsType<StromalCellMutationState>() == false) && (p_cell_B->GetMutationState()->IsType<StromalCellMutationState>() == true) )
    	    || ( (p_cell_A->GetMutationState()->IsType<StromalCellMutationState>() == true) && (p_cell_B->GetMutationState()->IsType<StromalCellMutationState>() == false) ) ) )
    {
        if (distance_between_nodes > rest_length)
        {
        	return zero_vector<double>(DIM); // c_vector<double,DIM>() is not guaranteed to be fresh memory
        }
    }

    if (rCellPopulation.IsMeshBasedCellPopulation())
    {
        return multiplication_factor * spring_stiffness * unit_difference * overlap;
    }
    else
    {
        // A reasonably stable simple force law
        if (distance_between_nodes > rest_length)
        {
            double alpha = 5;
            c_vector<double, DIM> temp = spring_stiffness * unit_difference * overlap * exp(-alpha * overlap);
            for (unsigned i=0; i<DIM; i++)
            {
                assert(!std::isnan(temp[i]));
            }
            return temp;
        }
        else
        {
            c_vector<double, DIM> temp = spring_stiffness * unit_difference * log(1 + overlap);
            for (unsigned i=0; i<DIM; i++)
            {
                assert(!std::isnan(temp[i]));
            }
            return temp;
        }
    }
}

template<unsigned DIM>
void PeriodicCryptModelInteractionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile <<  "\t\t\t<UseCellTypeDependentSprings>"<< mUseCellTypeDependentSprings << "</UseCellTypeDependentSprings> \n" ;
	*rParamsFile <<  "\t\t\t<TransitTransitMultiplier>"<< mTransitTransitMultiplier << "</TransitTransitMultiplier> \n" ;
	*rParamsFile <<  "\t\t\t<DifferentiatedDifferentiatedMultiplier>"<< mDifferentiatedDifferentiatedMultiplier << "</DifferentiatedDifferentiatedMultiplier> \n" ;
	*rParamsFile <<  "\t\t\t<TransitDifferentiatedMultiplier>"<< mTransitDifferentiatedMultiplier << "</TransitDifferentiatedMultiplier> \n" ;
	*rParamsFile <<  "\t\t\t<UseEpithelialStromalCellDependentSprings>"<< mUseEpithelialStromalCellDependentSprings << "</UseEpithelialStromalCellDependentSprings> \n" ;
	*rParamsFile <<  "\t\t\t<EpithelialEpithelialMultiplier>"<< mEpithelialEpithelialMultiplier << "</EpithelialEpithelialMultiplier> \n" ;
	*rParamsFile <<  "\t\t\t<StromalStromalMultiplier>"<< mStromalStromalMultiplier << "</StromalStromalMultiplier> \n" ;
	*rParamsFile <<  "\t\t\t<EpithelialStromalMultiplier>"<< mEpithelialStromalMultiplier << "</EpithelialStromalMultiplier> \n" ;
	*rParamsFile <<  "\t\t\t<ApcTwoHitStromalMultiplier>"<< mApcTwoHitStromalMultiplier << "</ApcTwoHitStromalMultiplier> \n" ;
	*rParamsFile <<  "\t\t\t<UseOneWaySprings>"<<  mUseOneWaySprings << "</mUseOneWaySprings> \n" ;
	*rParamsFile <<  "\t\t\t<UseEdgeBasedSpringConstant>"<<  mUseEdgeBasedSpringConstant << "</UseEdgeBasedSpringConstant> \n" ;

	// Call direct parent class
	LinearSpringWithVariableSpringConstantsForce<DIM>::OutputForceParameters(rParamsFile);
}



/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class PeriodicCryptModelInteractionForce<1>;
template class PeriodicCryptModelInteractionForce<2>;
template class PeriodicCryptModelInteractionForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PeriodicCryptModelInteractionForce)

