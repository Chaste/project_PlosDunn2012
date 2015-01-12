#include "CryptModelInteractionForce.hpp"
#include "Debug.hpp"
#include <cmath>
#include <list>
#include <fstream>

/**
 * To avoid warnings on some compilers, C++ style initialization of member
 * variables should be done in the order they are defined in the header file.
 */
template<unsigned DIM>
CryptModelInteractionForce<DIM>::CryptModelInteractionForce()
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
CryptModelInteractionForce<DIM>::~CryptModelInteractionForce()
{
//    mMeinekeOutputFile->close();
}

template<unsigned DIM>
void CryptModelInteractionForce<DIM>::SetCellTypeDependentSprings(bool useCellTypeDependentSprings,
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
void CryptModelInteractionForce<DIM>::SetEpithelialStromalCellDependentSprings(bool useEpithelialStromalCellDependentSprings,
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
void CryptModelInteractionForce<DIM>::SetEdgeBasedSpringConstant(bool useEdgeBasedSpringConstant)
{
    assert(DIM == 2);
    mUseEdgeBasedSpringConstant = useEdgeBasedSpringConstant;
}

template<unsigned DIM>
void CryptModelInteractionForce<DIM>::SetUseOneWaySprings(bool useOneWaySprings)
{
	mUseOneWaySprings = useOneWaySprings;
}

template<unsigned DIM>
void CryptModelInteractionForce<DIM>::SetPositionDependentSpringConstants(bool usePositionDependentSpringConstants, double springConstantsMultiplier)
{
	mUsePositionDependentSpringConstants = usePositionDependentSpringConstants;
	mSpringConstantsMultiplier = springConstantsMultiplier;
}

template<unsigned DIM>
double CryptModelInteractionForce<DIM>::GetPositionDependentSpringConstants()
{
	return mSpringConstantsMultiplier;
}

template<unsigned DIM>
double CryptModelInteractionForce<DIM>::VariableSpringConstantMultiplicationFactor(
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

		c_vector<double, 2> node_A_location = p_tissue->GetNode(nodeAGlobalIndex)->rGetLocation();;
		c_vector<double, 2> node_B_location = p_tissue->GetNode(nodeBGlobalIndex)->rGetLocation();;

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
c_vector<double,2> CryptModelInteractionForce<DIM>::GetCryptHeightExtremes(AbstractCellPopulation<DIM>& rCellPopulation)
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
void CryptModelInteractionForce<DIM>::RemoveDuplicates1D(std::vector<unsigned>& rVectorWithDuplicates)
{
    std::sort(rVectorWithDuplicates.begin(), rVectorWithDuplicates.end());
    rVectorWithDuplicates.erase(std::unique(rVectorWithDuplicates.begin(), rVectorWithDuplicates.end()), rVectorWithDuplicates.end());
}

/*
 * Method to determine whether an element contains ghost nodes
 */
template<unsigned DIM>
bool CryptModelInteractionForce<DIM>::DoesElementContainGhostNodes(AbstractCellPopulation<DIM>& rCellPopulation, unsigned elementIndex)
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
unsigned CryptModelInteractionForce<DIM>::GetNumContainingElementsWithoutGhostNodes(AbstractCellPopulation<DIM>& rCellPopulation, unsigned nodeIndex)
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
std::set<unsigned> CryptModelInteractionForce<DIM>::GetNeighbouringNodeIndices(AbstractCellPopulation<DIM>& rCellPopulation, unsigned nodeIndex)
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
bool CryptModelInteractionForce<DIM>::HasEpithelialCellDetachedFromBasementMembrane(AbstractCellPopulation<DIM>& rCellPopulation, unsigned nodeIndex)
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
void CryptModelInteractionForce<DIM>::AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                                                         AbstractCellPopulation<DIM>& rCellPopulation)
{	
//	MeshBasedCellPopulation<DIM>* p_tissue = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

//	// Find the epithelial nodes and apply a migration force
//	
//    for (typename MeshBasedCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
//         cell_iter != rCellPopulation.End(); ++cell_iter)
//    {    	
//		c_vector<double, 2> active_migration_direction;		
//		double active_migration_parameter = 1.5;
//		
//    	CellPropertyCollection collection = cell_iter->rGetCellPropertyCollection();
//    
//    	if(cell_iter->GetCellMutationState->IsType<StromalCellMutationState>() == false)		// If epithelial
//    	{
//    	    unsigned node_index = p_tissue->GetLocationIndexUsingCell(*cell_iter);
//    		
//    		// Get the vector between nearest epithelial neighbours
//    		
//    		std::set<unsigned> neighbours = GetNeighbouringNodeIndices(rCellPopulation, node_index);	// Might need to copy over this method
//
//    		std::vector<c_vector<double, 2> > epithelial_neighbours; 	// Vector of coordinates of only the neighbouring epithelial nodes
//
//    		for(std::set<unsigned>::iterator neighbour_iter = neighbours.begin();
//					         neighbour_iter != neighbours.end();
//						     ++neighbour_iter)
//		    {
//    			
//		    	if( (!p_tissue->IsGhostNode(*neighbour_iter)) && (p_tissue->GetCellUsingLocationIndex(*neighbour_iter)->GetCellMutationState->IsType<StromalCellMutationState>()==false) )
//		    	{
//			   		c_vector<double,2> coord_of_cell = p_tissue->GetNode(*neighbour_iter)->rGetLocation();
//
//			   		epithelial_neighbours.push_back(coord_of_cell);
//		    	}
//		    }
//    		
//    		if ( (epithelial_neighbours.size() != 2) )
//		    {
//			 	active_migration_direction[0] = 0.0;
//			 	active_migration_direction[1] = 0.0;
//		    }
//
//		    else 
//		    {    			
//	    		assert(epithelial_neighbours.size() == 2);
//	    		
//				// The direction of migration shall be according to the vector that connects the neighbouring epithelial nodes
//				// and will act in the upwards direction, i.e. the node will move to increase its y coordinate
//	
//				if(epithelial_neighbours[1][1] > epithelial_neighbours[0][1])
//				{
//					active_migration_direction = epithelial_neighbours[1] - epithelial_neighbours[0];
//				}
//				else
//				{
//					active_migration_direction = epithelial_neighbours[0] - epithelial_neighbours[1];
//				}
//	
//				double distance_between_nodes = norm_2(active_migration_direction);
//				assert(distance_between_nodes > 0);
//				assert(!isnan(distance_between_nodes));
//	
//				active_migration_direction /= distance_between_nodes;		// Normalise
//		    }
//
//			c_vector<double, 2> force_due_to_active_migration = active_migration_parameter*active_migration_direction;		
//    		
//			// Add the force due to the basal lamina to the forces acting on that epithelial node
//			rForces[node_index] += force_due_to_active_migration;
//    	}    
//    }
	
    /*
     * Here we want to deal with the case where we have a recent division, and hence a marked spring between two cells,
     * and avoid one of these new cells being forced out of the layer immediately and removed by anoikis
     * THIS NEEDS TO BE FIXED SO THAT YOU PUSH THE CELL THAT HAS POPPED OUT EITHER ANTICLOCKWISE OR CLOCKWISE, DEPENDING ON
     * WHICH SIDE OF THE PARENT CELL IT IS
     *
     *
     * 		X            X
     *       \          /
     *        \        /
     *         \      /
     * 			X    X
     *
     */
//    bool cell_detached_from_basement_membrane; 	// It's feeling melancholy
//    double cos_or_sin_of_one_over_root_two = 1.0/sqrt(2.0);
//
//	// Loop over all the springs and only work with those which are marked (between cells that have recently been born, so these
//	// are short springs)
//    for (typename MeshBasedCellPopulation<DIM>::SpringIterator spring_iterator=(static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))->SpringsBegin();
//            spring_iterator!=(static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))->SpringsEnd();
//            ++spring_iterator)
//    {
//        unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
//        unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
//
//        CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeA_global_index);
//        CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeB_global_index);
//
//        std::pair<CellPtr,CellPtr> cell_pair = p_tissue->CreateCellPair(p_cell_A, p_cell_B);
//
//		if ( (p_tissue->IsMarkedSpring(cell_pair)) && ( (this->HasEpithelialCellDetachedFromBasementMembrane(rCellPopulation, nodeA_global_index)==true)
//				|| (this->HasEpithelialCellDetachedFromBasementMembrane(rCellPopulation, nodeB_global_index)==true) ) )
//		{
//			PRINT_VARIABLE(SimulationTime::Instance()->GetTime());
//			TRACE("Applying torque force");
//			// Need to apply this force clockwise or anticlockwise depending on whether the detached cell lies to the left or the right
//			// of the attached cell
//
//		    double distance_between_nodes = norm_2(p_tissue->rGetMesh().GetVectorFromAtoB(p_tissue->GetNode(nodeA_global_index)->rGetLocation(),p_tissue->GetNode(nodeB_global_index)->rGetLocation()));
//
//			// Get vector between the springs and use the direction that is 45 degrees to this
//			c_vector<double, 2> vector_nodeA_to_nodeB = p_tissue->rGetMesh().GetVectorFromAtoB(p_tissue->GetNode(nodeA_global_index)->rGetLocation(),p_tissue->GetNode(nodeB_global_index)->rGetLocation());
//			c_vector<double, 2> vector_nodeB_to_nodeA = p_tissue->rGetMesh().GetVectorFromAtoB(p_tissue->GetNode(nodeB_global_index)->rGetLocation(),p_tissue->GetNode(nodeA_global_index)->rGetLocation());
//
//			// We want to rotate these vectors so that we apply a clockwise torque to each node (hence we need to rotate each anticlockwise)
//
//			c_vector<double, 2> torque_force_direction_A, torque_force_direction_B;
//
//			torque_force_direction_A[0] = cos_or_sin_of_one_over_root_two*vector_nodeA_to_nodeB[0] - cos_or_sin_of_one_over_root_two*vector_nodeA_to_nodeB[1];
//			torque_force_direction_A[1] = cos_or_sin_of_one_over_root_two*vector_nodeA_to_nodeB[0] + cos_or_sin_of_one_over_root_two*vector_nodeA_to_nodeB[1];
//
//			torque_force_direction_B[0] =  cos_or_sin_of_one_over_root_two*vector_nodeB_to_nodeA[0] - cos_or_sin_of_one_over_root_two*vector_nodeB_to_nodeA[1];
//			torque_force_direction_B[1] =  cos_or_sin_of_one_over_root_two*vector_nodeB_to_nodeA[0] + cos_or_sin_of_one_over_root_two*vector_nodeB_to_nodeA[1];
//
//			torque_force_direction_A /= distance_between_nodes;
//			torque_force_direction_B /= distance_between_nodes;
//
//			// What is the torque force? i.e. the magnitude of it? Need to think about that
//
//			double force_magnitude = fabs(distance_between_nodes - 1.0);	// Difference between current spring length and the full rest length of a spring
//
//			rForces[nodeA_global_index] = force_magnitude*torque_force_direction_A;
//			rForces[nodeB_global_index] = force_magnitude*torque_force_direction_B;
//		}
//	}

	for (typename MeshBasedCellPopulation<DIM>::SpringIterator spring_iterator=(static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))->SpringsBegin();
        spring_iterator!=(static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))->SpringsEnd();
        ++spring_iterator)
    {
        unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
        unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
        
        c_vector<double, DIM> force = CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index, rCellPopulation); 
        
        rForces[nodeB_global_index] -= force;
        rForces[nodeA_global_index] += force;
    }
}


template<unsigned DIM>
c_vector<double, DIM> CryptModelInteractionForce<DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
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

    // Get the unit vector parallel to the line joining the two nodes
    c_vector<double, DIM> unit_difference;
    
    if (rCellPopulation.IsMeshBasedCellPopulation())
    {
        /*
         * We use the mesh method GetVectorFromAtoB() to compute the direction of the
         * unit vector along the line joining the two nodes, rather than simply subtract
         * their positions, because this method can be overloaded (e.g. to enforce a
         * periodic boundary in Cylindrical2dMesh).
         */
        unit_difference = (static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))->rGetMesh().GetVectorFromAtoB(node_a_location, node_b_location);
    }
    else
    {
        unit_difference = node_b_location - node_a_location;
    }
   
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
void CryptModelInteractionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
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

template class CryptModelInteractionForce<1>;
template class CryptModelInteractionForce<2>;
template class CryptModelInteractionForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CryptModelInteractionForce)

