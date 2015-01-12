#include "BasementMembraneForce.hpp"
#include "Debug.hpp"

/*
 * Functions to allow the vectors of coordinates to be compared in
 * terms of their x-value and y-value accordingly.
 */
bool SortAccordingToXCoordinate(const c_vector<double, 2> lhs, const c_vector<double, 2> rhs)
{
    return lhs[0] < rhs[0];
}

/**
 * To avoid warnings on some compilers, C++ style initialization of member
 * variables should be done in the order they are defined in the header file.
 */
BasementMembraneForce::BasementMembraneForce()
   :  AbstractForce<2>(),
   mBasementMembraneParameter(DOUBLE_UNSET),
   mCryptBaseCurvature(DOUBLE_UNSET),
   mLeftBoundary(DOUBLE_UNSET),
   mRightBoundary(DOUBLE_UNSET),
   mUsePositionDependentMembraneForce(false),
   mMembraneForceMultiplier(DOUBLE_UNSET)
{
    // Sets up output file
//	OutputFileHandler output_file_handler("CurvatureData/", false);
//	mMeinekeOutputFile = output_file_handler.OpenOutputFile("results.curvature");
}

BasementMembraneForce::~BasementMembraneForce()
{
//    mMeinekeOutputFile->close();
}

void BasementMembraneForce::SetBasementMembraneParameter(double basementMembraneParameter)
{
	mBasementMembraneParameter = basementMembraneParameter;
}


double BasementMembraneForce::GetBasementMembraneParameter()
{
	return mBasementMembraneParameter;
}


void BasementMembraneForce::SetCryptBaseCurvature(double cryptBaseCurvature, double leftBoundary, double rightBoundary)
{
	mCryptBaseCurvature = cryptBaseCurvature;
	mLeftBoundary = leftBoundary;
	mRightBoundary = rightBoundary;
}


double BasementMembraneForce::GetCryptBaseCurvature()
{
	return mCryptBaseCurvature;
}


void BasementMembraneForce::SetPositionDependentMultiplier(bool usePositionDependentMembraneForce, double membraneForceMultiplier)
{
	mUsePositionDependentMembraneForce = usePositionDependentMembraneForce;
	mMembraneForceMultiplier = membraneForceMultiplier;
}


double BasementMembraneForce::GetPositionDependentMultiplier()
{
	return mMembraneForceMultiplier;
}


void BasementMembraneForce::RemoveDuplicates1D(std::vector<unsigned>& rVectorWithDuplicates)
{
    std::sort(rVectorWithDuplicates.begin(), rVectorWithDuplicates.end());
    rVectorWithDuplicates.erase(std::unique(rVectorWithDuplicates.begin(), rVectorWithDuplicates.end()), rVectorWithDuplicates.end());
}

/*
 * A method to find all the pairs of connections between healthy epithelial cells and labelled tissue cells.
 * Returns a vector of node pairings, without repeats. The first of each pair is the epithelial node index,
 * and the second is the tissue node index. Updating so that it also returns mutant-labelled cell pairs.
 */

std::vector<c_vector<unsigned, 2> > BasementMembraneForce::GetEpithelialTissuePairs(AbstractCellPopulation<2>& rCellPopulation)
{
    MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

    // Create a vector to record the pairs of nodes corresponding to *joined* epithelial and tissue nodes
    std::vector<c_vector<unsigned, 2> > node_pairs;
    c_vector<double, 2> pair;

    // We iterate over all cells in the tissue, and deal only with those that are epithelial cells
    for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
    	boost::shared_ptr<AbstractCellMutationState> p_state = cell_iter->GetMutationState();

    	// Need these to not be stromal cells (and not dead)
    	if ( (p_state->IsType<StromalCellMutationState>()==false) && (!cell_iter->IsDead()) )	// an epithelial cell
    	{
    		Node<2>* p_node = p_tissue->GetNodeCorrespondingToCell(*cell_iter);	// Pointer to node
    		unsigned node_index = p_node->GetIndex();

    		assert(!(p_tissue->IsGhostNode(node_index)));  // bit unnecessary at this stage but paranoia demands it

    		// ITERATE OVER CONTAINING ELEMENTS and only work with those that DO NOT contain ghost nodes

    		std::vector<unsigned> tissue_nodes;

    		for (Node<2>::ContainingElementIterator iter = p_node->ContainingElementsBegin();
		         iter != p_node->ContainingElementsEnd();
		         ++iter)
    		{
    			bool element_contains_ghost_nodes = false;

    			// Get a pointer to the element
    			Element<2,2>* p_element = p_tissue->rGetMesh().GetElement(*iter);

    			// ITERATE OVER NODES owned by this element
    			for (unsigned local_index=0; local_index<3; local_index++)
    			{
    				unsigned nodeBGlobalIndex = p_element->GetNodeGlobalIndex(local_index);

    				if (p_tissue->IsGhostNode(nodeBGlobalIndex) == true)
    				{
    					element_contains_ghost_nodes = true;
    					break; 				// This should break out of the inner for loop
    				}
    			}

				if (element_contains_ghost_nodes==false)
				{
                    // ITERATE OVER NODES owned by this element
                    for (unsigned local_index=0; local_index<3; local_index++)
                    {
                        unsigned nodeBGlobalIndex = p_element->GetNodeGlobalIndex(local_index);

                        CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

						if (p_cell->GetMutationState()->IsType<StromalCellMutationState>()==true)
						{
							// Store the index of each tissue node that is attached to the epithelial
							// node. There will be repetitions due to iterating over neighbouring elements
							tissue_nodes.push_back(nodeBGlobalIndex);
						}
                    }
				}
			}

    		// Remove any nodes that have been found twice
    		RemoveDuplicates1D(tissue_nodes);

    		// Now construct the vector of node pairs
    		for (unsigned i=0; i<tissue_nodes.size(); i++)
    		{
    			pair[0] = node_index;
    			pair[1] = tissue_nodes[i];
    			node_pairs.push_back(pair);

    			// Check that these node share a common element
				bool has_common_element = false;

				// The elements that contain this epithelial node:
				std::set<unsigned> epithelial_elements = rCellPopulation.GetNode(node_index)->rGetContainingElementIndices();
				assert(epithelial_elements.size() != 0);

				// The elements that contain the tissue node:
				std::set<unsigned> tissue_elements = rCellPopulation.GetNode(tissue_nodes[i])->rGetContainingElementIndices();
				assert(tissue_elements.size() != 0);

				// Loop over all elements that contain the tissue node
				for (Node<2>::ContainingElementIterator elt_it = rCellPopulation.GetNode(tissue_nodes[i])->ContainingElementsBegin();
				         elt_it != rCellPopulation.GetNode(tissue_nodes[i])->ContainingElementsEnd();
				         ++elt_it)
				{
					unsigned elt_index = *elt_it;

					bool elt_contains_ghost_nodes = DoesElementContainGhostNodes(rCellPopulation, elt_index);

					// Keep only those elements that also contain the epithelial node, but do not have ghost nodes
					if ( (elt_contains_ghost_nodes == false) && (epithelial_elements.find(elt_index) != epithelial_elements.end()) )
					{
						// Common element
						has_common_element = true;
						break;
					}
				}

				if (!has_common_element)
				{
					TRACE("No common element between:");
					PRINT_2_VARIABLES(node_index,tissue_nodes[i]);
				}
				assert(has_common_element);
    		}
    	}
    }

	return node_pairs;
}

/*
 * Method to determine whether an element contains ghost nodes
 */

bool BasementMembraneForce::DoesElementContainGhostNodes(AbstractCellPopulation<2>& rCellPopulation, unsigned elementIndex)
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

	bool element_contains_ghost_nodes = false;

	// Get a pointer to the element
	Element<2,2>* p_element = p_tissue->rGetMesh().GetElement(elementIndex);

	// ITERATE OVER NODES owned by this element
	for (unsigned local_index=0; local_index<3; local_index++)
	{
		if (p_tissue->IsGhostNode(p_element->GetNodeGlobalIndex(local_index)) == true)
		{
			element_contains_ghost_nodes = true;

		}
	}

	return element_contains_ghost_nodes;
}

/* Method to return the current coordinates of the crypt orifice and
 * crypt base - these can be used to accurately define the region of the
 * crypt base. (This will be the y-coordinate in 2D, or the z coordinate in 3D)
 * [0] - y or z-coordinate of orifice
 * [1] - y or z-coordinate of base
 */

c_vector<double,2> BasementMembraneForce::GetCryptHeightExtremes(AbstractCellPopulation<2>& rCellPopulation)
{
    MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

    // Create a vector to store the y-coordinates of the lowest point of the crypt base and the highest point of the
    // crypt orifice
    c_vector<double,2> height_extremes;

    double max_height = 0.0;
    double min_height = DBL_MAX;

    double current_height_coordinate;

    // We iterate over all cells in the tissue, and deal only with those that are epithelial cells
    for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
    	boost::shared_ptr<AbstractCellMutationState> p_state = cell_iter->GetMutationState();

	   	// Need these to not be labelled cells
	   	if ( (p_state->IsType<StromalCellMutationState>()==false) )
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

/*
 * Using the vector of node pairs found in GetEpithelialTissuePairs to determine the curvature
 * of the curve that passes through the midpoints of the neighbouring springs, given a epithelial-tissue
 * node pairing. (Note - it ignores the end pairs because one of the common elements will contain ghost nodes, but this
 * should only crop up if you don't have periodic boundary conditions)
 * Updating this so that it will still find the curvature if one of the epithelial cells is a mutant cell, eg. apc2 hit
 */

double BasementMembraneForce::GetCurvatureFromNodePair(AbstractCellPopulation<2>& rCellPopulation, unsigned epithelialNodeIndex,
																unsigned tissueNodeIndex)
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

    // Seeking the common elements that contain both the epithelial node and the tissue node
	// Note: Don't want to count any elements that have ghost nodes

	std::vector<unsigned> common_elements;	// Initialising

	// The elements that contain this epithelial node:
    std::set<unsigned> epithelial_elements = rCellPopulation.GetNode(epithelialNodeIndex)->rGetContainingElementIndices();
    assert(epithelial_elements.size() != 0);

    assert(rCellPopulation.GetNode(tissueNodeIndex)->GetNumContainingElements() != 0);

    // Loop over all elements that contain the tissue node
    for (Node<2>::ContainingElementIterator elt_it = rCellPopulation.GetNode(tissueNodeIndex)->ContainingElementsBegin();
         elt_it != rCellPopulation.GetNode(tissueNodeIndex)->ContainingElementsEnd();
         ++elt_it)
    {
    	unsigned elt_index = *elt_it;

    	bool elt_contains_ghost_nodes = DoesElementContainGhostNodes(rCellPopulation, elt_index);

    	// Keep only those elements that also contain the epithelial node, but do not have ghost nodes

        if ( (elt_contains_ghost_nodes == false) && (epithelial_elements.find(elt_index) != epithelial_elements.end()) )
        {
        	// Common element
           	common_elements.push_back(elt_index);
        }
    }

	assert(common_elements.size() != 0);		// This is bad - the nodes should be connected in the first place...

	// We iterate over these common elements to find the midpoints of the springs that
	// connect epithelial and  tissue nodes

	c_vector<double, 2> spring_midpoint_a, spring_midpoint_b, spring_midpoint_c;

	spring_midpoint_b = p_tissue->GetNode(epithelialNodeIndex)->rGetLocation() + 0.5*(p_tissue->rGetMesh().GetVectorFromAtoB(p_tissue->GetNode(epithelialNodeIndex)->rGetLocation(), p_tissue->GetNode(tissueNodeIndex)->rGetLocation()));

    // If there is only one such common element, then this epithelial node will be at either end and so we don't
    // consider the force along the very first / very last spring (only happens if you don't use a cylindrical
	// mesh!)
    if (common_elements.size() == 1)
    {
//    	std::cout << "Had a situation with a single common element. Check that you are using a periodic, cylindrical mesh. \n" << std::flush;
    	double curvature = 0.0;
    	return curvature;
    }

    else
    {
    	assert(common_elements.size() == 2);		// Should only be two common elements

    	for (std::vector<unsigned>::iterator elem_iter = common_elements.begin();
    									   	 elem_iter != common_elements.end();
    									   	 ++elem_iter)
    	{
    		// Need to determine the cell type of each local node in the element
    		// Want to find the midpoint between epithelial and tissue pairs

    		// Looking at the three nodes which form the vertices of the element
    		unsigned global_index_A = p_tissue->rGetMesh().GetElement(*elem_iter)->GetNodeGlobalIndex(0);
	   		unsigned global_index_B = p_tissue->rGetMesh().GetElement(*elem_iter)->GetNodeGlobalIndex(1);
	   		unsigned global_index_C = p_tissue->rGetMesh().GetElement(*elem_iter)->GetNodeGlobalIndex(2);

	   		unsigned E=UINT_MAX;  // E - the epithelial node we're interested in,
	   		unsigned T=UINT_MAX;  // T - the tissue node it's connected to, and
	   		unsigned P=UINT_MAX;  // P - the other node in the element.

	   	    assert(!(p_tissue->IsGhostNode(global_index_A)));
	   	    assert(!(p_tissue->IsGhostNode(global_index_B)));
	   	    assert(!(p_tissue->IsGhostNode(global_index_C)));

	   		if (global_index_A == epithelialNodeIndex)
	   		{
	   			E = global_index_A;

	   			if (global_index_B == tissueNodeIndex)
	   			{
	   				T = global_index_B;
	   				P = global_index_C;
	   			}
	   			else
	   			{
	   				P = global_index_B;
	   				T = global_index_C;
	   			}
	   		}
	   		else if (global_index_B == epithelialNodeIndex)
	   		{
	   			E = global_index_B;

	   			if (global_index_A == tissueNodeIndex)
	   			{
	   				T = global_index_A;
	   				P = global_index_C;
	   			}
	   			else
	   			{
	   				P = global_index_A;
	   				T = global_index_C;
	   			}
	   		}
	   		else if (global_index_C == epithelialNodeIndex)
	   		{
	   			E = global_index_C;

	   			if (global_index_A == tissueNodeIndex)
	   			{
	   				T = global_index_A;
	   				P = global_index_B;
	   			}
	   			else
	   			{
	   				P = global_index_A;
	   				T = global_index_B;
	   			}
	   		}

	   		assert(E<UINT_MAX);
	   		assert(T<UINT_MAX);
	   		assert(P<UINT_MAX);

	   		// Check that we have assigned the epithelial node correctly
	   		boost::shared_ptr<AbstractCellMutationState> p_E_state = p_tissue->GetCellUsingLocationIndex(E)->GetMutationState();
	   		assert( (p_E_state->IsType<WildTypeCellMutationState>()) || (p_E_state->IsType<ApcTwoHitCellMutationState>()) );
	   		assert(!p_E_state->IsType<StromalCellMutationState>());
			assert(E == epithelialNodeIndex);

	   		// Check that we have assigned the tissue node correctly
	   		CellPtr p_T_cell = rCellPopulation.GetCellUsingLocationIndex(T);
	   		boost::shared_ptr<AbstractCellMutationState> p_T_state = p_tissue->GetCellUsingLocationIndex(T)->GetMutationState();
	   		assert(p_T_state->IsType<StromalCellMutationState>());
	   		assert(T == tissueNodeIndex);

	   		/*
	   		 * Now we work with E (the epithelial node), T (the tissue node it's connected to) and P, the other node,
	   		 * which we will now assign as either P1 (det < 0) or P2 (det > 0).
	   		 *
	   		 * We also need to determine whether P1 and P2 are epithelial or tissue, as this will affect which vector we
	   		 * choose to take the spring midpoint from.
	   		 */
	   		c_vector<double, 2> vector_T_to_E = p_tissue->rGetMesh().GetVectorFromAtoB(p_tissue->GetNode(T)->rGetLocation(),p_tissue->GetNode(E)->rGetLocation());
	   		c_vector<double, 2> vector_E_to_T = p_tissue->rGetMesh().GetVectorFromAtoB(p_tissue->GetNode(E)->rGetLocation(),p_tissue->GetNode(T)->rGetLocation());
	   		c_vector<double, 2> vector_E_to_P = p_tissue->rGetMesh().GetVectorFromAtoB(p_tissue->GetNode(E)->rGetLocation(),p_tissue->GetNode(P)->rGetLocation());
	   		c_vector<double, 2> vector_T_to_P = vector_T_to_E + vector_E_to_P;
	   		c_vector<double, 2> vector_P_to_T = p_tissue->rGetMesh().GetVectorFromAtoB(p_tissue->GetNode(P)->rGetLocation(),p_tissue->GetNode(T)->rGetLocation());

	   		// Now calculate the determinant (TE x EP)

	   		double det = vector_T_to_E[0]*vector_E_to_P[1] - vector_T_to_E[1]*vector_E_to_P[0];
	  		assert(!isnan(det));
	   		assert(det != 0.0);

	   		/*
	   		 * If det < 0 then P = P1 and we can assign spring_midpoint_c
	   		 * If det > 0 then P = P2 and we can assign spring_midpoint_a
	   		 * Also need to take into account whether P is a epithelial or tissue node
	   		 * to choose the right spring
	   		 */
	   		boost::shared_ptr<AbstractCellMutationState> p_state = p_tissue->GetCellUsingLocationIndex(P)->GetMutationState();
	   		CellPtr p_cell_P = p_tissue->GetCellUsingLocationIndex(P);

	   		if ( (det < 0) && (p_state->IsType<StromalCellMutationState>()==false) )	// P = Epithelial, not labelled
	   		{
	   			spring_midpoint_c = p_tissue->GetNode(E)->rGetLocation() + vector_E_to_P + 0.5*vector_P_to_T;
	   		}
	   		else if ((det < 0) && (p_state->IsType<StromalCellMutationState>()==true))	// P = Tissue
	   		{
	   			spring_midpoint_c = p_tissue->GetNode(E)->rGetLocation() + 0.5*vector_E_to_P;
	   		}

	   		else if ( (det > 0) && (p_state->IsType<StromalCellMutationState>()==false) )	// P = Epithelial, not labelled
	   		{
	   			spring_midpoint_a = p_tissue->GetNode(E)->rGetLocation() + vector_E_to_T + 0.5*vector_T_to_P;
	   		}

	   		else if ((det > 0) && (p_state->IsType<StromalCellMutationState>()==true))	// P = Tissue
	   		{
	   			spring_midpoint_a = p_tissue->GetNode(E)->rGetLocation() + 0.5*vector_E_to_P;
	   		}
    	}

    	double curvature = FindParametricCurvature(spring_midpoint_a, spring_midpoint_b, spring_midpoint_c);

    	// Need to subtract the target curvature if the epithelial node lies in the crypt base
    	c_vector<double, 2> epithelial_location = p_tissue->GetNode(epithelialNodeIndex)->rGetLocation();

		// We need to use the current crypt height and base levels to determine where to apply the non-zero curvature
		// [0] - y-coordinate of orifice
		// [1] - y-coordinate of base
		c_vector<double, 2> height_extremes = GetCryptHeightExtremes(rCellPopulation);

		double top_of_crypt_base = height_extremes(1) + (height_extremes(0) - height_extremes(1))*0.2;

//		// For cross section crypt model:
		if ( (epithelial_location[1] <= top_of_crypt_base))
		{
			curvature -= mCryptBaseCurvature;
		}

		// For sandwich box model:
//		if ( (epithelial_location[1] <= 20.0) && (epithelial_location[0] > mLeftBoundary) && (epithelial_location[0] < mRightBoundary) )
//		{
//			curvature -= mCryptBaseCurvature;
//		}

    	assert(!isnan(curvature));
    	return curvature;
    }
}

/*
 * Function to return the curvature between three midpoints parametrically - in this case, we find the normal
 * to the vector joining the left and right midpoints, and then find the perpendicular distance of the centre midpoint
 * from the left->right vector
 */

double BasementMembraneForce::GetCurvatureFromMidpoints(AbstractCellPopulation<2>& rCellPopulation,
																c_vector<double, 2> leftMidpoint,
																c_vector<double, 2> centreMidpoint,
																c_vector<double, 2> rightMidpoint)
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

	// Firstly find the normal to the vector joining the left and right midpoints
	c_vector<double, 2>	vector_left_to_right = p_tissue->rGetMesh().GetVectorFromAtoB(leftMidpoint, rightMidpoint);
	c_vector<double, 2> normal_vector;
	normal_vector[0] = vector_left_to_right[1];
	normal_vector[1] = -vector_left_to_right[0];

	c_vector<double, 2> vector_left_to_centre = p_tissue->rGetMesh().GetVectorFromAtoB(leftMidpoint, centreMidpoint);

	double curvature = normal_vector[0]*vector_left_to_centre[0] + normal_vector[1]*vector_left_to_centre[1];

	return curvature;
}

/*
859	 * Function to return the curvature between three points parametrically - the midpoints of the springs connecting the
860	 * transit cells to the differentiated cells. NB. The input arguments need to be in order from either left to right
861	 * or right to left. If they are wrongly arranged (eg. middle, left, right) then you get a different curvature,
862	 * but left->right = -(right-> left).
863	 */

double BasementMembraneForce::FindParametricCurvature(c_vector<double, 2> leftMidpoint,
															c_vector<double, 2> centreMidpoint,
															c_vector<double, 2> rightMidpoint)
{
	// Firstly find the parametric intervals
	double left_s = sqrt(pow(centreMidpoint[0] - leftMidpoint[0],2) + pow(centreMidpoint[1] - leftMidpoint[1],2));
	double right_s = sqrt(pow(rightMidpoint[0] - centreMidpoint[0],2) + pow(rightMidpoint[1] - centreMidpoint[1],2));

	assert(left_s >= 0);
	assert(right_s >= 0);

	double sum_intervals = left_s + right_s;

	double x_prime = (rightMidpoint[0] - leftMidpoint[0])/sum_intervals;
	double y_prime = (rightMidpoint[1] - leftMidpoint[1])/sum_intervals;

	double x_double_prime = 2*(left_s*rightMidpoint[0] - sum_intervals*centreMidpoint[0] + right_s*leftMidpoint[0])/(left_s*right_s*sum_intervals);
	double y_double_prime = 2*(left_s*rightMidpoint[1] - sum_intervals*centreMidpoint[1] + right_s*leftMidpoint[1])/(left_s*right_s*sum_intervals);

	double curvature = (x_prime*y_double_prime - y_prime*x_double_prime)/pow((pow(x_prime,2) + pow(y_prime,2)),3/2);

	return curvature;
}


/*
 * A method to return the number of elements that contain a particular node,
 * excluding those elements that have ghost nodes
 */

unsigned BasementMembraneForce::GetNumContainingElementsWithoutGhostNodes(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex)
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

    // Get pointer to the node
    Node<2>* p_node = p_tissue->GetNode(nodeIndex);
    assert(!(p_tissue->IsGhostNode(nodeIndex)));

    unsigned num_elements_with_no_ghost_nodes = 0;		// Initialise

    // Iterate over containing elements
    for (Node<2>::ContainingElementIterator iter = p_node->ContainingElementsBegin();
         iter != p_node->ContainingElementsEnd(); ++iter)
    {
        bool element_contains_ghost_nodes = false;
        Element<2,2>* p_element = p_tissue->rGetMesh().GetElement(*iter);

        // Iterate over nodes owned by this element
        for (unsigned local_index=0; local_index<3; local_index++)
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

std::set<unsigned> BasementMembraneForce::GetNeighbouringNodeIndices(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex)
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

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

bool BasementMembraneForce::HasEpithelialCellDetachedFromBasementMembrane(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex)
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

	bool has_cell_detached = false;	// Initialising

   	std::set<unsigned> neighbours = GetNeighbouringNodeIndices(rCellPopulation, nodeIndex);

   	unsigned num_stromal_neighbours = 0;

   	// Iterate over the neighbouring cells to check the number of tissue cell neighbours

   	for(std::set<unsigned>::iterator neighbour_iter=neighbours.begin();
   							neighbour_iter != neighbours.end();
   							++neighbour_iter)
	{
   		boost::shared_ptr<AbstractCellMutationState> p_state = p_tissue->GetCellUsingLocationIndex(*neighbour_iter)->GetMutationState();
		if ( (!p_tissue->IsGhostNode(*neighbour_iter)) && (p_state->IsType<StromalCellMutationState>()==true))
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

void BasementMembraneForce::AddForceContribution(std::vector<c_vector<double, 2> >& rForces,
                                                         AbstractCellPopulation<2>& rCellPopulation)
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

	// First determine the force acting on each epithelial cell due to the basement membrane
	// Start by identifying the epithelial-tissue node pairs (now also returns any apc2hit-tissue pairs)
	std::vector<c_vector<unsigned, 2> > node_pairs = GetEpithelialTissuePairs(rCellPopulation);

	// We loop over the epithelial-tissue node pairs to find the force acting on that
	// epithelial node, and the direction in which it acts
	for (unsigned i=0; i<node_pairs.size(); i++)
	{
		unsigned epithelial_node_index = node_pairs[i][0];
		unsigned tissue_node_index = node_pairs[i][1];

   		CellPtr p_cell_epithelial = p_tissue->GetCellUsingLocationIndex(epithelial_node_index);
   		assert(p_cell_epithelial->GetMutationState()->IsType<StromalCellMutationState>() == false);

		CellPtr p_cell_tissue = p_tissue->GetCellUsingLocationIndex(tissue_node_index);
		assert(p_cell_tissue->GetMutationState()->IsType<StromalCellMutationState>() == true);

		c_vector<double, 2> epithelial_location = rCellPopulation.GetNode(epithelial_node_index)->rGetLocation();
		c_vector<double, 2> tissue_location = rCellPopulation.GetNode(tissue_node_index)->rGetLocation();

		// The force due to the basal lamina acts along the spring connecting the epithelial and tissue nodes, T->E direction
		c_vector<double, 2> curvature_force_direction = p_tissue->rGetMesh().GetVectorFromAtoB(tissue_location, epithelial_location);

		double distance_between_nodes = norm_2(curvature_force_direction);
		assert(distance_between_nodes > 0);
		assert(!isnan(distance_between_nodes));

		curvature_force_direction /= distance_between_nodes;

		double curvature = GetCurvatureFromNodePair(rCellPopulation, epithelial_node_index, tissue_node_index);

		double basement_membrane_parameter = GetBasementMembraneParameter();

		// Here we can apply a different basement membrane parameter to a specific region of the crypt. Currently
		// this is applying a multiplier to those cells in the rounded crypt base (but mMembraneForceMultiplier defaults
		// to 1 unless told otherwise).
		if (mUsePositionDependentMembraneForce)
		{
			// If the basement membrane force is defined to be position dependent
			c_vector<double, 2> height_extremes = GetCryptHeightExtremes(rCellPopulation);
			double top_of_crypt_base = height_extremes(1) + (height_extremes(0) - height_extremes(1))*0.25;

			if  (epithelial_location[1] <= top_of_crypt_base)
			{
				basement_membrane_parameter *= mMembraneForceMultiplier;
			}
		}

		c_vector<double, 2> force_due_to_basement_membrane = basement_membrane_parameter*curvature*curvature_force_direction;

		// Add the force due to the basal lamina to the forces acting on that epithelial node
		rForces[epithelial_node_index] += force_due_to_basement_membrane;
	}

//	// Find the epithelial nodes and apply a migration force
//
//    for (MeshBasedCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
//         cell_iter != rCellPopulation.End(); ++cell_iter)
//    {
//		c_vector<double, 2> active_migration_direction;
//		double active_migration_parameter = 1.5;
//
//    	if(cell_iter->GetCellMutationState()->IsType<StromalCellMutationState>() == false)		// If epithelial
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
//		    	if( (!p_tissue->IsGhostNode(*neighbour_iter))
//					&& ((p_tissue->GetCellUsingLocationIndex(*neighbour_iter)->GetCellMutationState()->IsType<StromalCellMutationState>() == false)) )
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
//    for (MeshBasedCellPopulation<2>::SpringIterator spring_iterator=(static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation))->SpringsBegin();
//            spring_iterator!=(static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation))->SpringsEnd();
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

}

void BasementMembraneForce::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile <<  "\t\t\t<BasementMembraneParameter>"<<  mBasementMembraneParameter << "</BasementMembraneParameter> \n" ;
	*rParamsFile <<  "\t\t\t<CryptBaseCurvature>"<<  mCryptBaseCurvature << "</CryptBaseCurvature> \n" ;
	*rParamsFile <<  "\t\t\t<LeftBoundary>"<<  mLeftBoundary << "</LeftBoundary> \n" ;
	*rParamsFile <<  "\t\t\t<RightBoundary>"<<  mRightBoundary << "</RightBoundary> \n" ;
	*rParamsFile <<  "\t\t\t<UsePositionDependentMembraneForce>"<<  mUsePositionDependentMembraneForce << "</UsePositionDependentMembraneForce> \n" ;
	*rParamsFile <<  "\t\t\t<MembraneForceMultiplier>"<<  mMembraneForceMultiplier << "</MembraneForceMultiplier> \n" ;

	// Call direct parent class
	AbstractForce<2>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(BasementMembraneForce)
