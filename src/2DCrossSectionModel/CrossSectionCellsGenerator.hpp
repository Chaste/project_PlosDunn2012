/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#ifndef CROSSSECTIONCELLSGENERATOR_HPP_
#define CROSSSECTIONCELLSGENERATOR_HPP_

#include <boost/mpl/integral_c.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/if.hpp>

#include "CellsGenerator.hpp"
#include "StromalCellMutationState.hpp"
#include "TetrahedralMesh.hpp"
#include "CellLabel.hpp"

///**
// *  Small helper method, which returns whether the two classes given as the template parameters
// *  are identical or not
// */
//template<class T1, class T2>
//bool ClassesAreSame()
//{
//    using namespace boost::mpl;
//    using namespace boost;
//    typedef typename if_< is_same<T1, T2>, integral_c<unsigned, 1>, integral_c<unsigned, 0> >::type selector_t;
//    return  (selector_t()==1);
//}


/**
 * A subclass of CellsGenerator that generates cells for 2D cross-section crypt simulations.
 *
 * It is templated over types of cell cycle model.
 */
template<class CELL_CYCLE_MODEL>
class CrossSectionCellsGenerator : public CellsGenerator<CELL_CYCLE_MODEL,2>
{
public:
    /**
     * Whether cells are able to fully differentiate.
     * Defaults to false unless overridden.
     */
    bool CanCellsTerminallyDifferentiate();

    /**
     * Generates cells of a specified cell cycle type under the correct
     * crypt conditions: we define a layer of transit cells around the crypt 
     * and on the top layer
     *
     * @param rCells  An empty cells vector for this function to fill up
     * @param rMesh  The crypt mesh (can be cylindrical)
     * @param cryptEdgeRight  The right hand boundary of the crypt
     * @param cryptEdgeLeft   The left hand boundary of the crypt
     * @param height   The height of the tissue box that contains the crypt 	
     * @param realNodeIndices The node indices corresponding to real cells
     * @param ghostNodeIndices   The node indices corresponding to ghost nodes (inside and around the crypt)
     * @param initialiseCells  Whether to initialise the cell cycle models as each
     *   cell is created
     */
    void Generate(std::vector<CellPtr>& rCells,
			 	MutableMesh<2,2>& rMesh,
			 	double cryptEdgeRight,
			 	double cryptEdgeLeft,
			 	double height,
			 	std::vector<unsigned> realNodeIndices,		
			 	std::vector<unsigned> ghostNodeIndices,
			 	bool initialiseCells);

    /** Similar generator to the above, except this is used for experiment to test the Kaur and Potten 1986 results,
     * where cell division is eliminated after a certain time. In this situation, the cell cycle model will be the
     * FixedDurationGenerationBasedCellCycleModel, and the number of generations set such that cell division will cease only
     * once the crypt has reached equilibrium. The proliferative and non-proliferative regions will have to be defined to
     * begin with. Will have to play around with this to account for the shrinking that will occur from the inital approximate
     * geometry.
     */
    void GenerateForElimatingDivisionExperiments(std::vector<CellPtr>& rCells,
			 	MutableMesh<2,2>& rMesh,
			 	double cryptEdgeRight,
			 	double cryptEdgeLeft,
			 	double height,
			 	std::vector<unsigned> realNodeIndices,
			 	std::vector<unsigned> ghostNodeIndices,
			 	bool initialiseCells);
};


template<class CELL_CYCLE_MODEL>
void CrossSectionCellsGenerator<CELL_CYCLE_MODEL>::Generate(std::vector<CellPtr>& rCells,
	 														  MutableMesh<2,2>& rMesh,
	 														  double cryptEdgeRight,
	 														  double cryptEdgeLeft,
	 														  double height,
	 														  std::vector<unsigned> realNodeIndices,		
	 														  std::vector<unsigned> ghostNodeIndices,
	 														  bool initialiseCells)
{
    if(initialiseCells)
    {
	    for (std::vector<unsigned>::iterator real_node_iter=realNodeIndices.begin();
	    real_node_iter != realNodeIndices.end();
	     ++real_node_iter)
	    {
	        // Create a set of neighbouring node indices
	    	std::set<unsigned> neighbours;
	        unsigned total_neighbouring_ghosts = 0;

	        // Find the indices of the elements (triangles!) owned by this node
	    	std::set<unsigned> containing_elem_indices = rMesh.GetNode(*real_node_iter)->rGetContainingElementIndices();

	        // Iterate over these elements
	        for (std::set<unsigned>::iterator elem_iter=containing_elem_indices.begin();
	             elem_iter != containing_elem_indices.end();
	             ++elem_iter)
	        {
	            // Get all the nodes contained in this element
	            // (note that we've also included the original node, but this doesn't matter too much)

	        	unsigned neighbour_global_index;

	            for (unsigned local_index=0; local_index<3; local_index++)
	            {
	                neighbour_global_index = rMesh.GetElement(*elem_iter)->GetNodeGlobalIndex(local_index);
	                neighbours.insert(neighbour_global_index);
	            }
	        }

			for(std::set<unsigned>::iterator neighbour_iter=neighbours.begin();
			         neighbour_iter != neighbours.end();
				     ++neighbour_iter)
			{
		        for (std::vector<unsigned>::iterator ghost_node_iter=ghostNodeIndices.begin();
		        ghost_node_iter != ghostNodeIndices.end();
			     ++ghost_node_iter)
		        {
		        	if(*neighbour_iter == *ghost_node_iter)
		        	{
		        		total_neighbouring_ghosts++;
		        	}
		        }
			}

			boost::shared_ptr<AbstractCellMutationState> p_epithelial_state(new WildTypeCellMutationState);
			boost::shared_ptr<AbstractCellMutationState> p_stromal_state(new StromalCellMutationState);

			// Any cells that are attached to ghost nodes and line the top and inner layer of the crypt are proliferating cells
			if ( ( (total_neighbouring_ghosts > 1e-6) &&
					( ((rMesh.GetNode(*real_node_iter)->rGetLocation()[1] > 0.5) &&
					(rMesh.GetNode(*real_node_iter)->rGetLocation()[0] < cryptEdgeRight+1)) &&
					(rMesh.GetNode(*real_node_iter)->rGetLocation()[0] > cryptEdgeLeft-1) ) ) ||
					(rMesh.GetNode(*real_node_iter)->rGetLocation()[1] > height-0.5) )
			{
	            CELL_CYCLE_MODEL* p_model = new CELL_CYCLE_MODEL;
                p_model->SetCellProliferativeType(TRANSIT);
                p_model->SetDimension(2);
//                p_model->SetMaxTransitGenerations(UINT_MAX);
//                p_model->SetWntTransitThreshold(0.45);

				CellPropertyCollection collection;

	            double birth_time = -RandomNumberGenerator::Instance()->ranf()*
		       	                     (p_model->GetStemCellG1Duration()
		       	                      + p_model->GetSG2MDuration());

                // Want to track a couple of cells
//                if (*real_node_iter == 108)
//                {
//                    boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
//
//                	CellPropertyCollection collection;
//                	collection.AddProperty(p_label);
//                    CellPtr p_epithelial_cell(new Cell(p_epithelial_state, p_model, false, collection));
//					p_epithelial_cell->InitialiseCellCycleModel();
//	                p_epithelial_cell->SetBirthTime(birth_time);
//					rCells.push_back(p_epithelial_cell);
//
//                }
//                else
//                {
					CellPtr p_epithelial_cell(new Cell(p_epithelial_state, p_model, false, collection));
					p_epithelial_cell->InitialiseCellCycleModel();
					rCells.push_back(p_epithelial_cell);
					p_epithelial_cell->SetBirthTime(birth_time);
//		            p_model->SetMaxTransitGenerations(UINT_MAX);
//                }
           	}
	    	else
	    	{
	            CELL_CYCLE_MODEL* p_model = new CELL_CYCLE_MODEL;
                p_model->SetCellProliferativeType(DIFFERENTIATED);
                p_model->SetDimension(2);
//                p_model->SetWntTransitThreshold(0.45);

				CellPropertyCollection collection;

	            CellPtr p_stromal_cell(new Cell(p_stromal_state, p_model, false, collection));
	            p_stromal_cell->InitialiseCellCycleModel();

	            rCells.push_back(p_stromal_cell);
	    	}

	    	for (unsigned j=0; j<rCells.size(); j++)
	    	{

	    	}
	    }
    }
}

template<class CELL_CYCLE_MODEL>
void CrossSectionCellsGenerator<CELL_CYCLE_MODEL>::GenerateForElimatingDivisionExperiments(std::vector<CellPtr>& rCells,
	 														  MutableMesh<2,2>& rMesh,
	 														  double cryptEdgeRight,
	 														  double cryptEdgeLeft,
	 														  double height,
	 														  std::vector<unsigned> realNodeIndices,
	 														  std::vector<unsigned> ghostNodeIndices,
	 														  bool initialiseCells)
{
    if(initialiseCells)
    {
	    for (std::vector<unsigned>::iterator real_node_iter=realNodeIndices.begin();
	    real_node_iter != realNodeIndices.end();
	     ++real_node_iter)
	    {
	        // Create a set of neighbouring node indices
	    	std::set<unsigned> neighbours;
	        unsigned total_neighbouring_ghosts = 0;

	        // Find the indices of the elements (triangles!) owned by this node
	    	std::set<unsigned> containing_elem_indices = rMesh.GetNode(*real_node_iter)->rGetContainingElementIndices();

	        // Iterate over these elements
	        for (std::set<unsigned>::iterator elem_iter=containing_elem_indices.begin();
	             elem_iter != containing_elem_indices.end();
	             ++elem_iter)
	        {
	            // Get all the nodes contained in this element
	            // (note that we've also included the original node, but this doesn't matter too much)

	        	unsigned neighbour_global_index;

	            for (unsigned local_index=0; local_index<3; local_index++)
	            {
	                neighbour_global_index = rMesh.GetElement(*elem_iter)->GetNodeGlobalIndex(local_index);
	                neighbours.insert(neighbour_global_index);
	            }
	        }

			for(std::set<unsigned>::iterator neighbour_iter=neighbours.begin();
			         neighbour_iter != neighbours.end();
				     ++neighbour_iter)
			{
		        for (std::vector<unsigned>::iterator ghost_node_iter=ghostNodeIndices.begin();
		        ghost_node_iter != ghostNodeIndices.end();
			     ++ghost_node_iter)
		        {
		        	if(*neighbour_iter == *ghost_node_iter)
		        	{
		        		total_neighbouring_ghosts++;
		        	}
		        }
			}

			boost::shared_ptr<AbstractCellMutationState> p_epithelial_state(new WildTypeCellMutationState);
			boost::shared_ptr<AbstractCellMutationState> p_stromal_state(new StromalCellMutationState);

			// Any cells that are attached to ghost nodes and line the top and inner layer of the crypt are proliferating cells
			if ( ( (total_neighbouring_ghosts > 1e-6) &&
					( ((rMesh.GetNode(*real_node_iter)->rGetLocation()[1] > 0.5) &&
					(rMesh.GetNode(*real_node_iter)->rGetLocation()[0] < cryptEdgeRight+1)) &&
					(rMesh.GetNode(*real_node_iter)->rGetLocation()[0] > cryptEdgeLeft-1) ) ) ||
					(rMesh.GetNode(*real_node_iter)->rGetLocation()[1] > height-0.5) )
			{
	            CELL_CYCLE_MODEL* p_model = new CELL_CYCLE_MODEL;
                p_model->SetCellProliferativeType(TRANSIT);
                p_model->SetDimension(2);
//                p_model->SetMaxTransitGenerations(UINT_MAX);
//                p_model->SetWntTransitThreshold(0.45);

				CellPropertyCollection collection;

	            double birth_time = -RandomNumberGenerator::Instance()->ranf()*
		       	                     (p_model->GetStemCellG1Duration()
		       	                      + p_model->GetSG2MDuration());

                // Want to track a couple of cells
//                if (*real_node_iter == 108)
//                {
//                    boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
//
//                	CellPropertyCollection collection;
//                	collection.AddProperty(p_label);
//                    CellPtr p_epithelial_cell(new Cell(p_epithelial_state, p_model, false, collection));
//					p_epithelial_cell->InitialiseCellCycleModel();
//	                p_epithelial_cell->SetBirthTime(birth_time);
//					rCells.push_back(p_epithelial_cell);
//
//                }
//                else
//                {
					CellPtr p_epithelial_cell(new Cell(p_epithelial_state, p_model, false, collection));
					p_epithelial_cell->InitialiseCellCycleModel();
					rCells.push_back(p_epithelial_cell);
					p_epithelial_cell->SetBirthTime(birth_time);
//		            p_model->SetMaxTransitGenerations(UINT_MAX);
//                }
           	}
	    	else
	    	{
	            CELL_CYCLE_MODEL* p_model = new CELL_CYCLE_MODEL;
                p_model->SetCellProliferativeType(DIFFERENTIATED);
                p_model->SetDimension(2);
//                p_model->SetWntTransitThreshold(0.45);
                
				CellPropertyCollection collection;

	            CellPtr p_stromal_cell(new Cell(p_stromal_state, p_model, false, collection));
	            p_stromal_cell->InitialiseCellCycleModel();
	            
	            rCells.push_back(p_stromal_cell);
	    	}
			
	    	for (unsigned j=0; j<rCells.size(); j++)
	    	{

	    	}					
	    }
    }
}

#endif /* CROSSSECTIONCELLSGENERATOR_HPP_ */
