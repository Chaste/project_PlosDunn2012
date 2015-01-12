#ifndef CROSSSECTIONMODEL_HPP_
#define CROSSSECTIONMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>

#include "AbstractMesh.hpp"
#include "CellBasedSimulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "AbstractVanLeeuwen2009WntSwatCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StromalCellMutationState.hpp"
#include "Cylindrical2dMesh.hpp"
#include "CryptModelInteractionForce.hpp"

class CrossSectionModel : public CellBasedSimulation<2>
{
    // Allow tests to access private members, in order to test computation of
    // private functions
    friend class TestCrossSectionModel;

protected:

    /** Define a stopping event which says stop if you have mismatched boundary elements (mMismatchedBoundaryElements == true)*/
    bool StoppingEventHasOccurred();

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the simulation and member variable.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & boost::serialization::base_object<CellBasedSimulation<2> >(*this);
        archive & mUseJiggledBottomCells;
        archive & mPinBottomCells;
        archive & mCellPopulationWidth;
        archive & mBasementMembraneParameter;
        archive & mMaxHeightForPinnedCells;
        archive & mRightHandBoundaryForPinnedCells;
        archive & mLeftHandBoundaryForPinnedCells;
    }

    // Whether to use a flat bottom surface or to jiggle the cells on the bottom surface
    bool mUseJiggledBottomCells;
    
    // The width of the original mesh of nodes (i.e. the width of the tissue box)
    double mCellPopulationWidth;
    
    // Whether to pin the bottom row of cells
    bool mPinBottomCells;

    // The file that the values of beta catenin is written out to
    out_stream mBetaCatResultsFile;

    // The file that the number of cells at each timestep are written out to
    out_stream mNumCellsResultsFile;

    // Helper member that is a static cast of the tissue
    MeshBasedCellPopulationWithGhostNodes<2>* mpStaticCastCellPopulation;

    // The parameter for the basal lamina force - multiplies the local curvature
    double mBasementMembraneParameter;

    // Parameter used to define the height below which cells should be pinned
    double mMaxHeightForPinnedCells;

    // Parameter used to define the height below which cells should be pinned
    double mRightHandBoundaryForPinnedCells;

    // Parameter used to define the height below which cells should be pinned
    double mLeftHandBoundaryForPinnedCells;

    // The output file directory for the simulation data - number of popped up cells,
    // average distance between cells etc.
    std::string mDataOutputFile;

    /**
     * Calculates the new locations of a dividing cell's cell centres.
     * There are two choices here: one returns a division vector that is parallel to 
     * the vector connecting the neighbouring epithelial cell centres, and the second 
     * chooses a random direction, as below.
     * 
     * Moves the dividing node a bit and returns co-ordinates for the new node.
     * It does this by picking a random direction (0->2PI) and placing the parent
     * and daughter in opposing directions on this axis.
     *
     * @param pParentCell pointer to the parent cell
     *
     * @return daughter_coords the coordinates for the daughter cell.
     */
    c_vector<double, 2> CalculateCellDivisionVector(CellPtr pParentCell);

    /**
     * Overridden WriteVisualizerSetupFile() method.
     *
     * Writes out special information about the mesh to the visualizer.
     */
    void WriteVisualizerSetupFile();

    /**
     * Get the data output directory of the simulation.
     */
     std::string GetDataOutputFile();

     /**
      * Overridden SetupSolve() method.
      *
      * Write results of interest to file if required.
      */
     void SetupSolve();

     /**
      * Overridden PostSolve() method.
      *
      * Write results of interest to file if required.
      */
     void PostSolve();


    /**
     * Overridden AfterSolve() method.
     *
     * Closes beta catenin results file if required, then calls
     * the base class method.
     */
    void AfterSolve();

public :

    /**
     *  Constructor.
     *
     *  @param rCellPopulation A tissue facade class (contains a mesh and cells)
     *  @param deleteCellPopulationAndForceCollection Whether to delete the tissue and force collection on destruction to free up memory
     *  @param initialiseCells whether to initialise cells (set to false when loading from an archive)
     *  @param tissueWidth Width of mesh of real nodes
     *  @param sloughTransitCells whether to slough off those transit cells that move beyond the boundaries
     */
	CrossSectionModel(AbstractCellPopulation<2>& rCellPopulation,
                      bool deleteCellPopulationAndForceCollection=false,
                      bool initialiseCells=true,
                      double mCellPopulationWidth = 0.1,
                      bool mPinBottomCells = true);

    /** Set method for mUseJiggledBottomCells. */
    void UseJiggledBottomCells();

    /**
     * Overridden ApplyCellPopulationBoundaryConditions() method.
     *
     * Cells lining the bottom of the crypt are pinned. Cells cannot move beyond the initial
     * boundaries at zero and mCellPopulationWidth
     *
     * @param rOldLocations: the node locations at the previous time step
     */
    void ApplyCellPopulationBoundaryConditions(const std::vector<c_vector<double,2> >& rOldLocations);

    /**
     * Given a node, find a set containing the indices of its neighbouring nodes.
     *
     * @param nodeIndex global index of the node
     * @return its neighbouring node indices
     */
    std::set<unsigned> GetNeighbouringNodeIndices(unsigned nodeIndex);

    void SetMaxHeightForPinnedCells(double maxHeightForPinnedCells);

    void SetRightHandBoundaryForPinnedCells(double rightHandBoundaryForPinnedCells);

    void SetLeftHandBoundaryForPinnedCells(double leftHandBoundaryForPinnedCells);

    double GetMaxHeightForPinnedCells();
    
    c_vector<double,2> OutputVectorOfBasalLaminaResults(double basalLaminaParameter);
    
//    void Solve();

    std::vector<c_vector<double,4> > OutputVectorOfForceResults();

};


// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CrossSectionModel)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a BoxModelSimulation
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const CrossSectionModel * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<2> * p_tissue = &(t->rGetCellPopulation());
    ar & p_tissue;
}

/**
 * De-serialize constructor parameters and initialise a CrossSectionModel
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, CrossSectionModel * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<2>* p_tissue;
    ar >> p_tissue;

    // Invoke inplace constructor to initialise instance
    ::new(t)CrossSectionModel(*p_tissue, true, false, 0.1, true);
}
}
} // namespace

#endif /*CROSSSECTIONMODEL_HPP_*/
