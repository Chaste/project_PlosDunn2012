#ifndef PERIODICCRYPTMODELINTERACTIONFORCE_HPP_
#define PERIODICCRYPTMODELINTERACTIONFORCE_HPP_

#include "LinearSpringWithVariableSpringConstantsForce.hpp"
#include "StromalCellMutationState.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

// Needed here to avoid serialization errors (on Boost<1.37)
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "AbstractVanLeeuwen2009WntSwatCellCycleModel.hpp"


/**
 * A subclass of GeneralisedLinearSpringForce with variable spring constants, that defines how the spring forces should
 * act on the 2D Cross Section Crypt Model. This is the periodic version of CryptModelInteractionForce, such that the cells
 * on the vertical edges experience forces due to image nodes created from the nodes on the opposite vertical edge. This
 * force class makes it possible to use periodic boundary conditions without imposing the cylindrical mesh. Hopefully.
 */

template<unsigned DIM>
class PeriodicCryptModelInteractionForce : public LinearSpringWithVariableSpringConstantsForce<DIM>
{
    friend class TestPeriodicCryptModelInteractionForce;

private :

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<GeneralisedLinearSpringForce<DIM> >(*this);
        archive & mUseCellTypeDependentSprings;
    	archive & mTransitTransitMultiplier;
    	archive & mDifferentiatedDifferentiatedMultiplier;
    	archive & mTransitDifferentiatedMultiplier;
    	archive & mUseEpithelialStromalCellDependentSprings;
    	archive & mEpithelialEpithelialMultiplier;
    	archive & mStromalStromalMultiplier;
    	archive & mEpithelialStromalMultiplier;
    	archive & mApcTwoHitStromalMultiplier;
        archive & mUseEdgeBasedSpringConstant;
        archive & mUseOneWaySprings;
        archive & mUsePositionDependentSpringConstants;
        archive & mSpringConstantsMultiplier;
    }

	// Whether to use different spring stiffness dependent on the types of neighbouring cells
	bool mUseCellTypeDependentSprings;

	// Multiplier for spring stiffness if cells are both transit type
	double mTransitTransitMultiplier;

	// Multiplier for spring stiffness if cells are both differentiated type
	double mDifferentiatedDifferentiatedMultiplier;

	// Multiplier for spring stiffness if cells are of different types
	double mTransitDifferentiatedMultiplier;

	// Whether to use different spring stiffness dependent on the mutation states of epithelial and stromal cells
	bool mUseEpithelialStromalCellDependentSprings;

	// Multiplier for spring stiffness if cells are both epithelial (wildtype for now)
	double mEpithelialEpithelialMultiplier;

	// Multiplier for spring stiffness if cells are both stromal
	double mStromalStromalMultiplier;

	// Multiplier for spring stiffness if one cell is epithelial and the other is stromal
	double mEpithelialStromalMultiplier;

	// Multiplier for spring stiffness if one cell is apc2hit and the other is stromal
	double mApcTwoHitStromalMultiplier;

    /** Whether to use spring constant proportional to cell-cell contact length/area (defaults to false) */
    bool mUseEdgeBasedSpringConstant;

    /** Use one way springs between epithelial and tissue nodes, that only act under repulsion */
    bool mUseOneWaySprings;

    /** Make the spring constants dependent on the position of a cell up the crypt */
    bool mUsePositionDependentSpringConstants;

    /** The multiplication factor for the spring constants */
    double mSpringConstantsMultiplier;

    // The output file directory for the simulation data that corresponds to the number of cells
    // killed by anoikis
    out_stream mMeinekeOutputFile;

//    std::string mOutputDirectory;

public :

    /**
     * Constructor.
     */
	PeriodicCryptModelInteractionForce();

    /**
     * Destructor.
     */
    ~PeriodicCryptModelInteractionForce();

//    void SetOutputDirectory(std::string outputDirectory);
//
//    std::string GetOutputDirectory();

    /* Set whether to use spring constants dependent on the types of neighbouring cells:
     * Transit-transit, transit-differentiated, differentiated-differentiated

     * @param useCellTypeDependentSprings  whether to use mutant springs
     * @param transitTransitMultiplier  the multiplier for springs connecting two transit cells
     * @param transitDifferentiatedMultiplier  the multiplier for springs connecting a transit cell with a differentiated cell
     * @param differentiatedDifferentiatedMultiplier the multiplier for springs connecting two differentiated cells
     */
    void SetCellTypeDependentSprings(bool useCellTypeDependentSprings, double transitTransitMultiplier = 1.0,
    		double differentiatedDifferentiatedMultiplier = 1.0,
    		double transitDifferentiatedMultiplier = 1.0);

    /* Set whether to use spring constants dependent on the mutation state of neighbouring cells:
     * Healthy-healthy, labelled-labelled and healthy-labelled
     * @param useCellLabelStateDependentSprings whether to use springs dependent on whether cells are labelled or not
     * @param healthyHealthyMultiplier the multiplier for springs connecting two healthy (epithelial) cells (not labelled)
     * @param healthyLabelledMultiplier the multiplier for springs connecting a healthy cell with a labelled cell
     * @param labelledLabelledMultiplier the multiplier for springs connecting two labelled cells
     * @param apcTwoHitLabelledMultiplier the multiplier for springs connecting an apc2hit cell and a labelled cell (This has been
     * hacked in here temporarily)
     */
    void SetEpithelialStromalCellDependentSprings(bool useEpithelialStromalCellDependentSprings = false,
    		double epithelialEpithelialMultiplier = 1.0,
    		double stromalStromalMultiplier = 1.0,
    		double epithelialStromalMultiplier = 1.0,
    		double apcTwoHitStromalMultiplier = 1.0);

    /**
     * Set whether to use an edge-based spring constant.
     *
     * @param useEdgeBasedSpringConstant
     */
    void SetEdgeBasedSpringConstant(bool useEdgeBasedSpringConstant);

    /* Set method for using one way springs between epithelial and tissue nodes
     */
    void SetUseOneWaySprings(bool useOneWaySprings = false);

    /* Set method for applying position-dependent spring strengths (i.e. for use at the base of the crypt)
     * @param useSpringConstantsMultiplier whether to multiply the spring constants by a factor
     * @param springConstantsMultiplier the multiplication factor for the spring constants
     */
    void SetPositionDependentSpringConstants(bool usePositionDependentSpringConstants = false, double springConstantsMultiplier = 1.0);

    /* Get method for spring constants multiplier
     */
    double GetPositionDependentSpringConstants();

    /* Removing duplicated entries of a vector
     */
    void RemoveDuplicates1D(std::vector<unsigned>& rVectorWithDuplicates);

    /* Returns a boolean for whether the element contains ghost nodes
     */

    bool DoesElementContainGhostNodes(AbstractCellPopulation<DIM>& rCellPopulation, unsigned elementIndex);

    /* Finding the number of elements that a node belongs to, which contain only real nodes
     * and not ghost nodes
     */
    unsigned GetNumContainingElementsWithoutGhostNodes(AbstractCellPopulation<DIM>& rCellPopulation, unsigned nodeIndex);

    /**
     * Return a multiplication factor for the spring constant, which
     * may depend on whether the given pair of neighbouring cells are
     * e.g. undergoing apoptosis, have mutations, or experience variable
     * levels of beta catenin.
     *
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeAGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the tissue
     * @param isCloserThanRestLength whether the neighbouring nodes lie closer than the rest length of their connecting spring
     *
     * @return the multiplication factor.
     */
    double VariableSpringConstantMultiplicationFactor(unsigned nodeAGlobalIndex,
                                                      unsigned nodeBGlobalIndex,
                                                      AbstractCellPopulation<DIM>& rCellPopulation,
                                                      bool isCloserThanRestLength);

    c_vector<double,2> GetCryptHeightExtremes(AbstractCellPopulation<DIM>& rCellPopulation);

    /*
     * Method to return the nodes connected to a particular node via the Delaunay
     * triangulation, excluding ghost nodes.
     */
    std::set<unsigned> GetNeighbouringNodeIndices(AbstractCellPopulation<DIM>& rCellPopulation, unsigned nodeIndex);

    /* Method to return a boolean that indicates whether this node/cell has detached from the basement membrane
     */
    bool HasEpithelialCellDetachedFromBasementMembrane(AbstractCellPopulation<DIM>& rCellPopulation, unsigned nodeIndex);

    /**
     * Overridden AddForceContribution method.
     *
     * @param rForces reference to vector of forces on nodes
     * @param rCellPopulation reference to the tissue
     */
    void AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                              AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden CalculateForceBetweenNodes() method.
     *
     * Calculates the force between two nodes.
     *
     * Note that this assumes they are connected and is called by rCalculateVelocitiesOfEachNode()
     *
     * Returns a zero force if one node is epithelial and the other is tissue and the spring is under tension. Returns
     * the force in the normal way if the spring is under compression.
     *
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the cell population
     *
     * @return The force exerted on Node A by Node B.
     */
    c_vector<double, DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Outputs force Parameters to file
	 *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PeriodicCryptModelInteractionForce)

#endif /* PERIODICCRYPTMODELINTERACTIONFORCE_HPP_ */
