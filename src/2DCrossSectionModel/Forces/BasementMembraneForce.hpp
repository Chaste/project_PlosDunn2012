#ifndef BASEMENTMEMBRANEFORCE_HPP_
#define BASEMENTMEMBRANEFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "StromalCellMutationState.hpp"
#include "MeshBasedCellPopulation.hpp"

#include <cmath>
#include <list>
#include <fstream>

/**
 * A force class that defines the force due to the basement membrane.
 */

class BasementMembraneForce : public AbstractForce<2>
{
    friend class TestCrossSectionModelInteractionForce;

private :

    /** Parameter that multiplies the curvature to give the basement membrane force */
    double mBasementMembraneParameter;

    /** Target curvature for the curved base of the crypt (the actin basket region) */
    double mCryptBaseCurvature;

    /** x-coordinate that encloses the region in which to apply a non-zero target curvature */
    double mLeftBoundary;

    /** x-coordinate that encloses the region in which to apply a non-zero target curvature */
    double mRightBoundary;

    /** Make the basement membrane force dependent on the position of a cell up the crypt */
    bool mUsePositionDependentMembraneForce;

    /** The multiplication factor for the basement membrane parameter */
    double mMembraneForceMultiplier;

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
        archive & boost::serialization::base_object<AbstractForce<2> >(*this);
        archive & mBasementMembraneParameter;
        archive & mCryptBaseCurvature;
        archive & mLeftBoundary;
        archive & mRightBoundary;
        archive & mUsePositionDependentMembraneForce;
        archive & mMembraneForceMultiplier;
    }

public :

    /**
     * Constructor.
     */
	BasementMembraneForce();

    /**
     * Destructor.
     */
    ~BasementMembraneForce();

    /* Set method for Basement Membrane Parameter
     */
    void SetBasementMembraneParameter(double basementMembraneParameter);

    /* Get method for Basement Membrane Parameter
     */
    double GetBasementMembraneParameter();

    /* Value of curvature at crypt base - left and right boundary used for the box model set up - don't set for the cross
     * section model */
    void SetCryptBaseCurvature(double cryptBaseCurvature = 0.0, double leftBoundary = 0.0, double rightBoundary = 100.0);

    double GetCryptBaseCurvature();

    /* Set method for position-dependent basement membrane force multiplier (i.e. if you want to apply a different
     * basement membrane parameter in the crypt base, or at the orifice)
     * @param usePositionDependentMembraneForce whether to multiply the basement membrane force by a factor
     * @param membraneForceMultiplier the multiplication factor for the basement membrane force
     */
    void SetPositionDependentMultiplier(bool usePositionDependentMembraneForce = false, double membraneForceMultiplier = 1.0);

    /* Get method for basement membrane force strength multiplier
     */
    double GetPositionDependentMultiplier();

    /* Removing duplicated entries of a vector
     */
    void RemoveDuplicates1D(std::vector<unsigned>& rVectorWithDuplicates);

    /* Returns a boolean for whether the element contains ghost nodes
     */
    bool DoesElementContainGhostNodes(AbstractCellPopulation<2>& rCellPopulation, unsigned elementIndex);

    /* Finding the connected pairs of epithelial-tissue nodes
     */
    std::vector<c_vector<unsigned, 2> > GetEpithelialTissuePairs(AbstractCellPopulation<2>& rCellPopulation);

    /* Takes an epithelial node index and a tissue node index and returns the curvature of
     * the curve passing through the midpoints of the epithelial-tissue springs of the
     * common elements
     */
    double GetCurvatureFromNodePair(AbstractCellPopulation<2>& rCellPopulation, unsigned epithelialNodeIndex,
    		unsigned tissueNodeIndex);

    /*
     * Finding the curvature between three midpoints parametrically - in this case, we find the normal
     * to the vector joining the left and right midpoints, and then find the perpendicular distance of
     * the centre midpoint from the left->right vector
     */
    double GetCurvatureFromMidpoints(AbstractCellPopulation<2>& rCellPopulation, c_vector<double, 2> leftMidpoint,
    														c_vector<double, 2> centreMidpoint,
    														c_vector<double, 2> rightMidpoint);

    double FindParametricCurvature(c_vector<double, 2> leftMidpoint,
									c_vector<double, 2> centreMidpoint,
									c_vector<double, 2> rightMidpoint);

    /* Finding the number of elements that a node belongs to, which contain only real nodes
     * and not ghost nodes
     */
    unsigned GetNumContainingElementsWithoutGhostNodes(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex);

    /* Finding the y-coordinates of the crypt orifice and crypt base
     * The first entry of the resulting vector is the orifice, the second is the base
     */
    c_vector<double,2> GetCryptHeightExtremes(AbstractCellPopulation<2>& rCellPopulation);

    /*
     * Method to return the nodes connected to a particular node via the Delaunay
     * triangulation, excluding ghost nodes.
     */
    std::set<unsigned> GetNeighbouringNodeIndices(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex);

    /* Method to return a boolean that indicates whether this node/cell has detached from the basement membrane
     */
    bool HasEpithelialCellDetachedFromBasementMembrane(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex);

    /**
     * Overridden AddForceContribution method.
     *
     * @param rForces reference to vector of forces on nodes
     * @param rCellPopulation reference to the tissue
     */
    void AddForceContribution(std::vector<c_vector<double, 2> >& rForces,
                              AbstractCellPopulation<2>& rCellPopulation);

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
CHASTE_CLASS_EXPORT(BasementMembraneForce)

#endif /*BASEMENTMEMBRANEFORCE_HPP_*/
