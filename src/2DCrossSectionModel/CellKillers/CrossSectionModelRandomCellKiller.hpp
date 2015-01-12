#ifndef CROSSSECTIONMODELRANDOMCELLKILLER_HPP_
#define CROSSSECTIONMODELRANDOMCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "OutputFileHandler.hpp"
#include "StromalCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApoptoticCellProperty.hpp"
#include "SimpleWntCellCycleModel.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/* A cell killer that kills transit cells by anoikis, and will slough cells
 * from the crypt orifice randomly
 */
class CrossSectionModelRandomCellKiller : public AbstractCellKiller<2>
{
private:

    /** Whether cells should be sloughed randomly from the top of the crypt. */
    bool mSloughOrifice;

    /**
     * Probability that in an hour's worth of trying, this cell killer
     * will have successfully killed a given cell.
     */
    double mProbabilityOfDeathInAnHour;

    unsigned mCellsRemovedByAnoikis;

    unsigned mCellsRemovedRandomly;

    std::vector<c_vector<double,3> > mLocationsOfAnoikisCells;

    // The output file directory for the simulation data that corresponds to the number of cells
    // killed by anoikis
    out_stream mAnoikisOutputFile;

    std::string mOutputDirectory;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * Serialization of singleton objects must be done with care.
     * Before the object is serialized via a pointer, it *MUST* be
     * serialized directly, or an assertion will trip when a second
     * instance of the class is created on de-serialization.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<2> >(*this);

        archive & mSloughOrifice; // done in load_construct_data
        archive & mProbabilityOfDeathInAnHour;
        archive & mCellsRemovedByAnoikis;
        archive & mCellsRemovedRandomly;
        archive & mOutputDirectory;
    }

public:

    /**
     * Default constructor.
     *
     * @param pCellPopulation pointer to a tissue
     * @param sloughOrifice whether to slough compressed cells at crypt orifice
     */
	CrossSectionModelRandomCellKiller(AbstractCellPopulation<2>* pCellPopulation, bool sloughOrifice=true, double probabilityOfDeathInAnHour = 0.1);

	// Destructor
	~CrossSectionModelRandomCellKiller();

    bool GetSloughOrifice() const;

    /* @return mProbabilityOfDeathInAnHour */
    double GetDeathProbabilityInAnHour() const;

    void SetOutputDirectory(std::string outputDirectory);

    std::string GetOutputDirectory();

    std::set<unsigned> GetNeighbouringNodeIndices(unsigned nodeIndex);

    bool HasCellPoppedUp(unsigned nodeIndex);

    std::vector<c_vector<unsigned,3> > RemoveByAnoikisOrRandomApoptosis();

    /**
     *  Loops over and kills cells by anoikis or at the orifice if instructed.
     */
    void TestAndLabelCellsForApoptosisOrDeath();

    /* After each event of cell killing in TestAndLabelCellsForApoptosisOrDeath(), the information of whether to kill each cell
     * or not is passed to this method which then increments the member variables corresponding to the total number of cells
     * killed by anoikis or apoptosis through compression
     */
    void SetNumberCellsRemoved(std::vector<c_vector<unsigned,3> > cellsRemoved);

    /* Returns the total number of cells removed by anoikis ([0]) and by compression ([1])
     *
     */
    c_vector<unsigned,2> GetNumberCellsRemoved();

    /* Storing the x-locations of those epithelial cells that get removed by anoikis
     *
     */
    void SetLocationsOfCellsRemovedByAnoikis(std::vector<c_vector<unsigned,3> > cellsRemoved);

    /* Returns the x-coordinates of those cells removed by anoikis
     *
     */
    std::vector<c_vector<double,3> > GetLocationsOfCellsRemovedByAnoikis();

    /* Method to return the current coordinates of the crypt orifice and
     * crypt base - these can be used to accurately define the region of the
     * crypt base. (This will be the y-coordinate in 2D, or the z coordinate in 3D)
     * [0] - y or z-coordinate of orifice
     * [1] - y or z-coordinate of base
     */
    c_vector<double,2> GetCryptHeightExtremes();

    /**
     * Outputs cell killer parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellKillerParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CrossSectionModelRandomCellKiller)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CrossSectionModelRandomCellKiller.
 */
template<class Archive >
inline void save_construct_data(
    Archive & ar, const CrossSectionModelRandomCellKiller * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<2>* const p_tissue = t->GetCellPopulation();
    ar << p_tissue;
    bool slough_orifice = t->GetSloughOrifice();
    ar << slough_orifice;
    double death_probability = t->GetDeathProbabilityInAnHour();
    ar << death_probability;
}

/**
 * De-serialize constructor parameters and initialise Crypt.
 */
template<class Archive >
inline void load_construct_data(
    Archive & ar, CrossSectionModelRandomCellKiller * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<2>* p_tissue;
    ar >> p_tissue;
    bool slough_orifice;
    ar >> slough_orifice;
    double death_probability = t->GetDeathProbabilityInAnHour();
    ar >> death_probability;

    // Invoke inplace constructor to initialise instance
    ::new(t)CrossSectionModelRandomCellKiller(p_tissue, slough_orifice, death_probability);
}
}
} // namespace ...

#endif /* CROSSSECTIONMODELRANDOMCELLKILLER_HPP_ */
