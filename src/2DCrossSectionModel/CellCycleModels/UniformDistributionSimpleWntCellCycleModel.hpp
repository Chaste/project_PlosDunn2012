#ifndef UNIFORMDISTRIBUTIONSIMPLEWNTCELLCYCLEMODEL_HPP_
#define UNIFORMDISTRIBUTIONSIMPLEWNTCELLCYCLEMODEL_HPP_

#include "AbstractSimpleCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "WntConcentration.hpp"

// Needed here to avoid serialization errors
#include "WildTypeCellMutationState.hpp"
#include "CellLabel.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "StromalCellMutationState.hpp"


/**
 *  Simple Wnt-dependent cell cycle model that samples from a UNIFORM distribution, rather than
 *  a normal distribution, as is the case in SimpleWntCellCycleModel
 */
class UniformDistributionSimpleWntCellCycleModel : public AbstractSimpleCellCycleModel
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell cycle model, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleCellCycleModel>(*this);

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        archive & *p_gen;

        archive & mUseCellProliferativeTypeDependentG1Duration;
        archive & mWntStemThreshold;
        archive & mWntTransitThreshold;
        archive & mWntLabelledThreshold;
    }

    /**
     * Whether to use different mean G1 durations for different cell types.
     * For use in SetG1Duration().
     */
    bool mUseCellProliferativeTypeDependentG1Duration;

    /**
     * Non-dimensionalized Wnt threshold, above which cells behave as stem cells.
     */
    double mWntStemThreshold;

    /**
     * Non-dimensionalized Wnt threshold, above which cells progress through the cell cycle.
     */
    double mWntTransitThreshold;

    /**
     * Non-dimensionalized Wnt threshold, above which labelled cells progress through the cell cycle.
     */
    double mWntLabelledThreshold;

protected:

    /**
     * Get the Wnt level experienced by the cell.
     */
    double GetWntLevel();

    /**
     * Get the type of Wnt concentration (LINEAR, RADIAL, or NONE).
     * This affects how the cell cycle phase is updated.
     */
    WntConcentrationType GetWntType();

    /**
     * Stochastically set the G1 duration. The G1 duration is taken
     * from a normal distribution whose mean is the G1 duration given
     * in AbstractCellCycleModel for the cell type and whose standard deviation
     * is 1.
     *
     * Called on cell creation at the start of a simulation, and for both
     * parent and daughter cells at cell division.
     */
    void SetG1Duration();

public:

    /**
     * Constructor - just a default, mBirthTime is now set in the AbstractCellCycleModel class.
     * mG1Duration is set very high, it is set for the individual cells when InitialiseDaughterCell is called.
     */
    UniformDistributionSimpleWntCellCycleModel();

    /**
     * Overridden UpdateCellCyclePhase() method.
     */
    virtual void UpdateCellCyclePhase();

    /**
     * Overridden InitialiseDaughterCell() method.
     */
    virtual void InitialiseDaughterCell();

    /**
     * Overridden builder method to create new copies of
     * this cell cycle model.
     */
    virtual AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Set whether Whether the duration of the G1 phase is dependent on cell type
     * @param useCellProliferativeTypeDependentG1Duration - boolean, defaults to true.
     */
    void SetUseCellProliferativeTypeDependentG1Duration(bool useCellProliferativeTypeDependentG1Duration=true);

    /**
     * Overridden CanCellTerminallyDifferentiate() method.
     */
    virtual bool CanCellTerminallyDifferentiate();

    /**
     * @return mWntStemThreshold
     */
    double GetWntStemThreshold();

    /**
     * Set mWntStemThreshold.
     *
     * @param wntStemThreshold the value of mWntStemThreshold
     */
    void SetWntStemThreshold(double wntStemThreshold);

    /**
     * @return mWntTransitThreshold
     */
    double GetWntTransitThreshold();

    /**
     * Set mWntTransitThreshold.
     *
     * @param wntTransitThreshold the value of mWntTransitThreshold
     */
    void SetWntTransitThreshold(double wntTransitThreshold);

    /**
     * @return mWntLabelledThreshold
     */
    double GetWntLabelledThreshold();

    /**
     * Set mWntLabelledThreshold.
     *
     * @param wntLabelledThreshold the value of mWntLabelledThreshold
     */
    void SetWntLabelledThreshold(double wntLabelledThreshold);

    /**
     * Outputs cell cycle model parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(UniformDistributionSimpleWntCellCycleModel)


#endif /* UNIFORMDISTRIBUTIONSIMPLEWNTCELLCYCLEMODEL_HPP_ */
