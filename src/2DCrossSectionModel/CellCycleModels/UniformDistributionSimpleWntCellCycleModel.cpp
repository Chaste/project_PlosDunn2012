#include "UniformDistributionSimpleWntCellCycleModel.hpp"
#include "PetscTools.hpp"

UniformDistributionSimpleWntCellCycleModel::UniformDistributionSimpleWntCellCycleModel()
    : mUseCellProliferativeTypeDependentG1Duration(false),
      mWntStemThreshold(0.8),
      mWntTransitThreshold(0.65),
      mWntLabelledThreshold(0.65)
{
}


AbstractCellCycleModel* UniformDistributionSimpleWntCellCycleModel::CreateCellCycleModel()
{
    // Create a new cell cycle model
    UniformDistributionSimpleWntCellCycleModel* p_model = new UniformDistributionSimpleWntCellCycleModel();

    // Set the values of the member variables
    p_model->SetDimension(mDimension);
    p_model->SetCellProliferativeType(mCellProliferativeType);

    return p_model;
}


void UniformDistributionSimpleWntCellCycleModel::SetUseCellProliferativeTypeDependentG1Duration(bool useCellProliferativeTypeDependentG1Duration)
{
    mUseCellProliferativeTypeDependentG1Duration = useCellProliferativeTypeDependentG1Duration;
}


void UniformDistributionSimpleWntCellCycleModel::SetG1Duration()
{
    assert(mpCell != NULL);

    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    switch (mCellProliferativeType)
    {
        case STEM:
            if (mUseCellProliferativeTypeDependentG1Duration)
            {
                mG1Duration = 1 + 4*p_gen->ranf(); // U[1,5] according to Meineke
            }
            else
            {
                // Normally stem cells should behave just like transit cells in a Wnt simulation
                mG1Duration = 1 + 2*p_gen->ranf(); // U[1,3] according to Meineke
            }
            break;
        case TRANSIT:
            mG1Duration = 1 + 2*p_gen->ranf(); // U[1,3] according to Meineke
            break;
        case DIFFERENTIATED:
            mG1Duration = DBL_MAX;
            break;
        default:
            NEVER_REACHED;
    }

    // Check that the normal random deviate has not returned a small or negative G1 duration
    if (mG1Duration < mMinimumGapDuration)
    {
        mG1Duration = mMinimumGapDuration;
    }
}


double UniformDistributionSimpleWntCellCycleModel::GetWntLevel()
{
    assert(mpCell != NULL);
    double level = 0;

    switch (mDimension)
    {
        case 1:
        {
            const unsigned DIM = 1;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(mpCell);
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(mpCell);
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(mpCell);
            break;
        }
        default:
            NEVER_REACHED;
    }
    return level;
}


WntConcentrationType UniformDistributionSimpleWntCellCycleModel::GetWntType()
{
    WntConcentrationType wnt_type;
    switch (mDimension)
    {
        case 1:
        {
            const unsigned DIM = 1;
            wnt_type = WntConcentration<DIM>::Instance()->GetType();
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            wnt_type = WntConcentration<DIM>::Instance()->GetType();
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            wnt_type = WntConcentration<DIM>::Instance()->GetType();
            break;
        }
        default:
            NEVER_REACHED;
    }
    return wnt_type;
}


void UniformDistributionSimpleWntCellCycleModel::UpdateCellCyclePhase()
{
    // The cell can divide if the Wnt concentration >= wnt_division_threshold
    double wnt_division_threshold = DBL_MAX;

    // Set up under what level of Wnt stimulus a cell will divide
    double healthy_threshold = mWntTransitThreshold;
    double labelled_threshold = mWntLabelledThreshold;

    if (mpCell->GetMutationState()->IsType<WildTypeCellMutationState>())
    {
        wnt_division_threshold = healthy_threshold;
    }
    else if (mpCell->GetMutationState()->IsType<ApcOneHitCellMutationState>())
    {
        // should be less than healthy values
        wnt_division_threshold = 0.77*healthy_threshold;
    }
    else if (mpCell->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>())
    {
        // less than above value
        wnt_division_threshold = 0.155*healthy_threshold;
    }
    else if (mpCell->GetMutationState()->IsType<ApcTwoHitCellMutationState>())
    {
        // should be zero (no Wnt-dependence)
        wnt_division_threshold = 0.0;
    }
    else if (mpCell->GetMutationState()->IsType<StromalCellMutationState>())
    {
        // should be zero (no Wnt-dependence)
        wnt_division_threshold = 0.0;
    }
    else
    {
        NEVER_REACHED;
    }

    if (mpCell->HasCellProperty<CellLabel>())
    {
        wnt_division_threshold = labelled_threshold;
    }

    double wnt_level = GetWntLevel();
    WntConcentrationType wnt_type = GetWntType();

    // Set the cell type to TRANSIT if the Wnt stimulus exceeds wnt_division_threshold
    if (wnt_level >= wnt_division_threshold)
    {
        CellProliferativeType cell_type = TRANSIT;

        // For a RADIAL Wnt type, override the cell type to STEM if the Wnt stimulus exceeds a higher threshold
        if ( (wnt_type == RADIAL) && (wnt_level > mWntStemThreshold) )
        {
            cell_type = STEM;
        }

        mCellProliferativeType = cell_type;
    }
    else
    {
        // The cell is DIFFERENTIATED and so in G0 phase
        mCellProliferativeType = DIFFERENTIATED;
    }
    AbstractSimpleCellCycleModel::UpdateCellCyclePhase();
}

void UniformDistributionSimpleWntCellCycleModel::InitialiseDaughterCell()
{
    WntConcentrationType wnt_type = GetWntType();

    if (wnt_type == RADIAL)
    {
        mCellProliferativeType = TRANSIT;
    }

    AbstractSimpleCellCycleModel::InitialiseDaughterCell();
}

bool UniformDistributionSimpleWntCellCycleModel::CanCellTerminallyDifferentiate()
{
    return false;
}

double UniformDistributionSimpleWntCellCycleModel::GetWntStemThreshold()
{
    return mWntStemThreshold;
}

void UniformDistributionSimpleWntCellCycleModel::SetWntStemThreshold(double wntStemThreshold)
{
    assert(wntStemThreshold <= 1.0);
    assert(wntStemThreshold >= 0.0);
    mWntStemThreshold = wntStemThreshold;
}

double UniformDistributionSimpleWntCellCycleModel::GetWntTransitThreshold()
{
    return mWntTransitThreshold;
}

void UniformDistributionSimpleWntCellCycleModel::SetWntTransitThreshold(double wntTransitThreshold)
{
    assert(wntTransitThreshold <= 1.0);
    assert(wntTransitThreshold >= 0.0);
    mWntTransitThreshold = wntTransitThreshold;
}

double UniformDistributionSimpleWntCellCycleModel::GetWntLabelledThreshold()
{
    return mWntLabelledThreshold;
}

void UniformDistributionSimpleWntCellCycleModel::SetWntLabelledThreshold(double wntLabelledThreshold)
{
    assert(wntLabelledThreshold <= 1.0);
    assert(wntLabelledThreshold >= 0.0);
    mWntLabelledThreshold = wntLabelledThreshold;
}

void UniformDistributionSimpleWntCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile <<  "\t\t\t<WntStemThreshold>"<< mWntStemThreshold << "</WntStemThreshold>\n";
    *rParamsFile <<  "\t\t\t<WntTransitThreshold>"<< mWntTransitThreshold << "</WntTransitThreshold>\n";
    *rParamsFile <<  "\t\t\t<WntLabelledThreshold>"<< mWntLabelledThreshold << "</WntLabelledThreshold>\n";

    // Call direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(UniformDistributionSimpleWntCellCycleModel)
