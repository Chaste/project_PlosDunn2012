#ifndef STROMALCELLMUTATIONSTATE_HPP_
#define STROMALCELLMUTATIONSTATE_HPP_

#include "AbstractCellMutationState.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * Subclass of AbstractCellMutationState defining a 'wild type' mutation state.
 */
class StromalCellMutationState : public AbstractCellMutationState
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell cycle model.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellMutationState>(*this);
    }

public:
    /**
     * Constructor.
     */
    StromalCellMutationState();

};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(StromalCellMutationState)

#endif /* STROMALCELLMUTATIONSTATE_HPP_ */
