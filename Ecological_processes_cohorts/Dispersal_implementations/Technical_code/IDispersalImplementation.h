#ifndef IDISPERSALIMPLEMENTATION_H
#define IDISPERSALIMPLEMENTATION_H
#include <ModelGrid.h>
#include <Cohort.h>

/** \file IDispersalImplementation.h
 * \brief the IDispersalImplementation header file
 */

/** \brief Interface for implementations of the ecological process of dispersal */
class IDispersalImplementation {
public:
    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------

    //----------------------------------------------------------------------------------------------
     /** \brief Run the dispersal implementation */
    virtual void RunDispersal(vector<unsigned>& cellIndex, ModelGrid& gridForDispersal, Cohort& cohortToDisperse, int actingCohortFunctionalGroup, int actingCohortNumber, unsigned currentMonth) {
        ;
    }
};
#endif
