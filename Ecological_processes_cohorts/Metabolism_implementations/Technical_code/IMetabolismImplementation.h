#ifndef IMETABOLISMIMPLEMENTATION_H
#define IMETABOLISMIMPLEMENTATION_H
#include <GridCellCohortHandler.h>
#include <GridCellStockHandler.h>
#include <FunctionalGroupDefinitions.h>
/** \file IMetabolismImplementation.h
 * \brief the IMetabolismImplementation header file
 */

/** \brief Interface for implementations of the ecological process of metabolism */
class IMetabolismImplementation {
public:
    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------

    //----------------------------------------------------------------------------------------------
    /** \brief Calculate the biomass lost through metabolism and update the relevant deltas for the acting cohort
    @param gridCellCohorts The cohorts in the current grid cell 
    @param gridCellStocks The stocks in the current grid cell 
    @param actingCohort The position of the acting cohort in the jagged array of grid cell cohorts 
    @param cellEnvironment The environment in the current grid cell 
    @param madingleyCohortDefinitions The definitions for cohort functional groups in the model 
    @param madingleyStockDefinitions The definitions for stock functional groups in the model 
    @param currentTimestep The current model time step 
    @param currentMonth The current month in the model */
    virtual void RunMetabolism(  Cohort& actingCohort,
            map<string, vector<double> >& cellEnvironment, 
            unsigned currentTimestep, unsigned currentMonth) {
        ;
    }
    //----------------------------------------------------------------------------------------------
};
#endif
