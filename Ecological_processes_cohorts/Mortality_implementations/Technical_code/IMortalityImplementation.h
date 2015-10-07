#ifndef IMORTALITYIMPLEMENTATION_H
#define IMORTALITYIMPLEMENTATION_H
#include <GridCellCohortHandler.h>
/** \file IMortalityImplementation.h
 * \brief the IMortalityImplementation header file
 */


//namespace Madingley
//{

/** \brief Interface for implementations of the ecological process of mortality */
class IMortalityImplementation {
public:
    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------
    /** \brief Calculate the proportion of individuals in a cohort that die through a particular type of mortality in a model time step
    @param gridCellCohorts The cohorts in the current grid cell 
    @param actingCohort The position of the acting cohort in the jagged array of grid cell cohorts 
    @param bodyMassIncludingChangeThisTimeStep The body mass that individuals in this cohort will have at the end of this time step 
    @param deltas The sorted list to track changes in biomass and abundance of the acting cohort in this grid cell 
    @param currentTimestep The current model time step 
    @return The number of individuals lost to a cohort through mortality*/
    virtual double CalculateMortalityRate(GridCellCohortHandler& gridCellCohorts, vector<int>& actingCohort, double bodyMassIncludingChangeThisTimeStep, map<string, map<string, double>>&deltas, unsigned currentTimestep) {
        ;
    }
};
//}
#endif
