#ifndef IREPRODUCTIONIMPLEMENTATION_H
#define IREPRODUCTIONIMPLEMENTATION_H
#include <GridCellCohortHandler.h>
#include <GridCellStockHandler.h>
#include <FunctionalGroupDefinitions.h>
#include <ProcessTracker.h>
#include <MadingleyModelInitialisation.h>
#include <map>
/** \file IReproductionImplementation.h
 * \brief the IReproductionImplementation header file
 */
/** \brief Interface for implementations of the ecological process of reproduction */
class IReproductionImplementation {
public:
    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------       
//Note emtpy function bodies - in principle it should be possible to make these pure virtual methods
    //linking under g++ is proving a bit tricky...
    //----------------------------------------------------------------------------------------------       
    /** \brief Generate new cohorts from reproductive potential mass
    @param gridCellCohorts The cohorts in the current grid cell 
    @param gridCellStocks The stocks in the current grid cell 
    @param actingCohort The position of the acting cohort in the jagged array of grid cell cohorts 
    @param cellEnvironment The environment in the current grid cell 
    @param madingleyCohortDefinitions The definitions for cohort functional groups in the model 
    @param madingleyStockDefinitions The definitions for stock functional groups in the model 
    @param currentTimestep The current model time step 
    @param trackProcesses An instance of ProcessTracker to hold diagnostics for eating 
    @param partial Thread-locked variables 
    @param iteroparous Whether the acting cohort is iteroparous, as opposed to semelparous 
    @param currentMonth The current model month */
    virtual void RunReproductionEvents(GridCell& gcl, Cohort& actingCohort,
            unsigned currentTimestep,ThreadLockedParallelVariables& partial, 
            bool iteroparous, unsigned currentMonth, MadingleyModelInitialisation& params) {
        cout<<"IReproductionImplementation RunReproductionEvents should be virtual: you probably don't want to be here"<<endl;
    }
    //----------------------------------------------------------------------------------------------       
    /** \brief Assigns surplus body mass to reproductive potential mass
    @param gridCellCohorts The cohorts in the current grid cell 
    @param gridCellStocks The stocks in the current grid cell 
    @param actingCohort The position of the acting cohort in the jagged array of grid cell cohorts 
    @param cellEnvironment The environment in the current grid cell 
    @param madingleyCohortDefinitions The definitions for cohort functional groups in the model 
    @param madingleyStockDefinitions The definitions for stock functional groups in the model 
    @param currentTimestep The current model time step 
    @param trackProcesses An instance of ProcessTracker to hold diagnostics for reproduction */
    virtual void RunReproductiveMassAssignment(GridCell& gcl, Cohort& actingCohort, unsigned currentTimestep, MadingleyModelInitialisation& params) {
          cout<<"IReproductionImplementation RunReproductiveMassAssignment should be virtual: you probably don't want to be here"<<endl;

    }
    //----------------------------------------------------------------------------------------------       
};
#endif
