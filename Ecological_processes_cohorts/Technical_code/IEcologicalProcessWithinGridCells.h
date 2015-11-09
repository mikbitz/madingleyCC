#ifndef IECOLOGICALPROCESSWITHINGRIDCELLS_H
#define IECOLOGICALPROCESSWITHINGRIDCELLS_H
#include <FunctionalGroupDefinitions.h>
#include <GridCellStockHandler.h>
#include <MadingleyModelInitialisation.h>
#include <ProcessTracker.h>
#include <ThreadLocked.h> 

/** \file IEcologicalProcessWithinGridCells.h
 * \brief the IEcologicalProcessWithinGridCells header file
 */

/** \brief Interface for ecological process code */


class IEcologicalProcessWithinGridCell {
public:
    //----------------------------------------------------------------------------------------------
    /** \brief Run the ecological process
    @param gridCellCohorts The cohorts in the current grid cell 
    @param gridCellStocks The stocks in the current grid cell 
    @param actingCohort The position of the acting cohort in the jagged array of grid cell cohorts 
    @param cellEnvironment The environment in the current grid cell 
    @param deltas The sorted list to track changes in biomass and abundance of the acting cohort in this grid cell 
    @param madingleyCohortHandler The definitions of cohort functional groups in the model 
    @param madingleyStockHandler The definitions of stock functional groups in the model 
    @param currentTimestep The current model time step 
    @param trackProcesses An instance of ProcessTracker to hold diagnostics for this ecological process 
    @param partial Thread-locked variables 
    @param outputDetail The level of output detail used for this model simulation 
    @param currentMonth The current model month 
    @param initialisation The instance of the MadingleyModelInitialisation class for this simulation 
     */
    virtual void RunEcologicalProcess(GridCellCohortHandler& gridCellCohorts,
            GridCellStockHandler& gridCellStocks,
            Cohort& actingCohort, map<string, vector<double> >& cellEnvironment,
            map<string, map<string, double>>&deltas,
            unsigned currentTimestep,
            ThreadLockedParallelVariables& partial,
            unsigned currentMonth, MadingleyModelInitialisation& params) {
        ;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Initialises an implementation of the ecological process
    @param gridCellCohorts The cohorts in the current grid cell 
    @param gridCellStocks The stocks in the current grid cell 
    @param madingleyCohortDefinitions The definitions for cohort functional groups in the model 
    @param madingleyStockDefinitions The definitions for stock functional groups in the model 
    @param implementationKey The name of the specific implementation of this process to initialize 
     */
    virtual void InitializeEcologicalProcess(GridCellCohortHandler& gridCellCohorts, GridCellStockHandler& gridCellStocks, FunctionalGroupDefinitions& madingleyCohortDefinitions, FunctionalGroupDefinitions& madingleyStockDefinitions, string implementationKey) {
        ;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Virtual destructor */
    virtual ~IEcologicalProcessWithinGridCell() {
        ;
    }
};

#endif
