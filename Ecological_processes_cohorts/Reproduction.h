#ifndef REPRODUCTION_H
#define REPRODUCTION_H
#include <IReproductionImplementation.h>
#include <IEcologicalProcessWithinGridCells.h>
#include <TReproductionBasic.h>

/** \file Reproduction.h
 * \brief the Reproduction header file
 */

/** \brief Performs reproduction */
class Reproduction : public IEcologicalProcessWithinGridCell {
public:
    //----------------------------------------------------------------------------------------------
    //Variables
    //----------------------------------------------------------------------------------------------
    /** \brief The available implementations of the reproduction process */
    map<string, IReproductionImplementation*> Implementations;
    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------
    
    //----------------------------------------------------------------------------------------------
    /**  \brief Constructor for Reproduction: fills the list of available implementations of reproduction */
    Reproduction(string globalModelTimeStepUnit, bool drawRandomly) {
        // Add the basic reproduction implementation to the list of implementations
        ReproductionBasic* ReproductionImplementation = new ReproductionBasic(globalModelTimeStepUnit, drawRandomly);
        Implementations["reproduction basic"] = ReproductionImplementation;
    }
    //----------------------------------------------------------------------------------------------
    /** Destructor ensure we tidy everything up */
    ~Reproduction() {
        delete Implementations["reproduction basic"];
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Initialize an implementation of reproduction. This is only in here to satisfy the requirements of IEcologicalProcessWithinGridCells

    @param gridCellCohorts The cohorts in the current grid cell 
    @param gridCellStocks The stocks in the current grid cell 
    @param madingleyCohortDefinitions The definitions for cohort functional groups in the model 
    @param madingleyStockDefinitions The definitions for stock functional groups in the model 
    @param implementationKey The name of the reproduction implementation to initialize 
     */
    void InitializeEcologicalProcess(GridCell& gcl, MadingleyModelInitialisation& params, string implementationKey){
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Run reproduction
    @param gridCellCohorts The cohorts in the current grid cell 
    @param gridCellStocks The stocks in the current grid cell 
    @param actingCohort The position of the acting cohort in the jagged array of grid cell cohorts 
    @param cellEnvironment The environment in the current grid cell 
    @param madingleyCohortDefinitions The definitions for cohort functional groups in the model 
    @param madingleyStockDefinitions The definitions for stock functional groups in the model 
    @param currentTimeStep The current model time step 
    @param processTracker An instance of ProcessTracker to hold diagnostics for eating 
    @param partial Thread-locked variables for the parallelised version 
    @param specificLocations Whether the model is being run for specific locations 
    @param outputDetail The level of output detail being used for this model run 
    @param currentMonth The current model month 
     */
    void RunEcologicalProcess(GridCell& gcl,
            Cohort& actingCohort, 
            unsigned currentTimestep,
            ThreadLockedParallelVariables& partial,
            unsigned currentMonth, MadingleyModelInitialisation& params) {

        // Holds the reproductive strategy of a cohort
        bool _Iteroparous = params.CohortFunctionalGroupDefinitions.GetTraitNames("reproductive strategy", actingCohort.FunctionalGroupIndex) == "iteroparity";

        // Assign mass to reproductive potential
        Implementations["reproduction basic"]->RunReproductiveMassAssignment(gcl, actingCohort, currentTimestep, params);

        // Run reproductive events. Note that we can't skip juveniles here as they could conceivably grow to adulthood and get enough biomass to reproduce in a single time step
        // due to other ecological processes
        Implementations["reproduction basic"]->RunReproductionEvents(gcl, actingCohort, currentTimestep, partial, _Iteroparous, currentMonth,params);
    }
    //----------------------------------------------------------------------------------------------
};
#endif
