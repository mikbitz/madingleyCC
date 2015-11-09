#ifndef ECOLOGYCOHORT_H
#define ECOLOGYCOHORT_H
#include <ApplyEcology.h>
#include <IEcologicalProcessWithinGridCells.h>
#include <IEatingImplementation.h>
#include <Eating.h>
#include <Reproduction.h>
#include <Mortality.h>
#include <Metabolism.h>
/** \file EcologyCohort.h
 * \brief the EcologyCohort header file
 */

//
//namespace Madingley
//{

/** \brief A class to specify, initalise and run ecological processes pertaining to cohorts */
class EcologyCohort {
public:
    //----------------------------------------------------------------------------------------------
    //Variables
    //----------------------------------------------------------------------------------------------
    /** \brief  A vector of stopwatch objects for timing the ecological processes*/
    vector<StopWatch> s2;
    /** \brief A sorted list of formulations of metabolism */
    map<string, IEcologicalProcessWithinGridCell*> MetabolismFormulations;
    /** \brief A sorted list of formulations of eating */
    map<string, IEcologicalProcessWithinGridCell*> EatingFormulations;
    /** \brief A sorted list of formulations of mortality */
    map<string, IEcologicalProcessWithinGridCell*> MortalityFormulations;
    /** \brief A sorted list of formulations of reproduction */
    map<string, IEcologicalProcessWithinGridCell*> ReproductionFormulations;
    /** \brief An instance of apply ecology */
    ApplyEcology ApplyEcologicalProcessResults;
    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------

    //----------------------------------------------------------------------------------------------
    /** \brief Initalise the ecological processes */
    void InitializeEcology(double cellArea, string globalModelTimeStepUnit, bool drawRandomly) {
        // Declare and attach eating formulations
        Eating *EatingFormulation = new Eating(cellArea, globalModelTimeStepUnit);
        EatingFormulations["Basic eating"] = EatingFormulation;
        // Declare and attach metabolism formulations
        Metabolism *MetabolismFormulation = new Metabolism(globalModelTimeStepUnit);
        MetabolismFormulations["Basic metabolism"] = MetabolismFormulation;
        // Declare and attach mortality formulations
        Mortality *MortalityFormulation = new Mortality(globalModelTimeStepUnit);
        MortalityFormulations["Basic mortality"] = MortalityFormulation;
        // Declare and attach mortality formulations
        Reproduction *ReproductionFormulation = new Reproduction(globalModelTimeStepUnit, drawRandomly);
        ReproductionFormulations["Basic reproduction"] = ReproductionFormulation;
    }
    //----------------------------------------------------------------------------------------------
    ~EcologyCohort() {
        delete EatingFormulations["Basic eating"];
        delete MetabolismFormulations["Basic metabolism"];
        delete MortalityFormulations["Basic mortality"];
        delete ReproductionFormulations["Basic reproduction"];
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Run ecological processes that operate on cohorts within a single grid cell
    @param gridCellCohorts The cohorts in the current grid cell 
    @param gridCellStocks The stocks in the current grid cell 
    @param actingCohort The acting cohort 
    @param cellEnvironment The environment in the current grid cell 
    @param deltas A sorted list of deltas to track changes in abundances and biomasses during the ecological processes 
    @param madingleyCohortDefinitions The definitions for cohort functional groups in the model 
    @param madingleyStockDefinitions The definitions for stock functional groups in the model 
    @param currentTimestep The current model time step 
    @param trackProcesses An instance of the process tracker 
    @param partial Thread-locked local variables 
    @param outputDetail The level of output detail being used for this model run 
    @param currentMonth The current model month */
    void RunWithinCellEcology(GridCell& gcl, Cohort& actingCohort, unsigned currentTimestep, 
            ThreadLockedParallelVariables& partial,  unsigned currentMonth, MadingleyModelInitialisation& params) {

        // RUN EATING
        EatingFormulations["Basic eating"]->RunEcologicalProcess(gcl.GridCellCohorts, gcl.GridCellStocks, actingCohort, gcl.CellEnvironment,
                gcl.Deltas,  currentTimestep, partial,currentMonth, params);


        // RUN METABOLISM - THIS TIME TAKE THE METABOLIC LOSS TAKING INTO ACCOUNT WHAT HAS BEEN INGESTED THROUGH EATING
        MetabolismFormulations["Basic metabolism"]->RunEcologicalProcess(gcl.GridCellCohorts,gcl.GridCellStocks, actingCohort,
                gcl.CellEnvironment, gcl.Deltas,  currentTimestep, partial,currentMonth, params);


        // RUN REPRODUCTION - TAKING INTO ACCOUNT NET BIOMASS CHANGES RESULTING FROM EATING AND METABOLISING
        ReproductionFormulations["Basic reproduction"]->RunEcologicalProcess(gcl.GridCellCohorts, gcl.GridCellStocks, actingCohort,
                gcl.CellEnvironment, gcl.Deltas, currentTimestep, partial,currentMonth, params);


        // RUN MORTALITY - TAKING INTO ACCOUNT NET BIOMASS CHANGES RESULTING FROM EATING, METABOLISM AND REPRODUCTION
        MortalityFormulations["Basic mortality"]->RunEcologicalProcess(gcl.GridCellCohorts, gcl.GridCellStocks, actingCohort,
                gcl.CellEnvironment, gcl.Deltas,  currentTimestep, partial, currentMonth, params);
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Update the properties of the acting cohort and of the environmental biomass pools after running the ecological processes for a cohort
    @param gridCellCohorts The cohorts in the current grid cell 
    @param gridCellStocks The stocks in the current grid cell 
    @param actingCohort The acting cohort 
    @param cellEnvironment The environment of the current grid cell 
    @param deltas The sorted list of deltas for the current grid cell 
    @param madingleyCohortDefinitions The definitions for cohort functional groups in the model 
    @param madingleyStockDefinitions The definitions for stock functional groups in the model 
    @param currentTimestep The current model time step 
    @param tracker A process tracker */
    void UpdateEcology(GridCell& gcl, Cohort& actingCohort, unsigned currentTimestep) {
        // Apply the results of within-cell ecological processes
        ApplyEcologicalProcessResults.UpdateAllEcology(gcl.GridCellCohorts, actingCohort, gcl.CellEnvironment, gcl.Deltas, currentTimestep);
    }
    //----------------------------------------------------------------------------------------------
};
#endif
