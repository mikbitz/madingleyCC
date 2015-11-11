#ifndef MORTALITY_H
#define MORTALITY_H
#include <IMortalityImplementation.h>
#include <TBackgroundMortality.h>
#include <TSenescenceMortality.h>
#include <TStarvationMortality.h>
#include <limits>
/** \file Mortality.h
 * \brief the Mortality header file
 */

/** \brief  Performs mortality */
class Mortality : public IEcologicalProcessWithinGridCell {
public:
    //----------------------------------------------------------------------------------------------
    //Variables
    //----------------------------------------------------------------------------------------------
    /** \brief The available implementations of the mortality process */
    map<string, IMortalityImplementation*> Implementations;
    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------
    
    //----------------------------------------------------------------------------------------------
     /** \brief Constructor for Mortality: fills the list with available implementations of mortality
    @param globalModelTimeStepUnit The time step for the global model */
    Mortality(string globalModelTimeStepUnit) {

        // Add the background mortality implementation to the list of implementations
        BackgroundMortality* BackgroundMortalityImplementation = new BackgroundMortality(globalModelTimeStepUnit);
        Implementations["basic background mortality"] = BackgroundMortalityImplementation;

        // Add the senescence mortality implementation to the list of implementations
        SenescenceMortality* SenescenceMortalityImplementation = new SenescenceMortality(globalModelTimeStepUnit);
        Implementations["basic senescence mortality"] = SenescenceMortalityImplementation;

        // Add the starvation mortality implementation to the list of implementations
        StarvationMortality* StarvationMortalityImplementation = new StarvationMortality(globalModelTimeStepUnit);
        Implementations["basic starvation mortality"] = StarvationMortalityImplementation;

    }
    //----------------------------------------------------------------------------------------------
    /** \brief Destructor cleans up pointers
    */
    ~Mortality() {
        delete Implementations["basic background mortality"];
        delete Implementations["basic senescence mortality"];
        delete Implementations["basic starvation mortality"];
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Initialize an implementation of mortality. This is only in here to satisfy the requirements of IEcologicalProcessAcrossGridCells
    @param gridCellCohorts The cohorts in the current grid cell 
    @param gridCellStocks The stocks in the current grid cell 
    @param madingleyCohortDefinitions The definitions for cohort functional groups in the model 
    @param madingleyStockDefinitions The definitions for stock functional groups in the model 
    @param implementationKey The name of the implementation of mortality to initialize */
    void InitializeEcologicalProcess(GridCell& gcl, MadingleyModelInitialisation& params, string implementationKey) {

    }
    //----------------------------------------------------------------------------------------------
    /** \brief Run mortality
    @param gridCellCohorts The cohorts in the current grid cell 
    @param gridCellStocks The stocks in the current grid cell 
    @param actingCohort The position of the acting cohort in the jagged array of grid cell cohorts 
    @param cellEnvironment The environment in the current grid cell 
    @param madingleyCohortDefinitions The definitions for cohort functional groups in the model 
    @param madingleyStockDefinitions The definitions for stock functional groups in the model 
    @param currentTimestep The current model time step 
    @param trackProcesses An instance of ProcessTracker to hold diagnostics for mortality 
    @param partial Thread-locked variables 
    @param specificLocations Whether the model is being run for specific locations 
    @param outputDetail The level output detail being used for the current model run 
    @param currentMonth The current model month */
    void RunEcologicalProcess(GridCell& gcl,
            Cohort& actingCohort, 
            unsigned currentTimestep,
            ThreadLockedParallelVariables& partial,
            unsigned currentMonth, MadingleyModelInitialisation& params) {

        // Variables to hold the mortality rates
        double MortalityRateBackground;
        double MortalityRateSenescence;
        double MortalityRateStarvation;

        // Variable to hold the total abundance lost to all forms of mortality
        double MortalityTotal;

        // Individual body mass including change this time step as a result of other ecological processes
        double BodyMassIncludingChangeThisTimeStep;

        // Individual reproductive mass including change this time step as a result of other ecological processes
        double ReproductiveMassIncludingChangeThisTimeStep;

        // Calculate the body mass of individuals in this cohort including mass gained through eating this time step, up to but not exceeding adult body mass for this cohort. 
        // Should be fine because these deductions are made in the reproduction implementation, but use Math.Min to double check.

        BodyMassIncludingChangeThisTimeStep = 0.0;

        // Loop over all items in the biomass deltas
        for (auto Biomass : Cohort::Deltas["biomass"]) {
            // Add the delta biomass to net biomass
            BodyMassIncludingChangeThisTimeStep += Biomass.second;
        }
        BodyMassIncludingChangeThisTimeStep = min(actingCohort.AdultMass, BodyMassIncludingChangeThisTimeStep + actingCohort.IndividualBodyMass);

        // Temporary variable to hold net reproductive biomass change of individuals in this cohort as a result of other ecological processes
        ReproductiveMassIncludingChangeThisTimeStep = 0.0;

        // Loop over all items in the biomass Cohort::Deltas
        for (auto Biomass : Cohort::Deltas["reproductivebiomass"]) {
            // Add the delta biomass to net biomass
            ReproductiveMassIncludingChangeThisTimeStep += Biomass.second;
        }

        ReproductiveMassIncludingChangeThisTimeStep += actingCohort.IndividualReproductivePotentialMass;

        // Check to see if the cohort has already been killed by predation etc
        if (BodyMassIncludingChangeThisTimeStep <= 1.e-15) //MB a small number ! maybe should be larger?
        {
            // If individual body mass is not greater than zero, then all individuals become extinct
            MortalityTotal = actingCohort.CohortAbundance;
        } else {
            // Calculate background mortality rate
            MortalityRateBackground = Implementations["basic background mortality"]->CalculateMortalityRate(
                    actingCohort, BodyMassIncludingChangeThisTimeStep,  currentTimestep);

            // If the cohort has matured, then calculate senescence mortality rate, otherwise set rate to zero
            if (actingCohort.MaturityTimeStep != std::numeric_limits<unsigned>::max()) {
                MortalityRateSenescence = Implementations["basic senescence mortality"]->CalculateMortalityRate(
                        actingCohort, BodyMassIncludingChangeThisTimeStep,  currentTimestep);
            } else {
                MortalityRateSenescence = 0.0;
            }

            // Calculate the starvation mortality rate based on individual body mass and maximum body mass ever
            // achieved by this cohort
            MortalityRateStarvation = Implementations["basic starvation mortality"]->CalculateMortalityRate(actingCohort, BodyMassIncludingChangeThisTimeStep,  currentTimestep);

            // Calculate the number of individuals that suffer mortality this time step from all sources of mortality
            MortalityTotal = (1 - exp(-MortalityRateBackground - MortalityRateSenescence -
                    MortalityRateStarvation)) * actingCohort.CohortAbundance;
        }

        // Remove individuals that have died from the delta abundance for this cohort
        Cohort::Deltas["abundance"]["mortality"] = -MortalityTotal;

        // Add the biomass of individuals that have died to the delta biomass in the organic pool (including reproductive 
        // potential mass, and mass gained through eating, and excluding mass lost through metabolism)
        Cohort::Deltas["organicpool"]["mortality"] = MortalityTotal * (BodyMassIncludingChangeThisTimeStep + ReproductiveMassIncludingChangeThisTimeStep);
    }
    //----------------------------------------------------------------------------------------------
};
#endif
