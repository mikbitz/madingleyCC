#ifndef TREVISEDREPRODUCTION_H
#define TREVISEDREPRODUCTION_H
/** \file TRevisedReproduction.h
 * \brief the TRevisedReproduction header file
 */
///MB this class seems incomplete at present...
/** \brief A formulation of the process of reproduction */
class RevisedReproduction : IReproductionImplementation {
    //----------------------------------------------------------------------------------------------
    //Variables
    //----------------------------------------------------------------------------------------------  
    /** \brief Scalar to convert from the time step units used by this formulation of reproduction to global model time step units */
    double DeltaT;
    /** \brief Include Utility class */
    UtilityFunctions Utilities;
    /** \brief An instance of the simple random number generator class */
    std::default_random_engine RandomNumberGenerator;
public:
    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------  

    //----------------------------------------------------------------------------------------------  
    /** \brief
    Constructor for reproduction: assigns all parameter values
     */
    RevisedReproduction() {
        // Initialise ecological parameters for reproduction
        InitializeReproductionParameters();

        // Calculate the scalar to convert from the time step units used by this implementation of herbivory to the global model time step units
        DeltaT = UtilityFunctions.ConvertTimeUnits(MadingleyModel.GlobalModelTimeStepUnit, TimeUnitImplementation);
    }
    //----------------------------------------------------------------------------------------------  
    /** \brief    Generate new cohorts from reproductive potential mass
    @param gridCellCohorts The cohorts in the current grid cell 
    @param gridCellStocks The stocks in the current grid cell 
    @param actingCohort The position of the acting cohort in the jagged array of grid cell cohorts 
    @param cellEnvironment The environment of the current grid cell 
    @param madingleyCohortDefinitions The definitions of cohort functional groups in the model 
    @param madingleyStockDefinitions The definitions of stock functional groups in the model 
    @param currentTimestep The current model time step 
    @param tracker An instance of ProcessTracker to hold diagnostics for reproduction 
    @param partial Thread-locked variables */
    void RunReproduction(GridCellCohortHandler& gridCellCohorts, GridCellStockHandler& gridCellStocks,
            Cohort& actingCohort, map<string, vector<double> >& cellEnvironment, map<string,
            FunctionalGroupDefinitions& madingleyCohortDefinitions, FunctionalGroupDefinitions& madingleyStockDefinitions,
            unsigned currentTimestep, ProcessTracker& tracker, ThreadLockedParallelVariables& partial) {
        // Check that the abundance in the cohort to produce is greater than or equal to zero
        assert(OffspringCohortAbundance >= 0.0 && "Offspring abundance < 0");

        // Get the adult and juvenile masses of the cohort to produce
        vector<double> OffspringProperties = GetOffspringCohortProperties(gridCellCohorts, actingCohort,
                madingleyCohortDefinitions);

        // Update cohort abundance in case juvenile mass has been altered
        OffspringCohortAbundance = (OffspringCohortAbundance * actingCohort.JuvenileMass) /
                OffspringProperties[0];

        //Create the offspring cohort
        Cohort OffspringCohort(actingCohort[0],
                OffspringProperties[0],
                OffspringProperties[1],
                OffspringProperties[0],
                OffspringCohortAbundance,
                currentTimestep, partial.NextCohortIDThreadLocked);

        // Add the offspring cohort to the grid cell cohorts array
        gridCellCohorts[actingCohort[0]].push_back(OffspringCohort);

        // Subtract all of the reproductive potential mass of the parent cohort, which has been used to generate the new
        // cohort, from the delta reproductive potential mass
        Cohort::Deltas["reproductivebiomass"]["reproduction"] -= (actingCohort.IndividualReproductivePotentialMass);

    }
    //----------------------------------------------------------------------------------------------  
    /** \brief    Assigns biomass from body mass to reproductive potential mass
    @param gridCellCohorts The cohorts in the current grid cell 
    @param gridCellStocks The stocks in the current grid cell 
    @param actingCohort The position of the acting cohort in the jagged array of grid cell cohorts 
    @param cellEnvironment The environment in the current grid cell 
    @param madingleyCohortDefinitions The definitions of cohort functional groups in the model 
    @param madingleyStockDefinitions The definitions of stock functional groups in the model 
    @param currentTimestep The current model time step 
    @param tracker An instance of ProcessTracker to hold diagnostics for reproduction */
    void AssignMassToReproductivePotential(GridCellCohortHandler gridCellCohorts, GridCellStockHandler gridCellStocks,
            vector<int> actingCohort, map<string, vector<double> > cellEnvironment, 
            FunctionalGroupDefinitions madingleyCohortDefinitions, FunctionalGroupDefinitions madingleyStockDefinitions,
            unsigned currentTimestep, ProcessTracker tracker) {
        // If this is the first time reproductive potential mass has been assigned for this cohort, 
        // then set the maturity time step for this cohort as the current model time step
        if (actingCohort.MaturityTimeStep == std::numeric_limits<unsigned>::max()) {
            actingCohort.MaturityTimeStep = currentTimestep;

        }

        // Assign the specified mass to reproductive potential mass and remove it from individual biomass
        Cohort::Deltas["reproductivebiomass"]["reproduction"] += BiomassToAssignToReproductivePotential;
        Cohort::Deltas["biomass"]["reproduction"] -= BiomassToAssignToReproductivePotential;

    }
    //----------------------------------------------------------------------------------------------  
};

#endif
