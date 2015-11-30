#ifndef TREVISEDREPRODUCTION_H
#define TREVISEDREPRODUCTION_H
#include <IReproductionimplementation.h>
#include <UtilityFunctions.h>
/** \file TRevisedReproduction.h
 * \brief the TRevisedReproduction header file
 */
///MB this class seems incomplete at present...
/** \brief A formulation of the process of reproduction */
class RevisedReproduction : public IReproductionImplementation {
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
    @param gcl The current grid cell
    @param actingCohort The acting cohort 
    @param currentTimestep The current model time step
    @param partial some thread locked variables?
    @param iteroparous breed mode of the cohort
    @param currentMonth as it says
    @param params the parameters 
 */
    void RunReproduction(GridCell& gcl, Cohort& actingCohort,
            unsigned currentTimestep,ThreadLockedParallelVariables& partial, 
            bool iteroparous, unsigned currentMonth, MadingleyModelInitialisation& params) {
        // Check that the abundance in the cohort to produce is greater than or equal to zero
        assert(OffspringCohortAbundance >= 0.0 && "Offspring abundance < 0");

        // Get the adult and juvenile masses of the cohort to produce
        vector<double> OffspringProperties = GetOffspringCohortProperties(gcl.GridCellCohorts, actingCohort,
                params.CohortFunctionalGroupDefinitions);

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
        gcl.GridCellCohorts[actingCohort.FunctionalGroupIndex].push_back(OffspringCohort);

        // Subtract all of the reproductive potential mass of the parent cohort, which has been used to generate the new
        // cohort, from the delta reproductive potential mass
        Cohort::Deltas["reproductivebiomass"]["reproduction"] -= (actingCohort.IndividualReproductivePotentialMass);

    }
    //----------------------------------------------------------------------------------------------  
    /** \brief    Assigns biomass from body mass to reproductive potential mass
    @param gridCell The current cell 
    @param actingCohort The acting cohort  
    @param currentTimestep The current model time step  */
    void AssignMassToReproductivePotential(GridCell& gcl,
            Cohort& actingCohort, unsigned currentTimestep) {
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
