#ifndef TRESPONSIVEDISPERSAL_H
#define TRESPONSIVEDISPERSAL_H
#include <IDispersalImplementation.h>
#include <UtilityFunctions.h>
#include <random>
#include <chrono>
#include <assert.h>
#include <math.h>
#include <vector>
/** \file TResponsiveDispersal.h
 * \brief the TResponsiveDispersal header file
 */

/** \brief A formulation of the process of responsive dispersal */
class ResponsiveDispersal : public IDispersalImplementation {
    //----------------------------------------------------------------------------------------------
    //Variables
    //----------------------------------------------------------------------------------------------
    /** \brief The time units associated with this implementation of dispersal */
    const string TimeUnitImplementation = "month";

    /** \brief Density threshold below which adult individuals may move to look for other adults of the same cohort

    \remarks The density scales with cohort weight via: Min Density = DensityThresholdScaling / MatureMass (g) */
    const double DensityThresholdScaling = 50000;
    /** \brief Scalar relating dispersal speed to individual body mass */
    const double DispersalSpeedBodyMassScalar = 0.0278;
    /** \brief Body-mass exponent of the relationship between disperal speed and individual body mass */
    const double DispersalSpeedBodyMassExponent = 0.48;

    /** \brief The proportion of body mass loss at which the cohort will try to disperse every time during a time step */
    const double StarvationDispersalBodyMassThreshold = 0.8;

    /** \brief Scalar to convert from the time step units used by this formulation of dispersal to global model time step units */
    double DeltaT;

public:
    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------

    //----------------------------------------------------------------------------------------------
    /** \brief Assigns all parameter values for repsonsive dispersal */
    ResponsiveDispersal(string globalModelTimeStepUnit, bool DrawRandomly) {

        // Calculate the scalar to convert from the time step units used by this implementation of dispersal to the global model time step units
        DeltaT = Utilities.ConvertTimeUnits(globalModelTimeStepUnit, TimeUnitImplementation);

        // Set the seed for the random number generator
        if (DrawRandomly) {
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            RandomNumberGenerator.seed(seed);
        } else {
            RandomNumberGenerator.seed(14141);
        }
    }

    //----------------------------------------------------------------------------------------------
    /** \brief Run responsive dispersal
    @param cellIndices The longitudinal and latitudinal indices of the current grid cell 
    @param gridForDispersal The model grid to run dispersal for 
    @param cohortToDisperse The cohort for which to apply the dispersal process 
    @param actingCohortFunctionalGroup The functional group index of the acting cohort 
    @param actingCohortNumber The position of the acting cohort within the functional group in the array of grid cell cohorts 
    @param currentMonth The current model month 
     */
    bool RunDispersal(vector<Cohort>& dispersers, ModelGrid& gridForDispersal, Cohort& cohortToDisperse,
              const unsigned& currentMonth) {
        // Starvation driven dispersal takes precedence over density driven dispersal (i.e. a cohort can't do both). Also, the delta 
        // arrays only allow each cohort to perform one type of dispersal each time step
        bool CohortDispersed = false;

        // Check for starvation-driven dispersal
        CohortDispersed = CheckStarvationDispersal(dispersers,gridForDispersal, cohortToDisperse);

        if (!CohortDispersed) {
            // Check for density driven dispersal
            CheckDensityDrivenDispersal(dispersers,gridForDispersal, cohortToDisperse);
        }
        return false;
    }

    //----------------------------------------------------------------------------------------------
    bool CheckStarvationDispersal(vector<Cohort>& dispersers,ModelGrid& gridForDispersal, Cohort& cohortToDisperse) {
        // A boolean to check whether a cohort has dispersed
        bool CohortHasDispersed = false;

        // Check for starvation driven dispersal
        // What is the present body mass of the cohort?
        // Note that at present we are just tracking starvation for adults
        double IndividualBodyMass = cohortToDisperse.IndividualBodyMass;
        double AdultMass = cohortToDisperse.AdultMass;

        // Assume a linear relationship between probability of dispersal and body mass loss, up to _StarvationDispersalBodyMassThreshold
        // at which point the cohort will try to disperse every time step
        if (IndividualBodyMass < AdultMass) {
            double ProportionalPresentMass = IndividualBodyMass / AdultMass;

            // If the body mass loss is greater than the starvation dispersal body mass threshold, then the cohort tries to disperse
            if (ProportionalPresentMass < StarvationDispersalBodyMassThreshold) {
                // Cohort tries to disperse

                CalculateDispersalProbability(gridForDispersal, cohortToDisperse, CalculateDispersalSpeed(AdultMass));

                // Update the cell to disperse to, if the cohort moves
                if (cohortToDisperse.origin != cohortToDisperse.destination) {
                    // Update the delta array of cohorts
                    dispersers.push_back(cohortToDisperse);
                }

                // Note that regardless of whether or not it succeeds, if a cohort tries to disperse, it is counted as having dispersed for the purposes of not then allowing it to disperse
                // based on its density.
                CohortHasDispersed = true;
                // Otherwise, the cohort has a chance of trying to disperse proportional to its mass lass
            } else {          
            
                // Cohort tries to disperse with a particular probability
                // Draw a random number
                std::uniform_real_distribution<double> randomNumber(0.0, 1.0);
                double RandomValue = randomNumber(RandomNumberGenerator);
                if (((1.0 - ProportionalPresentMass) / (1.0 - StarvationDispersalBodyMassThreshold)) > RandomValue) {

                    CalculateDispersalProbability(gridForDispersal, cohortToDisperse, CalculateDispersalSpeed(AdultMass));

                // Update the cell to disperse to, if the cohort moves
                if (cohortToDisperse.origin != cohortToDisperse.destination) {
                    // Update the delta array of cohorts
                    dispersers.push_back(cohortToDisperse);
                }


                    CohortHasDispersed = true;
                }
            }

        }
        return CohortHasDispersed;
    }
    //----------------------------------------------------------------------------------------------

    void CheckDensityDrivenDispersal(vector<Cohort>& dispersers, ModelGrid& gridForDispersal, Cohort& cohortToDisperse) {
        // Check the population density
        double NumberOfIndividuals = cohortToDisperse.CohortAbundance;

        // Get the cell area, in kilometres squared
        double CellArea = cohortToDisperse.origin->CellArea();

        // If below the density threshold
        if ((NumberOfIndividuals / CellArea) < DensityThresholdScaling / cohortToDisperse.AdultMass) {
            // Check to see if it disperses (based on the same movement scaling as used in diffusive movement)
            // Calculate dispersal speed for that cohort
            double DispersalSpeed = CalculateDispersalSpeed(cohortToDisperse.AdultMass);

            // Cohort tries to disperse
            CalculateDispersalProbability(gridForDispersal, cohortToDisperse, DispersalSpeed);

            // Update the cell to disperse to, if the cohort moves
            if (cohortToDisperse.origin != cohortToDisperse.destination) {
                // Update the delta array of cohorts
                dispersers.push_back(cohortToDisperse);
            }
        }
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Calculate the average diffusive dispersal speed of individuals in a cohort given their body mass
    @param bodyMass The current body mass of an individual in the cohort 
    @return The (average) dispersal speed in kilometres per month*/
    double CalculateDispersalSpeed(double bodyMass) {
        return DispersalSpeedBodyMassScalar * pow(bodyMass, DispersalSpeedBodyMassExponent);
    }
    //----------------------------------------------------------------------------------------------
    /** \brief   Calculates the probability of responsive dispersal given average individual dispersal speed and grid cell
    @param madingleyGrid The model grid 
    @param latIndex The latitude index of the grid cell to check for dispersal 
    @param lonIndex The longitude index of the grid cell to check for dispersal 
    @param dispersalSpeed The average dispersal speed of individuals in the acting cohort 

    Note that the second, third, and fourth elements are always positive; thus, they do not indicate 'direction' in terms of dispersal.*/
    void CalculateDispersalProbability(ModelGrid& madingleyGrid, Cohort& c, double dispersalSpeed) {
        double LatCellLength = madingleyGrid.CellHeightsKm[c.origin->LatIndex()];
        double LonCellLength = madingleyGrid.CellWidthsKm[c.origin->LatIndex()];

        // Pick a direction at random
        std::uniform_real_distribution<double> randomNumber(0.0, 1.0);
        double RandomDirection = randomNumber(RandomNumberGenerator)* 2 * acos(-1.);

        // Calculate the u and v components given the dispersal speed
        double uSpeed = dispersalSpeed * cos(RandomDirection);
        double vSpeed = dispersalSpeed * sin(RandomDirection);

        // Check that the whole cell hasn't moved out (i.e. that dispersal speed is not greater than cell length). 
        // This could happen if dispersal speed was high enough; indicates a need to adjust the time step, or to slow dispersal
        if (uSpeed > LonCellLength)cout<<"Dispersal Big U "<< uSpeed<<endl;
        if (vSpeed > LatCellLength)cout<<"Dispersal Big V "<< vSpeed<<endl;

        //assert(((uSpeed > LonCellLength) || (vSpeed > LatCellLength)) && "Dispersal probability should always be <= 1");

        GridCell* destination=newCell(madingleyGrid,uSpeed,vSpeed,LatCellLength,LonCellLength,c.origin);
        c.TryLivingAt(destination);

    }
    
};
#endif
