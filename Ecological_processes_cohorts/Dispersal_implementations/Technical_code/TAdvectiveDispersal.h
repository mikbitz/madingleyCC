#ifndef TADVECTIVEDISPERSAL_H
#define TADVECTIVEDISPERSAL_H
#include <IDispersalImplementation.h>
#include <UtilityFunctions.h>
#include <random>
#include <chrono>
#include <assert.h>
#include <math.h>
/** \file TAdvectiveDispersal.h
 * \brief the TAdvectiveDispersal header file
 */

/** \brief A formulation of the process of dispersal */
class AdvectiveDispersal : public IDispersalImplementation {
public:
    //----------------------------------------------------------------------------------------------
    //Variables
    //----------------------------------------------------------------------------------------------

    /** \brief The horizontal diffusivity parameter (m^2/s) */
    const double HorizontalDiffusivity = 100;
    /** \brief The length of the time-step for advective dispersal, in hours*/
    const unsigned AdvectiveModelTimeStepLengthHours = 18;
    /** \brief Horizontal diffusivity in km^2/advective-dispersal-time-step*/
    const double HorizontalDiffusivityKmSqPerADTimeStep = HorizontalDiffusivity / (1000 * 1000) * 60 * 60 * AdvectiveModelTimeStepLengthHours;
    /** \brief Time unit scalar to apply to advective dispersal*/
    double AdvectionTimeStepsPerModelTimeStep;
    /** \brief The time units associated with this implementation of dispersal*/
    const string TimeUnitImplementation = "month";
    /** \brief Factor to convert velocity from m/s to km/month*/
    double VelocityUnitConversion;
    /** \brief Scalar to convert from the time step units used by this formulation of dispersal to global model time step units */
    double DeltaT;
    /** \brief Include Utility class */
    UtilityFunctions Utilities;

    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------
    
    //----------------------------------------------------------------------------------------------
    /** \brief Constructor for dispersal: assigns all parameter values */
    AdvectiveDispersal(const string& globalModelTimeStepUnit, const bool& DrawRandomly) {
        // Calculate the scalar to convert from the time step units used by this implementation of dispersal to the global model time step units
        DeltaT = Utilities.ConvertTimeUnits(globalModelTimeStepUnit, TimeUnitImplementation);

        // Initialise the advective dispersal temporal scaling to adjust between time steps appropriately
        AdvectionTimeStepsPerModelTimeStep = Utilities.ConvertTimeUnits(globalModelTimeStepUnit, "day") * 24 / AdvectiveModelTimeStepLengthHours;

        // Convert velocity from m/s to km/month. Note that if the _TimeUnitImplementation changes, this will also have to change.
        VelocityUnitConversion = 60 * 60 * 24 * Utilities.ConvertTimeUnits(globalModelTimeStepUnit, "day") * DeltaT / 1000;

        // Set the seed for the random number generator
        if (DrawRandomly) {
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            RandomNumberGenerator.seed(seed);
        } else {
            RandomNumberGenerator.seed(14141);
        }

    }
    //----------------------------------------------------------------------------------------------
    /** \brief Run advective dispersal
    @param cellIndex The longitudinal and latitudinal indices of the focal grid cell 
    @param gridForDispersal The model grid to run dispersal for 
    @param cohortToDisperse The cohort to run dispersal for 
    @param actingCohortFunctionalGroup The functional group index of the acting cohort 
    @param actingCohortNumber The position of the acting cohort within the functional group in the array of grid cell cohorts 
    @param currentMonth The current model month 
     */
    bool RunDispersal(vector<Cohort>& disperseMonkeys,  ModelGrid& gridForDispersal, Cohort& cohortToDisperse , const unsigned& currentMonth) {
        // A double to indicate whether or not the cohort has dispersed, and if it has dispersed, where to
        double CohortDispersed = 0;

        // An array to hold the present cohort location for the intermediate steps that occur before the final dispersal this time step
        vector<unsigned> destination = {cohortToDisperse.origin[0], cohortToDisperse.origin[1]};

        vector<unsigned>currentCell=destination;
        // Loop through a number of times proportional to the rescaled dispersal
        for (int mm = 0; mm < AdvectionTimeStepsPerModelTimeStep; mm++) {
            // Get the probability of dispersal and return a candidate cell
            currentCell = CalculateDispersalProbability(gridForDispersal, destination[0], destination[1], currentMonth);

                if (currentCell[0] < 999999) {
                    destination = currentCell;
                }
        }
        // Update the dipersal deltas for this cohort, if necessary
        if ((cohortToDisperse.origin[0] != destination[0]) || (cohortToDisperse.origin[1] != destination[1])) {
                relocate(disperseMonkeys,cohortToDisperse,destination);//stores a copy of the cohort
            
            return true;
        }
        return false;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief    Convert dispersal speed from m per second to km per dispersal time step (currently 18h)
    @param dispersalSpeed The dispersal speed in m per second 
    @return The dispersal speed in kilometres per time step*/
    inline const double RescaleDispersalSpeed(const double& dispersalSpeed) const {
        //            // Units are metres per second; need to convert to kilometres per global time step (currently one month) - use VelocityUnitConversion for this.
        //            // Also rescale based on the time step of the advective dispersal model - currently 18h
        return dispersalSpeed * VelocityUnitConversion / AdvectionTimeStepsPerModelTimeStep;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief    Calculates the probability of advective dispersal given the grid cell
    @param madingleyGrid The model grid 
    @param latIndex The latitude index of the grid cell to check for dispersal 
    @param lonIndex The longitude index of the grid cell to check for dispersal 
    @param currentMonth The current model month 
    @return A six element array. 
    The first element is the probability of dispersal.
    The second element is the probability of dispersing in the u (longitudinal) direction
    The third element is the probability of dispersing in the v (latitudinal) direction
    The fourth element is the probability of dispersing in the diagonal direction
    The fifth element is the distance travelled in the u direction (u velocity modified by the random diffusion component)
    The sixth element is the distance travelled in the v direction (v velocity modified by the random diffusion component)
    Note that the second, third, and fourth elements are always positive; thus, they do not indicate 'direction' in terms of dispersal.*/
    const vector<unsigned> CalculateDispersalProbability(ModelGrid& madingleyGrid, const unsigned& latIndex,const unsigned& lonIndex,const unsigned& currentMonth) {
        // Advective speed in u (longitudinal) direction
        double uAdvectiveSpeed;

        // Advective speed in v (latitudinal) direction
        double vAdvectiveSpeed;

        // Distance travelled in u (longitudinal) direction
        double uDistanceTravelled;

        // Distance travelled in v (latitudinal) direction
        double vDistanceTravelled;

        // U and V components of the diffusive velocity
        vector<double> DiffusiveUandVComponents(2);

        // Length in km of a cell boundary latitudinally
        double LatCellLength;

        // Length in km of a cell boundary longitudinally
        double LonCellLength;

        // Cell area, in kilometres squared
        double CellArea;

        // A variable to track whether a named data layer exists
        bool varExists;

        // Get the u speed and the v speed from the cell data
        uAdvectiveSpeed = madingleyGrid.GetEnviroLayer("uVel", currentMonth, latIndex, lonIndex, varExists);
        assert(uAdvectiveSpeed != -9999);

        vAdvectiveSpeed = madingleyGrid.GetEnviroLayer("vVel", currentMonth, latIndex, lonIndex, varExists);
        assert(vAdvectiveSpeed != -9999);

        // Calculate the diffusive movement speed, with a direction chosen at random
        DiffusiveUandVComponents = CalculateDiffusion();

        // Calculate the distance travelled in this dispersal (not global) time step. both advective and diffusive speeds need to have been converted to km / advective model time step
        uDistanceTravelled = RescaleDispersalSpeed(uAdvectiveSpeed) + DiffusiveUandVComponents[0];
        vDistanceTravelled = RescaleDispersalSpeed(vAdvectiveSpeed) + DiffusiveUandVComponents[1];
        
        // Check that the u distance travelled and v distance travelled are not greater than the cell length

        LatCellLength = madingleyGrid.CellHeightsKm[latIndex];
        LonCellLength = madingleyGrid.CellWidthsKm[latIndex];

        assert(abs(uDistanceTravelled) <= LonCellLength && "u velocity greater than cell width");
        assert(abs(vDistanceTravelled) <= LatCellLength && "v velocity greater than cell width");

        return newCell(madingleyGrid,uDistanceTravelled,vDistanceTravelled,LatCellLength,LonCellLength,latIndex,lonIndex);

    }
    //----------------------------------------------------------------------------------------------
    /** \brief    Get a randomly directed diffusion vector. This is derived from the LTRANS model formulation, which itself is derived from Visser 1997 (MEPS)
    We assume that the standard deviation of the random draw is 1.0

    @return A two element array, where the first element is the diffusion component in the u direction, and the second component is the
    diffusion component in the v direction*/
    vector<double> CalculateDiffusion() {
        // Create the array with which to send the output
        vector<double> UandVOutputs(2);

        // Note that this formulation drops the delta t because we set the horizontal diffusivity to be at the same temporal
        // scale as the time step
        std::normal_distribution<double> randomNumber(0., 1.0);
        UandVOutputs[0] = randomNumber(RandomNumberGenerator) * sqrt((2.0 * HorizontalDiffusivityKmSqPerADTimeStep));
        UandVOutputs[1] = randomNumber(RandomNumberGenerator) * sqrt((2.0 * HorizontalDiffusivityKmSqPerADTimeStep));

        return UandVOutputs;
    }
    //----------------------------------------------------------------------------------------------
};

#endif
