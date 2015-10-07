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
    /** \brief An instance of the simple random number generator class */
    std::default_random_engine RandomNumberGenerator;

    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------
    
    //----------------------------------------------------------------------------------------------
    /** \brief Constructor for dispersal: assigns all parameter values */
    AdvectiveDispersal(string globalModelTimeStepUnit, bool DrawRandomly) {
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
    @param actingCohortNumber The position of the acting cohort wihtin the functional group in the array of grid cell cohorts 
    @param currentMonth The current model month 
     */
    void RunDispersal(vector<unsigned>& cellIndex, ModelGrid& gridForDispersal, Cohort& cohortToDisperse, int actingCohortFunctionalGroup,
            int actingCohortNumber, unsigned currentMonth) {
        // An array to hold the dispersal information
        vector<double> DispersalArray(6);

        // A double to indicate whether or not the cohort has dispersed, and if it has dispersed, where to
        double CohortDispersed = 0;

        // An array to hold the present cohort location for the intermediate steps that occur before the final dispersal this time step
        vector<unsigned> PresentLocation = {cellIndex[0], cellIndex[1]};

        // Loop through a number of times proportional to the rescaled dispersal
        for (int mm = 0; mm < AdvectionTimeStepsPerModelTimeStep; mm++) {
            // Get the probability of dispersal
            DispersalArray = CalculateDispersalProbability(gridForDispersal, PresentLocation[0], PresentLocation[1], currentMonth);

            // Check to see if it does disperse
            CohortDispersed = CheckForDispersal(DispersalArray[0]);

            // If it does, check to see where it will end up
            if (CohortDispersed > 0) {
                // Check to see if the direction is actually dispersable
                vector<unsigned> DestinationCell = CellToDisperseTo(gridForDispersal, PresentLocation[0], PresentLocation[1], DispersalArray, CohortDispersed, DispersalArray[4], DispersalArray[5]);

                // If it is, go ahead and update the cohort location
                if (DestinationCell[0] < 999999) {
                    PresentLocation = DestinationCell;
                }
            }
        }


        // Update the dipersal deltas for this cohort, if necessary
        if ((cellIndex[0] != PresentLocation[0]) || (cellIndex[1] != PresentLocation[1])) {
            // Update the delta array of cohorts - MB does this work?? where reset? every timestep??
            gridForDispersal.DeltaFunctionalGroupDispersalArray[cellIndex[0]][ cellIndex[1]].push_back((unsigned) actingCohortFunctionalGroup);
            gridForDispersal.DeltaCohortNumberDispersalArray[cellIndex[0]][cellIndex[1]].push_back((unsigned) actingCohortNumber);

            // Update the delta array of cells to disperse to
            gridForDispersal.DeltaCellToDisperseToArray[cellIndex[0]][cellIndex[1]].push_back(PresentLocation);
        }
    }
    //----------------------------------------------------------------------------------------------
    /** \brief    Convert dispersal speed from m per second to km per dispersal time step (currently 18h)
    @param dispersalSpeed The dispersal speed in m per second 
    @return The dispersal speed in kilometres per time step*/
    double RescaleDispersalSpeed(double dispersalSpeed) {
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
    vector<double> CalculateDispersalProbability(ModelGrid& madingleyGrid, unsigned latIndex, unsigned lonIndex, unsigned currentMonth) {
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

        // Area of the grid cell that is outside in the diagonal direction after dispersal, in kilometres squared
        double AreaOutsideBoth;

        // Area of the grid cell that is  outside in the u (longitudinal) direction after dispersal, in kilometres squared
        double AreaOutsideU;

        // Area of the grid cell that is  outside in the v (latitudinal) direction after dispersal, in kilometres squared
        double AreaOutsideV;

        // Cell area, in kilometres squared
        double CellArea;

        // Probability of dispersal
        double DispersalProbability;

        // A variable to track whether a named data layer exists
        bool varExists;

        // Get the u speed and the v speed from the cell data
        uAdvectiveSpeed = madingleyGrid.GetEnviroLayer("uVel", currentMonth, latIndex, lonIndex, varExists);
        assert(uAdvectiveSpeed != -9999);

        vAdvectiveSpeed = madingleyGrid.GetEnviroLayer("vVel", currentMonth, latIndex, lonIndex, varExists);
        assert(vAdvectiveSpeed != -9999);

        // These should be unnecessary if interpolation is working correctly and all grid cells have a velocity speed. If unsure, then uncomment.
        /*
        if (uSpeed == -9999)
        {
            uSpeed = 0;
        }
        if (vSpeed == -9999)
        {
            vSpeed = 0;
        }
         */

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

        // We assume that the whole grid cell moves at the given velocity and calculate the area that is then outside the original grid cell location. 
        // This then becomes the probability of dispersal

        // Calculate the area of the grid cell that is now outside in the diagonal direction. 
        AreaOutsideBoth = abs(uDistanceTravelled * vDistanceTravelled);

        // Calculate the area of the grid cell that is now outside in the u (longitudinal) direction (not including the diagonal)
        AreaOutsideU = abs(uDistanceTravelled * LatCellLength) - AreaOutsideBoth;

        // Calculate the proportion of the grid cell that is outside in the v (latitudinal) direction (not including the diagonal)
        AreaOutsideV = abs(vDistanceTravelled * LonCellLength) - AreaOutsideBoth;

        // Get the cell area, in kilometres squared
        CellArea = madingleyGrid.GetCellEnvironment(latIndex, lonIndex)["Cell Area"][0];

        // Convert areas to a probability
        DispersalProbability = (AreaOutsideU + AreaOutsideV + AreaOutsideBoth) / CellArea;

        // Check that the whole cell hasn't moved out. Could this happen for the fastest currents in a month? Definitely, 
        // if current speeds were not constrained
        assert(DispersalProbability <= 1 && "Dispersal probability in advection should always be <= 1");


        vector<double> NewArray = {DispersalProbability, AreaOutsideU / CellArea, AreaOutsideV / CellArea, AreaOutsideBoth / CellArea, uDistanceTravelled, vDistanceTravelled};
        return NewArray;
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
    /** \brief    Generate a random value to see if a cohort disperses
    @param dispersalProbability The probability of dispersal 
    @return Returns either the random value, if it less than dispersal probability, or -1*/
    double CheckForDispersal(double dispersalProbability) {
        // Randomly check to see if dispersal occurs
        std::uniform_real_distribution<double> randomNumber(0.0, 1.0);
        double RandomValue = randomNumber(RandomNumberGenerator);
        if (dispersalProbability >= RandomValue) {
            return RandomValue;
        } else {
            return -1.0;
        }


    }
             
    // Determine to which cell the cohort disperses

    vector<unsigned> CellToDisperseTo(ModelGrid& madingleyGrid, unsigned latIndex, unsigned lonIndex, vector<double>& dispersalArray, double RandomValue, double uSpeedIncDiffusion, double vSpeedIncDiffusion) {
        vector<unsigned> DestinationCell;

        // Check to see in which axis the cohort disperses

        // Note that the values in the dispersal array are the proportional area moved outside the grid cell in each direction; we simply compare the random draw to this
        // to determine the direction in which the cohort moves probabilistically

        // Longitudinally
        if (RandomValue <= dispersalArray[1]) {
            // Work out whether dispersal is to the cell to the E or the W
            if (uSpeedIncDiffusion > 0) {
                DestinationCell = madingleyGrid.CheckDispersalEast(latIndex, lonIndex);
            } else {
                DestinationCell = madingleyGrid.CheckDispersalWest(latIndex, lonIndex);
            }

        } else {
            // Latitudinally
            if (RandomValue <= (dispersalArray[1] + dispersalArray[2])) {
                // Work out whether dispersal is to the cell to the N or the S
                if (vSpeedIncDiffusion > 0) {
                    DestinationCell = madingleyGrid.CheckDispersalNorth(latIndex, lonIndex);
                } else {
                    DestinationCell = madingleyGrid.CheckDispersalSouth(latIndex, lonIndex);
                }

            } else {
                // Diagonally
                if (RandomValue <= (dispersalArray[1] + dispersalArray[2] + dispersalArray[3])) {
                    // Work out to which cell dispersal occurs
                    if (uSpeedIncDiffusion > 0) {
                        if (vSpeedIncDiffusion > 0) {
                            DestinationCell = madingleyGrid.CheckDispersalNorthEast(latIndex, lonIndex);
                        } else {
                            DestinationCell = madingleyGrid.CheckDispersalSouthEast(latIndex, lonIndex);
                        }

                    } else {
                        if (vSpeedIncDiffusion > 0) {
                            DestinationCell = madingleyGrid.CheckDispersalNorthWest(latIndex, lonIndex);
                        } else {
                            DestinationCell = madingleyGrid.CheckDispersalSouthWest(latIndex, lonIndex);
                        }
                    }
                } else {
                    //INSERT ERROR HANDLING CODE HERE("Error with advection when determining cell to disperse to");
                    // Should also be there when not running in debug mode
                    assert("Error when determining which cell to disperse to");
                    vector<unsigned>DestinationCell = {9999999, 9999999};
                }
            }

        }
        return DestinationCell;
    }
     //----------------------------------------------------------------------------------------------
};

#endif
