#ifndef TDIFFUSIVEDISPERSAL_H
#define TDIFFUSIVEDISPERSAL_H
#include <IDispersalImplementation.h>
#include <UtilityFunctions.h>
#include <random>
#include <chrono>
#include <assert.h>
#include <math.h>
/** \file TDiffusiveDispersal.h
 * \brief the TDiffusiveDispersal header file
 */






//
//namespace Madingley
//{
/** \brief A formulation of the process of dispersal */
class DiffusiveDispersal : public IDispersalImplementation
    {
 /** \brief The time units associated with this implementation of dispersal */
           const string TimeUnitImplementation = "month";
/** \brief Scalar relating dispersal speed to individual body mass */
           const double DispersalSpeedBodyMassScalar = 0.0278;
/** \brief Body-mass exponent of the relationship between disperal speed and individual body mass */
           const double DispersalSpeedBodyMassExponent = 0.48;              
/** \brief Scalar to convert from the time step units used by this formulation of dispersal to global model time step units */
           DoubleProperty DeltaT;

/** \brief Include Utility class */
           UtilityFunctions Utilities;
/** \brief An instance of the simple random number generator class */
           std::default_random_engine RandomNumberGenerator;

    public:
/** \brief
Constructor for dispersal: assigns all parameter values
*/
         DiffusiveDispersal(string globalModelTimeStepUnit, bool DrawRandomly)
        {

            // Calculate the scalar to convert from the time step units used by this implementation of dispersal to the global model time step units
            DeltaT = Utilities.ConvertTimeUnits(globalModelTimeStepUnit, TimeUnitImplementation);

           // Set the seed for the random number generator
           if (DrawRandomly)
           {
               unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
               RandomNumberGenerator.seed(seed);
           }
           else
           {
               RandomNumberGenerator.seed(14141);
           }
        }
//
/** \brief
Run diffusive dispersal

@param cellIndices List of indices of active cells in the model grid 
@param gridForDispersal The model grid to run dispersal for 
@param cohortToDisperse The cohort for which to run the dispersal process for 
@param actingCohortFunctionalGroup The functional group index of the acting cohort 
@param actingCohortNumber The position of the cohort within the functional group in the array of grid cell cohorts 
@param currentMonth The current model month */
         void RunDispersal(vector<unsigned>& cellIndices, ModelGrid& gridForDispersal, Cohort& cohortToDisperse, 
            int actingCohortFunctionalGroup, int actingCohortNumber, unsigned currentMonth)
        {
           // Calculate dispersal speed for the cohort         
           double DispersalSpeed = CalculateDispersalSpeed(cohortToDisperse.IndividualBodyMass());

           // A double to indicate whether or not the cohort has dispersed, and if it has dispersed, where to
           double CohortDispersed = 0;

           // Get the probability of dispersal
           vector<double> DispersalArray = CalculateDispersalProbability(gridForDispersal, cellIndices[0], cellIndices[1], DispersalSpeed);

           // Check to see if it does disperse
           CohortDispersed = CheckForDispersal(DispersalArray[0]);

           // If it does, check to see where it will end up
           if (CohortDispersed > 0)
           {
               // Check to see if the direction is actually dispersable
               vector<unsigned> DestinationCell = CellToDisperseTo(gridForDispersal, cellIndices[0], cellIndices[1], DispersalArray, CohortDispersed, DispersalArray[4], DispersalArray[5]);

               if (DestinationCell[0] < 999999)
               {
                   // Update the delta array of cohorts
                   gridForDispersal.DeltaFunctionalGroupDispersalArray[cellIndices[0]][cellIndices[1]].push_back((unsigned)actingCohortFunctionalGroup);
                   gridForDispersal.DeltaCohortNumberDispersalArray[cellIndices[0]][cellIndices[1]].push_back((unsigned)actingCohortNumber);

                   // Update the delta array of cells to disperse to
                   gridForDispersal.DeltaCellToDisperseToArray[cellIndices[0]][cellIndices[1]].push_back(DestinationCell);
               }
           }
        }
/** \brief Calculates the average diffusive dispersal speed of individuals in a cohort given their body mass
@param bodyMass The current body mass of individuals in the cohort 
@return The average dispersal speed, in km per month*/
double CalculateDispersalSpeed(double bodyMass)
       {
               return DispersalSpeedBodyMassScalar * pow(bodyMass,DispersalSpeedBodyMassExponent);
       }

/** \brief Calculates the probability of diffusive dispersal given average individual dispersal speed
@param madingleyGrid The model grid 
@param latIndex The latitude index of the grid cell to check for dispersal 
@param lonIndex The longitude index of the grid cell to check for dispersal 
@param dispersalSpeed The average speed at which individuals in this cohort move around their environment, in km per month 
@return A six element array. 
The first element is the probability of dispersal.
The second element is the probability of dispersing in the u (longitudinal) direction
The third element is the probability of dispersing in the v (latitudinal) direction
The fourth element is the probability of dispersing in the diagonal direction
The fifth element is the u velocity modified by the random diffusion component
The sixth element is the v velocity modified by the random diffusion component
Note that the second, third, and fourth elements are always positive; thus, they do not indicate 'direction' in terms of dispersal.*/
vector <double> CalculateDispersalProbability(ModelGrid& madingleyGrid, unsigned latIndex, unsigned lonIndex, double dispersalSpeed)
        {
           // Check that the u speed and v speed are not greater than the cell length. If they are, then rescale them; this limits the max velocity
           // so that cohorts cannot be advected more than one grid cell per time step
           double LatCellLength = madingleyGrid.CellHeightsKm[latIndex];
           double LonCellLength = madingleyGrid.CellWidthsKm[latIndex];

           // Pick a direction at random
           std::uniform_real_distribution<double> randomNumber(0.0,1.0);
           double RandomDirection = randomNumber(RandomNumberGenerator)* 2 * acos(-1.);


           // Calculate the u and v components given the dispersal speed
           double uSpeed = dispersalSpeed * cos(RandomDirection);
           double vSpeed = dispersalSpeed * sin(RandomDirection);
           
           // Calculate the area of the grid cell that is now outside in the diagonal direction
           double AreaOutsideBoth = abs(uSpeed * vSpeed);

           // Calculate the area of the grid cell that is now outside in the u direction (not including the diagonal)
           double AreaOutsideU = abs(uSpeed * LatCellLength) - AreaOutsideBoth;

           // Calculate the proportion of the grid cell that is outside in the v direction (not including the diagonal
           double AreaOutsideV = abs(vSpeed * LonCellLength) - AreaOutsideBoth;

           // Get the cell area, in kilometres squared
           double CellArea = madingleyGrid.GetCellEnvironment(latIndex, lonIndex)["Cell Area"][0];

           // Convert areas to a probability
           double DispersalProbability = (AreaOutsideU + AreaOutsideV + AreaOutsideBoth) / CellArea;
           
           // Check that the whole cell hasn't moved out. This could happen if dispersal speed was high enough

           assert(DispersalProbability <= 1 && "Dispersal probability in diffusion should always be <= 1");


           vector<double> NewArray = { DispersalProbability, AreaOutsideU / CellArea, AreaOutsideV / CellArea, AreaOutsideBoth / CellArea, uSpeed, vSpeed };

           return NewArray;
       }

/** \brief Generate a random value to see if a cohort disperses
@param dispersalProbability The probability of dispersal 
@return Returns either the random value, if it less than dispersal probability, or -1*/
double CheckForDispersal(double dispersalProbability)
        {
           // Randomly check to see if dispersal occurs
           std::uniform_real_distribution<double> randomNumber(0.0,1.0);
           double RandomValue = randomNumber(RandomNumberGenerator);
           if (dispersalProbability >= RandomValue)
           {
               return RandomValue;
           }
           else
           {
               return -1.0;
           }
        }
//
//        // Determine to which cell the cohort disperses
//        // Note that if the direction is not dispersable, then it doesn't happen
vector<unsigned> CellToDisperseTo(ModelGrid& madingleyGrid, unsigned latIndex, unsigned lonIndex, vector<double>& dispersalArray, double RandomValue, double uSpeedIncDiffusion, double vSpeedIncDiffusion)
        {
           vector<unsigned> DestinationCell;

           // Check to see in which axis the cohort disperses
           // Longitudinally
           if (RandomValue <= dispersalArray[1])
           {
               // Work out whether dispersal is to the cell to the E or the W
               if (uSpeedIncDiffusion > 0)
               {
                   DestinationCell = madingleyGrid.CheckDispersalEast(latIndex, lonIndex);
               }
               else
               {
                   DestinationCell = madingleyGrid.CheckDispersalWest(latIndex, lonIndex);
               }

           }
           else
           {
               // Latitudinally
               if (RandomValue <= (dispersalArray[1] + dispersalArray[2]))
               {
                   // Work out whether dispersal is to the cell to the N or the S
                   if (vSpeedIncDiffusion > 0)
                   {
                       DestinationCell = madingleyGrid.CheckDispersalNorth(latIndex, lonIndex);
                   }
                   else
                   {
                       DestinationCell = madingleyGrid.CheckDispersalSouth(latIndex, lonIndex);
                   }

               }
               else
               {
                   // Diagonally
                   if (RandomValue <= (dispersalArray[1] + dispersalArray[2] + dispersalArray[3]))
                   {
                       // Work out to which cell dispersal occurs
                       if (uSpeedIncDiffusion > 0)
                       {
                           if (vSpeedIncDiffusion > 0)
                           {
                               DestinationCell = madingleyGrid.CheckDispersalNorthEast(latIndex, lonIndex);
                           }
                           else
                           {
                               DestinationCell = madingleyGrid.CheckDispersalSouthEast(latIndex, lonIndex);
                           }

                       }
                       else
                       {
                           if (vSpeedIncDiffusion > 0)
                           {
                               DestinationCell = madingleyGrid.CheckDispersalNorthWest(latIndex, lonIndex);
                           }
                           else
                           {
                               DestinationCell = madingleyGrid.CheckDispersalSouthWest(latIndex, lonIndex);
                           }
                       }
                   }
                   else
                   {
                       assert("Error with advection when determining cell to disperse to in TDiffusiveDispersal.h");
                       DestinationCell = { 9999999, 9999999 };
                   }
               }

           }
           return DestinationCell;
        }
    };
//}
#endif
