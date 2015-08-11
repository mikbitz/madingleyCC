#ifndef TRESPONSIVEDISPERSAL_H
#define TRESPONSIVEDISPERSAL_H
#include <IDispersalImplementation.h>
#include <UtilityFunctions.h>
#include <random>
#include <chrono>
#include <assert.h>
#include <math.h>
/** \file TResponsiveDispersal.h
 * \brief the TResponsiveDispersal header file
 */


//namespace Madingley
//{
/** \brief A formulation of the process of responsive dispersal */
class ResponsiveDispersal : public IDispersalImplementation
    {

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
           DoubleProperty DeltaT;

/** \brief Include Utility class */
           UtilityFunctions Utilities;
/** \brief An instance of the simple random number generator class */
           std::default_random_engine RandomNumberGenerator;
public:
/** \brief Assigns all parameter values for repsonsive dispersal */
        ResponsiveDispersal(string globalModelTimeStepUnit, bool DrawRandomly)
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

/** \brief Run responsive dispersal
@param cellIndices The longitudinal and latitudinal indices of the current grid cell 
@param gridForDispersal The model grid to run dispersal for 
@param cohortToDisperse The cohort for which to apply the dispersal process 
@param actingCohortFunctionalGroup The functional group index of the acting cohort 
@param actingCohortNumber The position of the acting cohort within the functional group in the array of grid cell cohorts 
@param currentMonth The current model month 
*/
        void RunDispersal(vector<unsigned> cellIndices, ModelGrid gridForDispersal, Cohort cohortToDisperse, 
           int actingCohortFunctionalGroup, int actingCohortNumber, unsigned currentMonth)
       {
           // Starvation driven dispersal takes precedence over density driven dispersal (i.e. a cohort can't do both). Also, the delta 
           // arrays only allow each cohort to perform one type of dispersal each time step
           bool CohortDispersed = false;

           // Check for starvation-driven dispersal
           CohortDispersed = CheckStarvationDispersal(gridForDispersal, cellIndices[0], cellIndices[1], cohortToDisperse, actingCohortFunctionalGroup, actingCohortNumber);

           if (!CohortDispersed)
           {
               // Check for density driven dispersal
               CheckDensityDrivenDispersal(gridForDispersal, cellIndices[0], cellIndices[1], cohortToDisperse, actingCohortFunctionalGroup, actingCohortNumber);
           }
        }
bool CheckStarvationDispersal(ModelGrid gridForDispersal, unsigned latIndex, unsigned lonIndex, Cohort cohortToDisperse, int functionalGroup, int cohortNumber)
        {
           // A boolean to check whether a cohort has dispersed
           bool CohortHasDispersed = false;

           // Check for starvation driven dispersal
           // What is the present body mass of the cohort?
           // Note that at present we are just tracking starvation for adults
           double IndividualBodyMass = cohortToDisperse.IndividualBodyMass();
           double AdultMass = cohortToDisperse.AdultMass;

           // Assume a linear relationship between probability of dispersal and body mass loss, up to _StarvationDispersalBodyMassThreshold
           // at which point the cohort will try to disperse every time step
           if (IndividualBodyMass < AdultMass)
           {
               double ProportionalPresentMass = IndividualBodyMass / AdultMass;

               // If the body mass loss is greater than the starvation dispersal body mass threshold, then the cohort tries to disperse
               if (ProportionalPresentMass < StarvationDispersalBodyMassThreshold)
               {
                   // Cohort tries to disperse
                   vector<double> DispersalArray = CalculateDispersalProbability(gridForDispersal, latIndex, lonIndex, CalculateDispersalSpeed(AdultMass));
                   double CohortDispersed = CheckForDispersal(DispersalArray[0]);
                   if (CohortDispersed > 0)
                   {
                       vector<unsigned> DestinationCell = CellToDisperseTo(gridForDispersal, latIndex, lonIndex, DispersalArray, DispersalArray[0], DispersalArray[4], DispersalArray[5]);
                       
                       // Update the delta array of cells to disperse to, if the cohort moves
                       if (DestinationCell[0] < 999999)
                       {
                           // Update the delta array of cohorts
                           gridForDispersal.DeltaFunctionalGroupDispersalArray[latIndex][lonIndex].push_back((unsigned)functionalGroup);
                           gridForDispersal.DeltaCohortNumberDispersalArray[latIndex][lonIndex].push_back((unsigned)cohortNumber);
                       
                           // Update the delta array of cells to disperse to
                           gridForDispersal.DeltaCellToDisperseToArray[latIndex][lonIndex].push_back(DestinationCell);
                       }
                   }

                   // Note that regardless of whether or not it succeeds, if a cohort tries to disperse, it is counted as having dispersed for the purposes of not then allowing it to disperse
                   // based on its density.
                   CohortHasDispersed = true;
               }

               // Otherwise, the cohort has a chance of trying to disperse proportional to its mass lass
               else
               {
                   // Cohort tries to disperse with a particular probability
                   // Draw a random number
                   std::uniform_real_distribution<double> randomNumber(0.0,1.0);
                   double RandomValue=randomNumber(RandomNumberGenerator);
                   if (((1.0 - ProportionalPresentMass) / (1.0 - StarvationDispersalBodyMassThreshold)) > RandomValue)
                   {
                       // Cohort tries to disperse
                       vector<double> DispersalArray = CalculateDispersalProbability(gridForDispersal, latIndex, lonIndex, CalculateDispersalSpeed(AdultMass));
                       double CohortDispersed = CheckForDispersal(DispersalArray[0]);
                       if (CohortDispersed > 0)
                       {
                           vector<unsigned> DestinationCell = CellToDisperseTo(gridForDispersal, latIndex, lonIndex, DispersalArray, DispersalArray[0], DispersalArray[4], DispersalArray[5]);

                           // Update the delta array of cells to disperse to, if the cohort moves
                           if (DestinationCell[0] < 999999)
                           {
                               // Update the delta array of cohorts
                               gridForDispersal.DeltaFunctionalGroupDispersalArray[latIndex][lonIndex].push_back((unsigned)functionalGroup);
                               gridForDispersal.DeltaCohortNumberDispersalArray[latIndex][lonIndex].push_back((unsigned)cohortNumber);

                               // Update the delta array of cells to disperse to
                               gridForDispersal.DeltaCellToDisperseToArray[latIndex][lonIndex].push_back(DestinationCell);
                           }
                       }

                       CohortHasDispersed = true;
                   }
               }

           }
           return CohortHasDispersed;
        }

void CheckDensityDrivenDispersal(ModelGrid gridForDispersal, unsigned latIndex, unsigned lonIndex, Cohort cohortToDisperse, int functionalGroup, int cohortNumber)
        {
           // Check the population density
           double NumberOfIndividuals = cohortToDisperse.CohortAbundance;

           // Get the cell area, in kilometres squared
           double CellArea = gridForDispersal.GetCellEnvironment(latIndex, lonIndex)["Cell Area"][0];

           // If below the density threshold
           if ((NumberOfIndividuals / CellArea) < DensityThresholdScaling / cohortToDisperse.AdultMass)
           {
               // Check to see if it disperses (based on the same movement scaling as used in diffusive movement)
               // Calculate dispersal speed for that cohort
               double DispersalSpeed = CalculateDispersalSpeed(cohortToDisperse.IndividualBodyMass());

               // Cohort tries to disperse
               vector<double> DispersalArray = CalculateDispersalProbability(gridForDispersal, latIndex, lonIndex, CalculateDispersalSpeed(cohortToDisperse.AdultMass));
               
               double CohortDispersed = CheckForDispersal(DispersalArray[0]);
               
               if (CohortDispersed > 0)
               {
                   vector<unsigned> DestinationCell = CellToDisperseTo(gridForDispersal, latIndex, lonIndex, DispersalArray, DispersalArray[0], DispersalArray[4], DispersalArray[5]);

                   // Update the delta array of cells to disperse to, if the cohort moves
                   if (DestinationCell[0] < 999999)
                   {
                       // Update the delta array of cohorts
                       gridForDispersal.DeltaFunctionalGroupDispersalArray[latIndex][lonIndex].push_back((unsigned)functionalGroup);
                       gridForDispersal.DeltaCohortNumberDispersalArray[latIndex][lonIndex].push_back((unsigned)cohortNumber);

                       // Update the delta array of cells to disperse to
                       gridForDispersal.DeltaCellToDisperseToArray[latIndex][lonIndex].push_back(DestinationCell);
                   }
               }
           }
        }
//
/** \brief Calculate the average diffusive dispersal speed of individuals in a cohort given their body mass
@param bodyMass The current body mass of an individual in the cohort 
@return The (average) dispersal speed in kilometres per month*/
double CalculateDispersalSpeed(double bodyMass)
       {
           return DispersalSpeedBodyMassScalar * pow(bodyMass, DispersalSpeedBodyMassExponent);
       }

/** \brief
Calculates the probability of responsive dispersal given average individual dispersal speed and grid cell

@param madingleyGrid The model grid 
@param latIndex The latitude index of the grid cell to check for dispersal 
@param lonIndex The longitude index of the grid cell to check for dispersal 
@param dispersalSpeed The average dispersal speed of individuals in the acting cohort 
@return A six element array. 
The first element is the probability of dispersal.
The second element is the probability of dispersing in the u (longitudinal) direction
The third element is the probability of dispersing in the v (latitudinal) direction
The fourth element is the probability of dispersing in the diagonal direction
The fifth element is the u velocity
The sixth element is the v velocity
Note that the second, third, and fourth elements are always positive; thus, they do not indicate 'direction' in terms of dispersal.*/
vector<double> CalculateDispersalProbability(ModelGrid madingleyGrid, unsigned latIndex, unsigned lonIndex, double dispersalSpeed)
        {
           double LatCellLength = madingleyGrid.CellHeightsKm[latIndex];
           double LonCellLength = madingleyGrid.CellWidthsKm[latIndex];

           // Pick a direction at random
           std::uniform_real_distribution<double> randomNumber(0.0,1.0);
           double RandomDirection = randomNumber(RandomNumberGenerator)* 2 * acos(-1.);

           // Calculate the u and v components given the dispersal speed
           double uSpeed = dispersalSpeed * cos(RandomDirection);
           double vSpeed = dispersalSpeed * sin(RandomDirection);

           // Check that the whole cell hasn't moved out (i.e. that dispersal speed is not greater than cell length). 
           // This could happen if dispersal speed was high enough; indicates a need to adjust the time step, or to slow dispersal
           assert(((uSpeed > LonCellLength) || (vSpeed > LatCellLength)) && "Dispersal probability should always be <= 1");
           

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

           // Check that we don't have any issues
           assert (DispersalProbability >= 1 && "Dispersal probability should always be <= 1");

           vector<double> NewArray = { DispersalProbability, AreaOutsideU / CellArea, AreaOutsideV / CellArea, AreaOutsideBoth / CellArea, uSpeed, vSpeed };

           return NewArray;
        }
//
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
//        // Here we make the assumption that if the cell in the direction chosen is unsuitable, that the dispersal does not happen 
//        // (i.e. that the cohorts 'bumps up' against unsuitable habitat
vector<unsigned> CellToDisperseTo(ModelGrid madingleyGrid, unsigned latIndex, unsigned lonIndex, vector<double> dispersalArray, double RandomValue, double uSpeedIncDiffusion, double vSpeedIncDiffusion)
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
                   // Diagonally. Also allow for rounding errors here, otherwise we can get sent to the debug.fail incorrectly
                   // Only an issue here, because otherwise the random value is always larger than the sum of the elements of the dispersal array
                   if (RandomValue <= (dispersalArray[1] + dispersalArray[2] + dispersalArray[3]) + 0.0000000001)
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
                       assert("Error selecting cell to disperse to in TResponsivedispersal.h");
                       DestinationCell = { 9999999, 9999999 };
                   }
               }

           }
           return DestinationCell;
        }

    };
//}
#endif
