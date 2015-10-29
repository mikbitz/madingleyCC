#ifndef IDISPERSALIMPLEMENTATION_H
#define IDISPERSALIMPLEMENTATION_H
#include <ModelGrid.h>
#include <Cohort.h>
#include <random>

/** \file IDispersalImplementation.h
 * \brief the IDispersalImplementation header file
 */

/** \brief Base Class for implementations of the ecological process of dispersal */
class IDispersalImplementation {
public:
    //----------------------------------------------------------------------------------------------
    //Variables
    //----------------------------------------------------------------------------------------------
    /** \brief Include Utility class */
    UtilityFunctions Utilities;
    /** \brief An instance of the simple random number generator class */
    std::default_random_engine RandomNumberGenerator;
    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------

    //----------------------------------------------------------------------------------------------
     /** \brief Run the dispersal implementation */
    virtual bool RunDispersal(vector<Cohort>& disperseMonkeys, vector<unsigned>& cellIndex, ModelGrid& gridForDispersal, Cohort& cohortToDisperse, int actingCohortFunctionalGroup, int actingCohortNumber, unsigned currentMonth) {
        ;
    }
    //----------------------------------------------------------------------------------------------
    void relocate(vector<Cohort>&disperseMonkeys, Cohort& cohortToDisperse,  const vector<unsigned>& destination) {
            cohortToDisperse.destination[0] = destination[0];
            cohortToDisperse.destination[1] = destination[1];
            disperseMonkeys.push_back(cohortToDisperse);
    }
    //----------------------------------------------------------------------------------------------
    vector<unsigned> newCell(ModelGrid& madingleyGrid,double& uSpeed,double& vSpeed,double & LatCellLength,double & LonCellLength,const unsigned& latIndex,const unsigned& lonIndex){
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
        assert(DispersalProbability <= 1 && "Dispersal probability should always be <= 1");

       vector<unsigned> DestinationCell= {9999999, 9999999};

        // Check to see in which axis the cohort disperses

        // Note that the values in the dispersal array are the proportional area moved outside the grid cell in each direction; we simply compare the random draw to this
        // to determine the direction in which the cohort moves probabilistically
        std::uniform_real_distribution<double> randomNumber(0.0, 1.0);
        double RandomValue = randomNumber(RandomNumberGenerator);
        if (DispersalProbability >= RandomValue) {
            int signu = (uSpeed > 0) - (uSpeed < 0);
            int signv = (vSpeed > 0) - (vSpeed < 0);
            // Longitudinally
            if (RandomValue <= AreaOutsideU / CellArea) {
                signv = 0;
            } else {
                //Latitudinally
                if (RandomValue <= (AreaOutsideU / CellArea + AreaOutsideV / CellArea)) {
                    signu = 0;
                }
            }
            // try to get a cell - only use it if realm is the same as the origin cell.
            
            vector<unsigned> FreshCell = madingleyGrid.getNewCell(latIndex, lonIndex, signv, signu);
            if ((FreshCell[0]<9999999) &&
                (madingleyGrid.InternalGrid[latIndex    ][lonIndex    ].CellEnvironment["Realm"][0] == 
                 madingleyGrid.InternalGrid[FreshCell[0]][FreshCell[1]].CellEnvironment["Realm"][0]) )
                    DestinationCell=FreshCell;
        }

        return DestinationCell;
    }
    //----------------------------------------------------------------------------------------------
};
#endif
