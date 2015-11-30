#ifndef MODELGRID_H
#define MODELGRID_H
using namespace std;
#include <GridCell.h>
#include <map>
#include <limits>
#include <string>
class MadingleyModelInitialisation;
/** \brief A class containing the model grid (composed of individual grid cells) along with grid attributes.
           The model grid is referenced by [Lat index, Lon index]
 */
class ModelGrid {
public:
    //----------------------------------------------------------------------------------------------
    //Variables
    //----------------------------------------------------------------------------------------------
    
    // Variable to make sure that not more than one grid is instantiated
    unsigned NumGrids = 0;
    /** \brief   minimum latitude of the grid */
    float MinLatitude;
    /** \brief leftmost longitude of the leftmost cell of the grid */
    float MinLongitude;
    /** \brief lowest latitude of the highest cell in the grid */
    float MaxLatitude;
    /** \brief leftmost longitude of the rightmost cell in the grid */
    float MaxLongitude;
    /** \brief latitudinal length of each grid cell. Currently assumes all cells are equal sized. */
    float LatCellSize;
    /** \brief the longitudinal length of each grid cell. Currently assumes all cells are equal sized.  */
    float LonCellSize;
    /** \brief The number of latitudinal cells in the model grid */
    long NumLatCells;
    /** \brief The number of longitudinal cells in the model grid */
    long NumLonCells;
    /** \brief The bottom (southern-most) latitude of each row of grid cells */
    vector<float> Lats;
    /** \brief The left (western-most) longitude of each column of grid cells */
    vector<float> Lons;
    /** \brief Array of grid cells */
    vector< vector <GridCell> > InternalGrid;

    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------
     //empty placeholder constructor just to get started.

    ModelGrid(){    ;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Constructor for model grid to construct the grid for specific locations
    @param minLat Minimum grid latitude (degrees) 
    @param minLon Minimum grid longitude (degrees, currently -180 to 180) 
    @param maxLat Maximum grid latitude (degrees) 
    @param maxLon Maximum grid longitude (degrees, currently -180 to 180) 
    @param latCellSize Latitudinal size of grid cells 
    @param lonCellSize Longitudinal size of grid cells 
     */
    void SetUpGrid(float minLat, float minLon, float maxLat, float maxLon, float latCellSize, float lonCellSize) {
        // Add one to the counter of the number of grids. If there is more than one model grid, exit the program with a debug crash.
        NumGrids = NumGrids + 1;
        assert(NumGrids < 2 && "You have initialised more than one grid on which to apply models. At present, this is not supported");
        cout << "Initialising grid cell environment:" << endl;

        // CURRENTLY DEFINING MODEL CELLS BY BOTTOM LEFT CORNER
        MinLatitude = minLat;
        MinLongitude = minLon;
        MaxLatitude = maxLat;
        MaxLongitude = maxLon;
        LatCellSize = latCellSize;
        LonCellSize = lonCellSize;


        // Check to see if the number of grid cells is an integer
        double u = (MaxLatitude - MinLatitude) / LatCellSize - floor((MaxLatitude - MinLatitude) / LatCellSize);
        assert((u <= 1.e-15) && "Error: number of grid cells is non-integer: check cell size");

        NumLatCells = (long) ((MaxLatitude - MinLatitude) / LatCellSize);
        NumLonCells = (long) ((MaxLongitude - MinLongitude) / LonCellSize);
        Lats.resize(NumLatCells);
        Lons.resize(NumLonCells);

        // Set up latitude and longitude vectors - lower left
        for (int ii = 0; ii < NumLatCells; ii++) {
            Lats[ii] = MinLatitude + ii * LatCellSize;
        }
        for (int jj = 0; jj < NumLonCells; jj++) {
            Lons[jj] = MinLongitude + jj * LonCellSize;
        }

        // Instantiate a grid of grid cells
        InternalGrid.resize(NumLatCells);
        for (auto &g : InternalGrid)g.resize(NumLonCells);


        setCellCoords();
        cout << "\n" << endl;

    }
    //----------------------------------------------------------------------------------------------
    void setCellCoords() {
        // Loop over cells to set up the model grid
        for (unsigned ii = 0; ii < InternalGrid.size(); ii++) {
            for (unsigned jj = 0; jj < InternalGrid[ii].size(); jj++) {
                // Create the grid cell at the specified position
                InternalGrid[ii][jj].setCellCoords(Lats[ii], ii,
                        Lons[jj], jj, LatCellSize, LonCellSize, Environment::MissingValue);
            }
        }
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Remove an individual cohort from a functional group; necessary due to dispersal moving cohorts from one cell to another
     * NB work on a copy here not a reference as the latter would store the next cell ref. into c in the calling function
    @param functionalGroup Cohort functional group 
     */
    void DeleteGridCellIndividualCohort(Cohort c) {
        vector<Cohort>& z=c.location->GridCellCohorts[c.FunctionalGroupIndex];
        auto h=find_if(z.begin(),z.end(),[c](Cohort& k){return c.ID==k.ID;});
        if (c.ID != (*h).ID)cout<<"Strange things happening in grid delete? "<<c.ID<<" "<< (*h).ID<<endl;
        z.erase(h);
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Add a new cohort to an existing list of cohorts in the grid cell
    @param c The cohort object to add 
     */
    void AddNewCohortToGridCell(Cohort c) {
        c.location=c.destination;
        c.destination->GridCellCohorts[c.FunctionalGroupIndex].push_back(c);
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Move a new cohort to another grid cell
    @param c The cohort object to move 
     */
    void Move(Cohort c) {
        vector<Cohort>& z=c.location->GridCellCohorts[c.FunctionalGroupIndex];
        auto h=find_if(z.begin(),z.end(),[c](Cohort& k){return c.ID==k.ID;});
        if (c.ID != (*h).ID)cout<<"Strange things happening in grid move? "<<c.ID<<" "<< (*h).ID<<endl;
        z.erase(h);
        c.location=c.destination;
        c.destination->GridCellCohorts[c.FunctionalGroupIndex].push_back(c);
    }
     //----------------------------------------------------------------------------------------------
/** \brief  Get pointer to a viable cell to move to
    @param  gcl Pointer to the focal grid cell
    @param  v latitudinal displacement
    @param  u longitudinal displacement
    @return Pointer to cell that lies at displacement u,v from the current cell
    @remark Currently assumes wrapping in longitude
     */
    GridCell* getNewCell(GridCell* gcl, const int& v, const int& u) {
        GridCell* Cell = 0;
        unsigned latCell=gcl->LatIndex(), lonCell=gcl->LonIndex();
        if (latCell + v >= 0 && latCell + v < NumLatCells) {
            int lnc = lonCell + u;
            while (lnc < 0)lnc += NumLonCells;
            while (lnc >= NumLonCells)lnc -= NumLonCells;
                Cell=&(InternalGrid[latCell + v][lnc]);
        }
        return Cell;
    }
    //----------------------------------------------------------------------------------------------
    //Apply any function that operates on a grid cell to all cells in the grid
    template <typename F>
    void ask(F f) {
        for (int j = 0; j < NumLonCells; j++) {
            for (int i = 0; i < NumLatCells; i++) {
                f(InternalGrid[i][j]);
            }
        }
    }
};
#endif