#include <Cohort.h>
#include <limits.h>
#include <GridCell.h>

/** \file Cohort.cc
 * \brief the Cohort implementation file
 */
//----------------------------------------------------------------------------------------------

    Cohort::Cohort(GridCell& gcl, unsigned functionalGroupIndex, double juvenileBodyMass, double adultBodyMass, double initialBodyMass, double initialAbundance, double optimalPreyBodySizeRatio, unsigned short birthTimeStep, double proportionTimeActive, long long &nextCohortID) {

        FunctionalGroupIndex = functionalGroupIndex;
        JuvenileMass = juvenileBodyMass;
        AdultMass = adultBodyMass;
        IndividualBodyMass = initialBodyMass;
        CohortAbundance = initialAbundance;
        BirthTimeStep = birthTimeStep;
        MaturityTimeStep = std::numeric_limits<unsigned>::max();
        LogOptimalPreyBodySizeRatio = log(optimalPreyBodySizeRatio);
        MaximumAchievedBodyMass = juvenileBodyMass;
        Merged = false;
        ProportionTimeActive = proportionTimeActive;
        positionInList=gcl.GridCellCohorts[FunctionalGroupIndex].size();
        origin[0]=gcl.CellEnvironment["LatIndex"][0];
        origin[1]=gcl.CellEnvironment["LonIndex"][0];        
        destination[0]=origin[0];
        destination[1]=origin[1];
        ID=nextCohortID;//added MB to track this object.
        nextCohortID++;
    }
    //----------------------------------------------------------------------------------------------

    Cohort::Cohort(Cohort& actingCohort, unsigned p, double juvenileBodyMass, double adultBodyMass, double initialBodyMass, double initialAbundance, unsigned birthTimeStep, long long& nextCohortID) {

        FunctionalGroupIndex = actingCohort.FunctionalGroupIndex;
        JuvenileMass = juvenileBodyMass;
        AdultMass = adultBodyMass;
        IndividualBodyMass = initialBodyMass;
        CohortAbundance = initialAbundance;
        BirthTimeStep = birthTimeStep;
        MaturityTimeStep = std::numeric_limits<unsigned>::max();
        LogOptimalPreyBodySizeRatio = actingCohort.LogOptimalPreyBodySizeRatio;
        MaximumAchievedBodyMass = juvenileBodyMass;
        Merged = false;
        ProportionTimeActive = actingCohort.ProportionTimeActive;
        positionInList = p;
        origin[0] = actingCohort.origin[0];
        origin[1] = actingCohort.origin[1];
        destination[0] = origin[0];
        destination[1] = origin[1];
        ID = nextCohortID; //added MB to track this object.
        nextCohortID++;
    }
    //----------------------------------------------------------------------------------------------
    bool Cohort::isMature(){
        return (MaturityTimeStep < std::numeric_limits<unsigned>::max());
    }
