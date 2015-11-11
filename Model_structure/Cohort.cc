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
        origin=&gcl;
        destination=origin;
        ID=nextCohortID;//added MB to track this object.
        nextCohortID++;
    }
    //----------------------------------------------------------------------------------------------

    Cohort::Cohort(Cohort& actingCohort, double juvenileBodyMass, double adultBodyMass, double initialBodyMass, double initialAbundance, unsigned birthTimeStep, long long& nextCohortID) {

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
        positionInList = 0;
        origin = actingCohort.origin;
        destination = origin;
        ID = nextCohortID; //added MB to track this object.
        nextCohortID++;

    }
    //----------------------------------------------------------------------------------------------
    bool Cohort::isMature(){
        return (MaturityTimeStep < std::numeric_limits<unsigned>::max());
    }
    //----------------------------------------------------------------------------------------------
    vector<Cohort> Cohort::newCohorts;
    map<string,map<string,double>>Cohort::Deltas;
    //----------------------------------------------------------------------------------------------
    void Cohort::zeroDeltas() {
        // Initialize delta abundance sorted list with appropriate processes

        Deltas["abundance"]["mortality"] = 0.0;

        // Initialize delta biomass sorted list with appropriate processes
        Deltas["biomass"]["metabolism"] = 0.0;
        Deltas["biomass"]["predation"] = 0.0;
        Deltas["biomass"]["herbivory"] = 0.0;
        Deltas["biomass"]["reproduction"] = 0.0;

        // Initialize delta reproductive biomass vector with appropriate processes

        Deltas["reproductivebiomass"]["reproduction"] = 0.0;

        // Initialize organic pool delta vector with appropriate processes
        Deltas["organicpool"]["herbivory"] = 0.0;
        Deltas["organicpool"]["predation"] = 0.0;
        Deltas["organicpool"]["mortality"] = 0.0;

        // Initialize respiratory CO2 pool delta vector with appropriate processes
        Deltas["respiratoryCO2pool"]["metabolism"] = 0.0;
    }