#ifndef COHORTMERGE_H
#define COHORTMERGE_H
#include <Cohort.h>
#include <GridCellCohortHandler.h>
/** \file CohortMerge.h
 * \brief The CohortMerge header file
 */
//namespace Madingley
//{

/** \brief Merges cohorts with similar properties
 */
class CohortMerge {
public:
    //----------------------------------------------------------------------------------------------
    //Variables
    //----------------------------------------------------------------------------------------------
    /** \brief An instance of the simple random number generator class */
    std::default_random_engine RandomNumberGenerator;

    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------

    //----------------------------------------------------------------------------------------------
    /** \brief Set the seed for the random number generator
    @param DrawRandomly  
     */
    void SetRandom(bool DrawRandomly) {

        if (DrawRandomly) {
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            RandomNumberGenerator.seed(seed);
        } else {
            RandomNumberGenerator.seed(4000);
        }
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Calculate the distance between two cohorts in multi-dimensional trait space (body mass, adult mass, juvenile mass)
    @param Cohort1 The first cohort to calculate distance to 
    @param Cohort2 The cohort to compare to 
    @returns The relative distance in trait space
     */
    double CalculateDistance(Cohort Cohort1, Cohort Cohort2) {
        double AdultMassDistance = abs(Cohort1.AdultMass - Cohort2.AdultMass) / Cohort1.AdultMass;
        double JuvenileMassDistance = abs(Cohort1.JuvenileMass - Cohort2.JuvenileMass) / Cohort1.JuvenileMass;
        double CurrentMassDistance = abs(Cohort1.IndividualBodyMass - Cohort2.IndividualBodyMass) / Cohort1.IndividualBodyMass;

        return sqrt((AdultMassDistance * AdultMassDistance) + (JuvenileMassDistance * JuvenileMassDistance) +
                (CurrentMassDistance * CurrentMassDistance));
    }
    //----------------------------------------------------------------------------------------------
    /** \brief  Merge cohorts until below a specified threshold number of cohorts in each grid cell
    @param gridCellCohorts The cohorts within this grid cell 
    @param TotalNumberOfCohorts The total number of cohorts in this grid cell 
    @param TargetCohortThreshold The target threshold to reduce the number of cohorts to 
    @return The number of cohorts that have been merged
     */
    struct tupleComparator {
        bool operator()(tuple<double, int, vector<int>> i, tuple<double, int, vector<int>> j) {
            double u, v;
            u = get<0>(i);
            v = get<0>(j);
            return (u < v);
        }
    } tupleCompare;
    //----------------------------------------------------------------------------------------------
    int MergeToReachThreshold(GridCellCohortHandler gridCellCohorts, int TotalNumberOfCohorts, int TargetCohortThreshold) {
        // A list of shortest distances between pairs of cohorts
        vector<tuple<double, int, vector<int>>> ShortestDistances;
        //A holding list
        vector<tuple<double, int, vector<int>>> Holdingvector;
        //Vector of lists of shortest distances in each functional group
        vector<vector < tuple<double, int, vector<int>>>> ShortestDistancesPerFunctionalGroup(gridCellCohorts.size());
        // Temporary
        vector<vector < tuple<double, int, vector<int>>>> ShortestDistancesPerFunctionalGroup2(gridCellCohorts.size());
        // How many cohorts to remove to hit the threshold
        int NumberToRemove = TotalNumberOfCohorts - TargetCohortThreshold;

        // Holds the pairwise distances between two cohorts; the functional group of the cohorts; the cohort IDs of each cohort
        tuple<double, int, vector<int>> PairwiseDistance;

        // Loop through functional groups
        for (int ff = 0; ff < gridCellCohorts.size(); ff++) {

            // Loop through cohorts within functional groups
            for (int cc = 0; cc < gridCellCohorts[ff].size() - 1; cc++) {

                // Loop through comparison cohorts
                for (int dd = cc + 1; dd < gridCellCohorts[ff].size(); dd++) {
                    vector<int> v;
                    // Randomly select which cohort is to be merged to & calculate distance between cohort pair

                    std::uniform_real_distribution<double> randomNumber(0.0, 1.0);
                    double RandomValue = randomNumber(RandomNumberGenerator);
                    if (RandomValue < 0.5) {
                        v = {cc, dd};
                    } else {
                        v = {dd, cc};
                    }
                    tuple<double, int, vector<int>>PairwiseDistance(CalculateDistance(gridCellCohorts[ff][cc], gridCellCohorts[ff][dd]), ff, v);

                    // Temporary
                    ShortestDistancesPerFunctionalGroup2[ff].push_back(PairwiseDistance);
                }


            }

            // Temporary
            sort(ShortestDistancesPerFunctionalGroup2[ff].begin(), ShortestDistancesPerFunctionalGroup2[ff].end(), tupleCompare);
        }


        // Hold the current position in the shortest distance list
        int CurrentvectorPosition = 0;

        // Now that the shortest distances have been calculated, do the merging execution                
        int FunctionalGroup;
        int CohortToMergeFrom;
        int CohortToMergeTo;

        // Temporary
        for (int gg = 0; gg < gridCellCohorts.size(); gg++) {
            CurrentvectorPosition = 0;
            while (CurrentvectorPosition < ShortestDistancesPerFunctionalGroup2[gg].size()) {
                CohortToMergeFrom = get<2>(ShortestDistancesPerFunctionalGroup2[gg][CurrentvectorPosition])[1];
                CohortToMergeTo = get<2>(ShortestDistancesPerFunctionalGroup2[gg][CurrentvectorPosition])[0];

                for (int cc = ShortestDistancesPerFunctionalGroup2[gg].size() - 1; cc > CurrentvectorPosition; cc--) {
                    if (get<2>(ShortestDistancesPerFunctionalGroup2[gg][cc])[0] == CohortToMergeFrom ||
                            get<2>(ShortestDistancesPerFunctionalGroup2[gg][cc])[1] == CohortToMergeFrom) {
                        ShortestDistancesPerFunctionalGroup2[gg].erase(ShortestDistancesPerFunctionalGroup2[gg].begin() + cc);
                    }

                }

                CurrentvectorPosition++;
            }
        }

        // Compile all shortest distances into a single list for merging purposes - note that we only need to do a limited number of merges
        for (int gg = 0; gg < gridCellCohorts.size(); gg++) {
            for (auto distance : ShortestDistancesPerFunctionalGroup2[gg]) {
                ShortestDistances.push_back(distance);
            }
        }

        sort(ShortestDistances.begin(), ShortestDistances.end(), tupleCompare);
        // Counts the number of merges that have happened
        int MergeCounter = 0;
        CurrentvectorPosition = 0;

        // While merging does not reach threshold, and while there are still elements in the list
        while ((MergeCounter < NumberToRemove) && (CurrentvectorPosition < ShortestDistances.size())) {
            // Get pairwise traits
            FunctionalGroup = get<1>(ShortestDistances[CurrentvectorPosition]);
            CohortToMergeFrom = get<2>(ShortestDistances[CurrentvectorPosition])[1];
            CohortToMergeTo = get<2>(ShortestDistances[CurrentvectorPosition])[0];

            // Check whether either cohort has already merged to something else this timestep merge
            //  execution, and hence is empty
            //                   if ((gridCellCohorts[FunctionalGroup][CohortToMergeFrom].CohortAbundance > 0) ||
            //                       (gridCellCohorts[FunctionalGroup][CohortToMergeTo].CohortAbundance > 0))
            //                   {

            // Add the abundance of the second cohort to that of the first
            gridCellCohorts[FunctionalGroup][CohortToMergeTo].CohortAbundance += (gridCellCohorts[FunctionalGroup][CohortToMergeFrom].CohortAbundance * gridCellCohorts[FunctionalGroup][CohortToMergeFrom].IndividualBodyMass) / gridCellCohorts[FunctionalGroup][CohortToMergeTo].IndividualBodyMass;
            // Add the reproductive potential mass of the second cohort to that of the first
            gridCellCohorts[FunctionalGroup][CohortToMergeTo].IndividualReproductivePotentialMass += (gridCellCohorts[FunctionalGroup][CohortToMergeFrom].IndividualReproductivePotentialMass * gridCellCohorts[FunctionalGroup][CohortToMergeFrom].CohortAbundance) / gridCellCohorts[FunctionalGroup][CohortToMergeTo].CohortAbundance;
            // Set the abundance of the second cohort to zero
            gridCellCohorts[FunctionalGroup][CohortToMergeFrom].CohortAbundance = 0.0;
            // Designate both cohorts as having merged
            gridCellCohorts[FunctionalGroup][CohortToMergeTo].Merged = true;
            gridCellCohorts[FunctionalGroup][CohortToMergeFrom].Merged = true;

            MergeCounter++;
            CurrentvectorPosition++;

            //                   }
            //                   else
            //                   {
            //                       CurrentvectorPosition++;
            //                   }
        }

        return MergeCounter;

    }
    //----------------------------------------------------------------------------------------------
    /** \brief
    Merge cohorts until below a specified threshold number of cohorts in each grid cell

    @param gridCellCohorts The cohorts within this grid cell 
    @param TotalNumberOfCohorts The total number of cohorts in this grid cell 
    @param TargetCohortThreshold The target threshold to reduce the number of cohorts to 
    @return The number of cohorts that have been merged
     */
    int MergeToReachThresholdFast(GridCellCohortHandler& gridCellCohorts, int TotalNumberOfCohorts, int TargetCohortThreshold) {

        // A list of shortest distances between pairs of cohorts
        vector<tuple<double, int, vector<int>>> ShortestDistances;

        // Vector of lists of shortest distances in each functional group
        vector<vector < tuple<double, int, vector<int>>>> ShortestDistancesPerFunctionalGroup(gridCellCohorts.size());
        // How many cohorts to remove to hit the threshold
        int NumberToRemove = TotalNumberOfCohorts - TargetCohortThreshold;

        // Holds the pairwise distances between two cohorts; the functional group of the cohorts; the cohort IDs of each cohort
        tuple<double, int, vector<int>> PairwiseDistance;

        //Loop through functional groups
        for (int ff = 0; ff < gridCellCohorts.size(); ff++) {
            // Loop through cohorts within functional groups
            for (int cc = 0; cc < gridCellCohorts[ff].size() - 1; cc++) {
                // A holding list
                vector<tuple<double, int, vector<int>>> Holdingvector;
                // Loop through comparison cohorts
                for (int dd = cc + 1; dd < gridCellCohorts[ff].size(); dd++) {
                    vector<int> v;
                    // Randomly select which cohort is to be merge to & calculate distance between cohort pair
                    std::uniform_real_distribution<double> randomNumber(0.0, 1.0);
                    double RandomValue = randomNumber(RandomNumberGenerator);
                    if (RandomValue < 0.5) {
                        v = {cc, dd};
                    } else {
                        v = {dd, cc};
                    }
                    // Randomly select which cohort is to be merge to & calculate distance between cohort pair

                    tuple<double, int, vector<int>>PairwiseDistance(CalculateDistance(gridCellCohorts[ff][cc], gridCellCohorts[ff][dd]), ff, v);

                    Holdingvector.push_back(PairwiseDistance);

                }

                sort(Holdingvector.begin(), Holdingvector.end(), tupleCompare);

                // Sort through and only keep those cohorts which are necessary

                // The value to which to compare
                int ValueToCompareTo = cc;
                if (Holdingvector.size()>0)ShortestDistancesPerFunctionalGroup[ff].push_back(Holdingvector[0]);

                // Only add to main list those that are valid. Note that this doesn't catch everything (because we don't yet know the full ordering), 
                // but clears out a lot of redundant information from the list
                int position = 0;

                while (position < Holdingvector.size()) {
                    if (get<2>(Holdingvector[position])[1] == ValueToCompareTo) {
                        ShortestDistancesPerFunctionalGroup[ff].push_back(Holdingvector[position]);
                        break;
                    } else
                        ShortestDistancesPerFunctionalGroup[ff].push_back(Holdingvector[position]);

                    position++;
                }


            }
            sort(ShortestDistancesPerFunctionalGroup[ff].begin(), ShortestDistancesPerFunctionalGroup[ff].end(), tupleCompare);

        }


        //Hold the current position in the shortest distance list
        int CurrentvectorPosition = 0;

        int FunctionalGroup;
        int CohortToMergeFrom;
        int CohortToMergeTo;

        for (int ff = 0; ff < gridCellCohorts.size(); ff++) {
            CurrentvectorPosition = 0;
            while (CurrentvectorPosition < ShortestDistancesPerFunctionalGroup[ff].size()) {
                CohortToMergeFrom = get<2>(ShortestDistancesPerFunctionalGroup[ff][CurrentvectorPosition])[1];
                CohortToMergeTo = get<2>(ShortestDistancesPerFunctionalGroup[ff][CurrentvectorPosition])[0];

                for (int cc = ShortestDistancesPerFunctionalGroup[ff].size() - 1; cc > CurrentvectorPosition; cc--) {
                    if (get<2>(ShortestDistancesPerFunctionalGroup[ff][cc])[0] == CohortToMergeFrom ||
                            get<2>(ShortestDistancesPerFunctionalGroup[ff][cc])[1] == CohortToMergeFrom) {
                        ShortestDistancesPerFunctionalGroup[ff].erase(ShortestDistancesPerFunctionalGroup[ff].begin() + cc);
                        ;
                    }

                }

                CurrentvectorPosition++;
            }
        }



        // Compile all shortest distances into a single list for merging purposes - note that we only need to do a limited number of merges
        for (int ff = 0; ff < gridCellCohorts.size(); ff++) {
            for (auto distance : ShortestDistancesPerFunctionalGroup[ff]) {
                ShortestDistances.push_back(distance);
            }
        }

        sort(ShortestDistances.begin(), ShortestDistances.end(), tupleCompare);


        //Counts the number of merges that have happened
        int MergeCounter = 0;
        CurrentvectorPosition = 0;

        // While merging does not reach threshold, and while there are still elements in the list
        while ((MergeCounter < NumberToRemove) && (CurrentvectorPosition < ShortestDistances.size())) {
            // Get pairwise traits
            FunctionalGroup = get<1>(ShortestDistances[CurrentvectorPosition]);
            CohortToMergeFrom = get<2>(ShortestDistances[CurrentvectorPosition])[1];
            CohortToMergeTo = get<2>(ShortestDistances[CurrentvectorPosition])[0];

            // Check whether either cohort has already merged to something else this timestep merge
            // execution, and hence is empty
            //          if ((gridCellCohorts[FunctionalGroup][CohortToMergeFrom].CohortAbundance > 0) ||
            //              (gridCellCohorts[FunctionalGroup][CohortToMergeTo].CohortAbundance > 0))
            //          {

            // Add the abundance of the second cohort to that of the first
            gridCellCohorts[FunctionalGroup][CohortToMergeTo].CohortAbundance += (gridCellCohorts[FunctionalGroup][CohortToMergeFrom].CohortAbundance * gridCellCohorts[FunctionalGroup][CohortToMergeFrom].IndividualBodyMass) / gridCellCohorts[FunctionalGroup][CohortToMergeTo].IndividualBodyMass;
            // Add the reproductive potential mass of the second cohort to that of the first
            gridCellCohorts[FunctionalGroup][CohortToMergeTo].IndividualReproductivePotentialMass += (gridCellCohorts[FunctionalGroup][CohortToMergeFrom].IndividualReproductivePotentialMass * gridCellCohorts[FunctionalGroup][CohortToMergeFrom].CohortAbundance) / gridCellCohorts[FunctionalGroup][CohortToMergeTo].CohortAbundance;
            // Set the abundance of the second cohort to zero
            gridCellCohorts[FunctionalGroup][CohortToMergeFrom].CohortAbundance = 0.0;
            // Designate both cohorts as having merged
            gridCellCohorts[FunctionalGroup][CohortToMergeTo].Merged = true;
            gridCellCohorts[FunctionalGroup][CohortToMergeFrom].Merged = true;

            MergeCounter++;
            CurrentvectorPosition++;

            //                //       }
            //                //   else
            //                //   {
            //                //        CurrentvectorPosition++;
            //                //      }
        }

        return MergeCounter;

    }
    //----------------------------------------------------------------------------------------------
    /** \brief  Merge cohorts for responsive dispersal only; merges identical cohorts, no matter how many times they have been merged before

    @param gridCellCohorts The cohorts in the current grid cell 
    @return Number of cohorts merged
     */
    int MergeForResponsiveDispersalOnly(GridCellCohortHandler gridCellCohorts) {
        // Variable to track the total number of cohorts merged
        int NumberCombined = 0;

        //Loop over all functional groups
        for (int i = 0; i < gridCellCohorts.size(); i++) {
            // Loop over each cohort in each functional group
            for (int j = 0; j < gridCellCohorts[i].size(); j++) {
                // If that cohort has abundance greater than zero  then check if there are similar cohorts that could be merged with it
                if (gridCellCohorts[i][j].CohortAbundance > 0) {
                    // Loop over all cohorts above the jth in the cohort list
                    for (int k = j + 1; k < gridCellCohorts[i].size(); k++) {
                        // Check that kth cohort has abunance and that the two cohorts being compared do not represent a juvenile adult pairing
                        if (gridCellCohorts[i][k].CohortAbundance > 0 &&
                                ((gridCellCohorts[i][j].MaturityTimeStep == std::numeric_limits<unsigned>::max() && gridCellCohorts[i][k].MaturityTimeStep == std::numeric_limits<unsigned>::max()) ||
                                (gridCellCohorts[i][j].MaturityTimeStep < std::numeric_limits<unsigned>::max() && gridCellCohorts[i][k].MaturityTimeStep < std::numeric_limits<unsigned>::max()))) {
                            //Check that the individual masses are widentical
                            if (gridCellCohorts[i][j].IndividualBodyMass == gridCellCohorts[i][k].IndividualBodyMass) {
                                //Check that the adult masses are similar
                                if (gridCellCohorts[i][j].AdultMass == gridCellCohorts[i][k].AdultMass) {
                                    //Check that the juvenile masses are similar
                                    if (gridCellCohorts[i][j].JuvenileMass == gridCellCohorts[i][k].JuvenileMass) {
                                        //Check that the Maximum achieved mass is similar
                                        if (gridCellCohorts[i][j].MaximumAchievedBodyMass == gridCellCohorts[i][k].MaximumAchievedBodyMass) {
                                            // In half of cases, add the abundance of the second cohort to that of the first and maintain the properties of the first
                                            std::uniform_real_distribution<double> randomNumber(0.0, 1.0);
                                            double RandomValue = randomNumber(RandomNumberGenerator);
                                            if (RandomValue < 0.5) {
                                                // Add the abundance of the second cohort to that of the first
                                                gridCellCohorts[i][j].CohortAbundance += (gridCellCohorts[i][k].CohortAbundance * gridCellCohorts[i][k].IndividualBodyMass) / gridCellCohorts[i][j].IndividualBodyMass;
                                                // Set the abundance of the second cohort to zero
                                                gridCellCohorts[i][k].CohortAbundance = 0.0;
                                                // Add the reproductive potential mass of the second cohort to that of the first
                                                gridCellCohorts[i][j].IndividualReproductivePotentialMass += (gridCellCohorts[i][k].IndividualReproductivePotentialMass * gridCellCohorts[i][k].CohortAbundance) / gridCellCohorts[i][j].CohortAbundance;
                                                // Designate both cohorts as having merged
                                                gridCellCohorts[i][j].Merged = true;
                                                gridCellCohorts[i][k].Merged = true;
                                            }                                                // In all other cases, add the abundance of the first cohort to that of the second and maintain the properties of the second
                                            else {
                                                // Add the abundance of the first cohort to that of the second
                                                gridCellCohorts[i][k].CohortAbundance += (gridCellCohorts[i][j].CohortAbundance * gridCellCohorts[i][j].IndividualBodyMass) / gridCellCohorts[i][k].IndividualBodyMass;
                                                // Set the abundance of the second cohort to zero
                                                gridCellCohorts[i][j].CohortAbundance = 0.0;
                                                // Add the reproductive potential mass of the second cohort to that of the first
                                                gridCellCohorts[i][k].IndividualReproductivePotentialMass += (gridCellCohorts[i][j].IndividualReproductivePotentialMass * gridCellCohorts[i][j].CohortAbundance) / gridCellCohorts[i][k].CohortAbundance;
                                                // Designate both cohorts as having merged
                                                gridCellCohorts[i][j].Merged = true;
                                                gridCellCohorts[i][k].Merged = true;
                                            }
                                            // Increment the number of cohorts combined
                                            NumberCombined += 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return NumberCombined;

    }
    //----------------------------------------------------------------------------------------------
};

#endif
