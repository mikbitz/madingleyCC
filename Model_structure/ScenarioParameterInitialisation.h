#ifndef SCENARIOPARAMETERINITIALISATION_H
#define SCENARIOPARAMETERINITIALISATION_H
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <assert.h>

using namespace std;
/**
 * \file ScenarioParameterInitialisation.h
 * \brief The ScenarioParameterInitialisation header file
 */

/**
 * \class  ScenarioParameterInitialisation
 * \brief Reads the file specifying which scenarios will be run, and stores this information
 */
class ScenarioParameterInitialisation {
public:
    //----------------------------------------------------------------------------------------------
    //Variables
    //----------------------------------------------------------------------------------------------

    /** The number of scenarios to be run*/
    int scenarioNumber;
    /** \brief String parameters for the different scenarios to run */
    map<string, vector<string> > scenarioParameters;
    /** \brief The number of simulations to run for each of the scenarios */
    vector<int> scenarioSimulationsNumber;

    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------

    //----------------------------------------------------------------------------------------------
    /** \brief
     * Constructor for ScenarioParameterInitialisation: reads in scenario parameters from a specified file
     * 
     * @param scenarioParameterFile The name of the scenario parameters file, which must be in the 'Model setup' directory</param>
     * @param outputPath The directory to write output files to</param>
     */
    ScenarioParameterInitialisation(string scenarioParameterFile, string outputPath) {
        cout << "Reading scenario parameters file..." << endl;
        ifstream infile(scenarioParameterFile.c_str());
        if (infile.is_open()) {

            string l, header[2];
            getline(infile, l);
            //trim off newline character
            l.pop_back();
            istringstream s(l);
            for (unsigned i = 0; i < 2; i++) {
                getline(s, header[i], ',');
                transform(header[i].begin(), header[i].end(), header[i].begin(), ::tolower);
            }
            if (!(((header[0] == "number of simulations") && (header[1] == "human npp extraction")) ||
                    ((header[1] == "number of simulations") && (header[0] == "human npp extraction")))) {
                cout << "Bad header in scenario file " << scenarioParameterFile << endl;
                exit(1);
            }

            while (infile.good()) {
                string l, data;
                getline(infile, l);
                if (infile.good())l.pop_back();
                if (l.length() > 1) {
                    istringstream s(l);
                    for (unsigned i = 0; i < 2; i++) {
                        getline(s, data, ',');
                        if (header[i] == "human npp extraction")scenarioParameters["Human NPP Extraction"].push_back(data);
                        if (header[i] == "number of simulations")scenarioSimulationsNumber.push_back(atoi(data.c_str()));
                    }
                }
                scenarioNumber++;
            }

        } else {
            cout << "Something wrong with scenario parameter file " << scenarioParameterFile << endl;
        }
        infile.close();
    }
    //----------------------------------------------------------------------------------------------

};
#endif