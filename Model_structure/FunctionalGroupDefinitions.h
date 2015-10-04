#ifndef FUNCTIONALGROUPDEFINITIONS_H
#define FUNCTIONALGROUPDEFINITIONS_H
#include <string>
using namespace std;
/** \file FunctionalGroupDefinitions.h
 * \brief the FunctionalGroupDefinitions header file
 */

//    /// <summary>
//    /// Reads in and performs look-ups on functional group definitions
//    /// </summary>
//    /// <remarks>Mass bins values currently defined as middle of each mass bins</remarks>
//    /// <todoM>Throw error if there are any blanks in csv file</todoM>
class FunctionalGroupDefinitions
    {
        public:
/** \brief

//        //A lookup device: sorted dictionary keyed by Functional Trait and valued by a sorted dictionary itself keyed by Unique Functional Trait Values and 
//        //valued by an integer array of functional group indices corresponding to each functional trait value

//        /// A dictionary to allow functional group indices to be looked up based on trait values
*/
        map<string, map<string, vector<int>>> IndexLookupFromTrait;
//
/** \brief
//        /// A sorted list of all of the properties of functional groups and their values
*/
        map<string,vector<double> > FunctionalGroupProperties;
/** \brief
//        /// Dictionary to allow traits of functional groups to be looked up based on the functional group index
*/
        map<string, vector<string> > TraitLookupFromIndex;
//
/** \brief
//        /// A list of the indices of all functional groups in the model
*/
        vector<int> AllFunctionalGroupsIndex;

FunctionalGroupDefinitions(){;}
/** \brief
//        /// Constructor for the functional group definitions: reads in the specified functional group definition file, 
//        /// constructs lookup tables, mass ranges and initial cohort numbers in each functional group
@param fileName The name of the functional group definition file to be read in
@param outputPath The path to the output folder, in which to copy the functional group definitions file
*/
FunctionalGroupDefinitions(string fileName, string outputPath)
       {
        cout<<"Reading those functional group definitions"<<endl;
        fileName="input/Model setup/"+fileName;
        ifstream infile(fileName.c_str());
        if(infile.is_open()){
            
            string l;
            vector<string>header,category;
            getline(infile,l);
            //trim off newline character
            l.pop_back();
            istringstream s(l);
            //split out the comma-separated header
            while(s.good()){
                string tmp;
                getline(s,tmp,',');
                transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
                //split out the header category (definition.property or note)
                istringstream splt(tmp);
                string dp,op;
                getline(splt,dp,'_');
                category.push_back(dp);
                getline(splt,op,'_');
                header.push_back(op);
            }
            int count=0;
            //retrieve the lines defining each functional group
            while (infile.good()){
                AllFunctionalGroupsIndex.push_back(count);
                count++;
                string l,data;
                getline(infile,l);
                if (infile.good())l.pop_back();
                if (l.length()>1){
                    istringstream s(l);
                    //step through the columns for this functional group
                    for (unsigned i=0;i<header.size();i++){
                        getline(s,data,',');
                        transform(data.begin(), data.end(), data.begin(), ::tolower);

                        if(category[i]=="definition"){
                            //for each trait, store the value for a given functional group
                            //indexed by functional group number
                            TraitLookupFromIndex[header[i]].push_back(data);
                            //for a given trait, store the functional group number
                            //which has a given value for that trait
                            IndexLookupFromTrait[header[i]][data].push_back(i);
                        }
                        //Otherwise get the value for the given property
                        //for this functional group
                        if (category[i]=="property"){
                          FunctionalGroupProperties[header[i]].push_back(atof(data.c_str()));
                        }
                    }
                }
            }

            
        }else{
            cout<<"Something wrong with functional group definitions file "<<fileName<<endl;  
        }
        infile.close();

       }

/** \brief
//        /// Return the value of a biological parameter for a given parameter and functional group

@param propertyName The name of the biological parameter
@param functionalGroup Functional group index
*/
//        /// <returns>The value of the biological parameter for the specified functional group</returns>
double GetBiologicalPropertyOneFunctionalGroup(string propertyName, int functionalGroup)
        {
          transform(propertyName.begin(), propertyName.end(), propertyName.begin(), ::tolower);

            return FunctionalGroupProperties[propertyName][functionalGroup];
       }
//
/** \brief
//        /// Get values of a functional group property for all functional groups

@param propertyName The name of the property to get values for
*/
//        /// <returns>The values of a functional group property for all functional groups</returns>
vector<double> GetBiologicalPropertyAllFunctionalGroups(string propertyName)
       {
          transform(propertyName.begin(), propertyName.end(), propertyName.begin(), ::tolower);
          return FunctionalGroupProperties[propertyName];
       }

/** \brief
Retrieves the values for all traits defined in the model

@returns String array of traits defined for the model
*/
       vector<string> GetTraits()
       {
           vector<string> Traits ;

           for(auto var : TraitLookupFromIndex)
           {
               Traits.push_back(var.first);
           }

           return Traits;
       }

/** \brief
Retrieves the trait values for all traits defined in the model

@param Trait The trait for which trait values are to be found

@returns String array of trait values for the specifiec trait
*/
vector<string> GetUniqueTraitValues(string Trait)
       {
           vector<string> TraitValues;
                     transform(Trait.begin(), Trait.end(), Trait.begin(), ::tolower);

           map<string, vector<int> > temp = IndexLookupFromTrait[Trait];

           for (auto var : temp)
           {
               //if(!TraitValues.Contains(var.first)) TraitValues.push_back(var.first);
           }

           return TraitValues;
       }



/** \brief
//        /// Returns a string of Trait Names associated with the specified search trait and functional group index value

@param searchTrait The name of the trait to get values for
@param functionalGroupIndex The functional group index to return the trait value for
*/
//        /// <returns>The value of the specified trait for the specified functional group</returns>
string GetTraitNames(string searchTrait, int functionalGroupIndex)
        {
             transform(searchTrait.begin(), searchTrait.end(), searchTrait.begin(), ::tolower);

             return TraitLookupFromIndex[searchTrait][functionalGroupIndex];   
        }
//
//
/** \brief
//        /// Get the values of a set of specified traits for a specified functional group

@param searchTraits A vector of trait names to get values for
@param functionalGroupIndex The functional group index to return trait values for
*/
//        /// <returns>A vector of values of the specified traits for a specified functional group</returns>
vector<string> GetTraitNames(vector<string> searchTraits, int functionalGroupIndex)
        {
//            string[] TraitNames = new string[searchTraits.Length];
//
//            for (int nn = 0; nn < searchTraits.Length; nn++)
//            {
//                TraitNames[nn] = TraitLookupFromIndex[searchTraits[nn].ToLower()].GetValue(functionalGroupIndex).ToString();
//            }
//
//            return TraitNames;
        }
//
//
//
/** \brief
//        /// Get the functional group indices that have specified values of specified traits

@param searchTraits Vector of trait names to search for
@param searchTraitValues Vector of trait values to search for
@param intersection Whether the intersection of the indices for the traits should be returned, otherwise return the union of the indices

@return A vector of functional group indices with the specified values of the specified traits
*/
vector<int> GetFunctionalGroupIndex(vector<string> searchTraits, vector<string> searchTraitValues, bool intersection)
        {
//            // Check that the numbers of traits and of trait values specified are equal
//            Debug.Assert((searchTraits.Length == searchTraitValues.Length), "Unequal search string arrays");
//
//            // List to hold the functional group indices for each trait-trait value pair
//            List<int[]> IndexList = new List<int[]>();
//            
//            int[] TempIndexList;
//            
//            //Sorted dictionary to hold the trait value index list sorted dictionary from the lookup table
//			SortedDictionary<string, int[]> TraitIndexList;
//
//            //Loop over the number of trait name and trait value pairs
//            for (int nn = 0; nn < searchTraits.Length; nn++)
//			{ 
//                //Check if the trait name is in the lookup table and if so pull out the <trait value, index vector> sorted dictionary for it
//                if (IndexLookupFromTrait.TryGetValue(searchTraits[nn].ToLower(), out TraitIndexList))
//                {
//                    //Check if the trait value string is found in the lookup table and if found pull out the index vector for it
//                    //and add it to the List of these for processing - intersection of union
//                    if (TraitIndexList.TryGetValue(searchTraitValues[nn], out TempIndexList))
//                    {
//                        IndexList.Add(TempIndexList);
//                    }
//                    //If trait value string not found then show error message
//                    else
//                    {
//                        Debug.Fail("Trait Value to search for not found in lookup tables");
//                    }
//                }
//                //If trait name string not found then show error message
//                else
//                {
//                    Debug.Fail("Trait to search for not found in lookup tables");
//                }
//			}
//
//            //If we are only searching for one traitname and trait value pair then return the index vector
//            if (searchTraits.Length == 1)
//            {
//                return IndexList[0];
//            }
//            //Otherwise process the List of index vectors
//            else
//            {
//                //Object to hold the array of index values found by the intersection method
//                IEnumerable<int> ReturnList;
//                //If intersection true then find the index values common to all traitname and trait value pairs
//                if (intersection)
//                {
//                    ReturnList = IndexList[0].Intersect(IndexList[1]);
//
//                    for (int nn = 2; nn < IndexList.Count; nn++)
//                    {
//                        ReturnList = ReturnList.Intersect(IndexList[nn]);
//                    }
//                }
//                //If intersection false then return all the index values found above
//                else
//                {
//                    ReturnList = IndexList[0].Union(IndexList[1]);
//
//                    for (int nn = 2; nn < IndexList.Count; nn++)
//                    {
//                        ReturnList = ReturnList.Union(IndexList[nn]);
//                    }
//
//                }
//                return ReturnList.ToArray();
//            }
        }
//
/** \brief
//        /// Function to return the integer index values for functional groups corresponding to given trait and trait value pair combinations.
//        /// Overloaded to accept a single string rather than an array in the traits to search and the trait values - both must be single strings

@param searchTraits String of Trait names to search for trait values within
@param searchTraitValues String of string Trait Values to find functional group indices for
@param intersection Boolean statement indicating if you want the intersection of the indices. Only valid if more than one Trait and Trait Value pair.
//        /// True means give intersection. False means give the union of indices
@return Int array containing functional group indices corresponding to the given search conditions
*/
vector<int> GetFunctionalGroupIndex(string searchTraits, string searchTraitValues, bool intersection)
       {
 
            //List to hold the index vectors for each trait trait value pair
            vector<int> IndexList;
            transform(searchTraits.begin(), searchTraits.end(), searchTraits.begin(), ::tolower);
            transform(searchTraitValues.begin(), searchTraitValues.end(), searchTraitValues.begin(), ::tolower);
            //Sorted dictionary to hold the trait value index list sorted dictionary from the lookup table
            map<string, vector<int>> TraitIndexList;

            //Check if the trait name is in the lookup table and if so pull out the <trait value, index vector> sorted dictionary for it
            if (IndexLookupFromTrait.count(searchTraits)!=0)
            {
                //Check if the trait value string is found in the lookup table and if found pull out the index vector for it
                //and add it to the List of these for processing - intersection of union
                if (IndexLookupFromTrait[searchTraits].count(searchTraitValues)!=0)
                {
                    return IndexLookupFromTrait[searchTraits][searchTraitValues];
                }
                //If trait value string not found then show error message
                else
                {
                    cout<<"Trait Value to search for not found in lookup tables"<<endl;
                    exit(1);
                }
            }
            //If trait name string not found then show error message
            else
            {
                cout<<"Trait to search for not found in lookup tables"<<endl;
                exit(1);
            }

            return IndexList;
        }
//
/** \brief Returns number of functional groups */
//        /// <returns>Number of functional groups</returns>
int GetNumberOfFunctionalGroups()
       {
            return (AllFunctionalGroupsIndex.size());
       }
//    }
};
#endif