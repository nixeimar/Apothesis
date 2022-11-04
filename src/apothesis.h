    //============================================================================
//    Apothesis: A kinetic Monte Calro (KMC) code for deposition processes.
//    Copyright (C) 2019  Nikolaos (Nikos) Cheimarios
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.

//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.
//============================================================================

#ifndef KMC_H
#define KMC_H

#include <map>
#include <list>
#include <vector>
#include <string>
#include <functional>
#include <set>
#include <valarray>
#include "species_new.h"
#include "species/species.h"

#define EXIT { printf("Apothesis terminated. \n"); exit( EXIT_FAILURE ); }

using namespace std;

/** The basic class of the kinetic monte carlo code. */

namespace Utils{ class ErrorHandler; class Parameters; class Properties; }
namespace SurfaceTiles{ class Site; }
namespace MicroProcesses { class Process; class Adsorption; class Desorption; class Diffusion; class SurfaceReaction; }
namespace RandomGen { class RandomGenerator; }

class Lattice;
class IO;
class Reader;

class Apothesis
{
public:
    Apothesis( int argc, char* argv[] );
    virtual ~Apothesis();

    /// Pointers to the classes that will share the common space i.e. the "pointer"
    /// Pointer to the input/output class
    IO* pIO;

    /// Ponter to the read class
    Reader* pReader;

    /// Pointer to the lattice class
    Lattice* pLattice;

    /// Pointer to the error class
    Utils::ErrorHandler* pErrorHandler;

    /// Pointer to the paramters class
    Utils::Parameters* pParameters;

    /// Pointer to the properties class
    Utils::Properties* pProperties;

    /// Random generator
    RandomGen::RandomGenerator *pRandomGen;

    /// Intialization of the KMC method. For example here the processes to be performed
    /// as these are written in the input file are constcucted through the factory method
    void init();

    /// Perform the KMC iteratios
    void exec();

    /// Function to log to output file whether a parameter is properly read
    void logSuccessfulRead(bool read, string parameter);

    /// Return access to list of species
    map<string, Species*> getAllSpecies();

    // Return species
    Species* getSpecies(string species);

    /// Return normalized probabilities of each process
    vector<double> calculateProbabilities(vector<MicroProcesses::Process*>);

    /// Return access to IO pointer
    inline IO* getIOPointer() { return pIO; }

    inline void setDebugMode(bool ifDebug) { m_debugMode = ifDebug;}
    bool getDebugMode() { return m_debugMode; }

    /// Return number of species
    int getNumSpecies();

private:
    /// The process map which holds all the processes and the sites that each can be performed.
    map< MicroProcesses::Process*, set< SurfaceTiles::Site* > > m_processMap;

    //This holds the process and the list of sites that can be performed.
    //This should replace m_processMap.
    //Have to check if string is the name of the process or should we put
    //set is log(n) in insert and delete and log(1) for delete
    // "Adsorption:CuAMD" must be a name of a process because we want to be able to do
    // m_procMap["Adsorption:CuAMD"]->addSite(Site);
    map< string, set< int > > m_procMap;

    //The map holding all tbe species. The int is the same as the ID of the species.
    // e.g Assuming the user has procided CuAMD, H2, SiH4, SiH2
    // then m_speciesMAP[ 0 ]= > CuAMD
    //      m_speciesMAP[ 1 ]= > H2
    //      m_speciesMAP[ 2 ]= > SiH4
    //      m_speciesMAP[ 3 ]= > SiH2
    map< int, species_new* > m_speciesMap;

    // Hold the name and the stichiometric coeficient of the reactants in the surface reactions
    // The elements of the vector can be accessed through thr id of the species.

    // e.g. 2CuAMD + 2H2 => ...
    // e.g. CuAMD + S => ...
    // CuAMD: 0
    // H2: 1
    // S: 2
    // e.g. Reaction_0: 2 2 0 0 0 0 0 0 etc.
    //      Reaction_1: 1 0 1 0 0 0 0 0 etc.
    // m_surfReactionsMap[ 0 ] -> will return the stoichiometric coeff for this reactions.
    map< string, valarray<int> > m_surfReacMap;

    // Same as m_surfReacMap holding the procudts. The string should be the same.
    map< string, valarray<int> > m_surfProdMap;

    /// The number of flags given by the user
    int m_iArgc;

    /// The flags given by the user
    char** m_vcArgv;

    // map of species
    map<string, Species*> m_species;

    // Set debug mode
    bool m_debugMode;

    /// number of species
    int m_nSpecies;

    //
    double m_dRTot;
    double m_dEndTime; // = 0.01;
    double m_dProcTime;
    double m_dProcRate; // = 0.0;
    double m_dt; // = 0.0;
    double m_dRandom; // = 0;
    double m_dSum; // = 0.0;
    int m_iSiteNum; // = 0;
    int n;
};

#endif // KMC_H
