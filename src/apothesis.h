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
    unordered_map< MicroProcesses::Process*, set< SurfaceTiles::Site* > > m_processMap;

    /// The number of flags given by the user
    int m_iArgc;

    /// The flags given by the user
    char** m_vcArgv;

    // Set debug mode
    bool m_debugMode;

    /// number of species
    int m_nSpecies;

    /// Analyzes the process and returns its type: Adsorption, Desorption, Diffusion or Reaction
    string mf_analyzeProc(string);

    double m_dRTot;
    double m_dEndTime;
    double m_dProcTime;
    double m_dProcRate;
    double m_dt;
    double m_iRandom;
    double m_dSum;
    int m_iSiteNum;
    bool m_bReportCoverages;
    bool m_bHasGrowth;
};

#endif // KMC_H
