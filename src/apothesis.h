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

class Apothesis
{
public:
    Apothesis();
    virtual ~Apothesis();

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

    /// Update apothesis from a previous run (used in ALD and ALE)
    void update( Utils::Parameters* parameters, Lattice* lattice );

    /// Return normalized probabilities of each process
    vector<double> calculateProbabilities(vector<MicroProcesses::Process*>);

    inline void setDebugMode(bool ifDebug) { m_debugMode = ifDebug;}
    bool getDebugMode() { return m_debugMode; }

    /// Return number of species
    int getNumSpecies();

    /// Given a process of the form A + * -> A(s) or A(s) -> A + * or A(s) -> A(s) or A(s) + B(s) -> AB(s) returns the reactants with the "*" included.
    vector<string> getReactants( string process );

    /// Given a process of the form A + * -> A(s) or A(s) -> A + * or A(s) -> A(s) or A(s) + B(s) -> AB(s) returns the products with the "*" included.
    vector<string> getProducts( string process );

    /// Given a reactant e.g. 2A it returns the 2 as stoichiometric coefficient and the A as the reactant
    pair<string, double> analyzeCompound( string reactant );

    IO *getIO() const;
    void setIO(IO *newPIO);

    void setLattice(Lattice* lattice);
    inline Lattice* getLattice() { return pLattice; } const;

    void buildLattice();

    void buildMicroProcesses();

    double getStartTime() const;
    void setStartTime(double newDStartTime);

    double getEndTime() const;
    void setEndTime(double newDEndTime);

private:
    IO* pIO;

    /// The process map which holds all the processes and the sites that each can be performed.
    map< MicroProcesses::Process*, set< SurfaceTiles::Site* > > m_processMap;

    // Set debug mode
    bool m_debugMode;

    /// number of species
    int m_nSpecies;

    /// Analyzes the process and returns its type: Adsorption, Desorption, Diffusion or Reaction
    string mf_analyzeProc(string);

    double m_dRTot;
    double m_dStartTime;
    double m_dEndTime;
    double m_dProcTime;
    double m_dProcRate;
    double m_dt;
    double m_iRandom;
    double m_dSum;
    int m_iSiteNum;
    bool m_bReportCoverages;
    bool m_bHasGrowth;
    bool m_bHasEtching;

    void mf_createWorkingDir( const string &name );
};

#endif // KMC_H
