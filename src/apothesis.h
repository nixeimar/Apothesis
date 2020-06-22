//============================================================================
//    Apothesis: A kinetic Monte Calro (KMC) code for deposotion processes.
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
#include "species.h"

#define EXIT { printf("Apothesis terminated. \n"); exit( EXIT_FAILURE ); }

using namespace std;

/** The basic class of the kinetic monte carlo code. */

namespace Utils{ class ErrorHandler; class Parameters;}
namespace SurfaceTiles{ class Site; }
namespace MicroProcesses { class Process; }
class Lattice;
class IO;
class Read;

class Apothesis
{
public:
    Apothesis( int argc, char* argv[] );
    virtual ~Apothesis();

    /// Pointers to the classes that will share the common space i.e. the "pointer"

    /// Ponter to the input/output class
    IO* pIO;

    /// Ponter to the read class
    Read* pRead;

    /// Pointer to the lattice class
    Lattice* pLattice;

    /// Pointer to the error class
    Utils::ErrorHandler* pErrorHandler;

    /// Pointer to the paramters class
    Utils::Parameters* pParameters;

    /// Intialization of the KMC method. For example here the processes to be performed
    /// as these are written in the input file are constcucted through the factory method
    void init();

    /// Perform the KMC iteratios
    void exec();

    /// Add a process
    void addProcess(string process);

    /// Return access to list of species
    vector<Species *> getSpecies();

private:
    /// The process map which holds all the processes and the sites that each can be performed.
    // Not to handy. Re-think... I have found another way... Implement it
    map< MicroProcesses::Process*, list< SurfaceTiles::Site* >* > m_processMap;

    /// Vector holding the processes to be performed.
    vector< MicroProcesses::Process* > m_vProcesses;

    /// Vector holding the name of the processes (string)
    vector<string> m_processes;

    /// Function to log to output file whether a parameter is properly read
    void logSuccessfulRead(bool read, string parameter);

    /// The number of flags given by the user
    int m_iArgc;

    /// The flags given by the user
    char** m_vcArgv;

    // List of species
    vector<Species *> m_species;
};

#endif // KMC_H
