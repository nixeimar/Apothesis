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

#ifndef READ_H
#define READ_H

#include <fstream>


#include <string>
#include <vector>
#include <map>
#include <algorithm>


#include "pointers.h"
#include "apothesis.h"
#include "lattice.h"
#include "FCC.h"
#include "species.h"

#include "errorhandler.h"
#include "parameters.h"

// Include RapidJSON header files
#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"
#include "rapidjson/pointer.h"
#include "rapidjson/rapidjson.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/writer.h"
#include "rapidjson/filewritestream.h"
#include <iostream>


#if defined( _WIN32) || defined( _WIN64)
#include <direct.h>
#include <io.h>
#define GetCurrentDir _getcwd
#define BAC '\\'
#elif ( __linux__)
#include <unistd.h>
#include <dirent.h>
#define GetCurrentDir getcwd
#define BAC '/'
#endif

using namespace std;
using namespace Utils;

/** Tha class for handling input/output operations */

class Pointers;

typedef rapidjson::Document Document;

class Read: public Pointers
  {
  public:
    enum CASE{ Sensitive, Insensitive };

    Read();
    Read( Apothesis *apothesis);

    virtual ~Read();

    /// Initialization of the files: Assigns paths and file names
    /// If they are not given by the user, the default values are used
    /// Dafault: The path of the executable
    /// Input file name: input.kmc
    /// Output file name: output.log

    /// Reads the input file " .kmc".
    Document readInputFile(string filename);

    Document& getDoc();

    //Document getInput() const;

    // Return vector of species names
    vector<string> getSpeciesNames();

    // Return vector of molecular weights
    vector<double> getMWs();

  protected:

    // Parsed input file
    Document m_input;

    /// The type of lattice
    string m_sLatticeType;

    /// Supported lattice types
    map< string, Lattice::Type> m_LatticeType;

    /// The processes to be constructed with their parameters (per process)
    map< string, list< string > >  m_mProc;

    /// The input file
    ifstream m_InputFile;

    /// The output file
    ofstream m_OutFile;

    /// Keywords:
    /// Process keyword
    string m_sProcess;

    /// Lattice keyword
    string m_sLattice;

    /// Temperature keyword
    string m_sTemperature;

    /// Pressure keyword
    string m_sPressure;

    /// Comment
    string m_sCommentLine;

    ///  Num of iterations keyword.
    string m_sIterations;

    // Dimensions in the lattice
    int m_dimensions;

    // Name of species
    vector<string> m_speciesName;

    // List of molecular weights
    vector<double> m_MWs;

  };

#endif // READ_H
